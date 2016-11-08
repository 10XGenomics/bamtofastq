extern crate docopt;
extern crate rust_htslib;
extern crate flate2;
extern crate shardio;
extern crate bincode;
extern crate itertools;
extern crate rustc_serialize;
extern crate regex;
extern crate tempfile;
extern crate tempdir;

use std::io::{Write, BufWriter, Result};
use std::fs::File;
use std::fs::create_dir;
use std::path::{PathBuf, Path};
use std::hash::{Hash, SipHasher, Hasher};

use std::collections::HashMap;
use itertools::Itertools;

use tempfile::NamedTempFile;


use flate2::write::GzEncoder;
use flate2::Compression;

use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::{Record};

use bincode::rustc_serialize::{encode_into, decode};
use shardio::shard::{Serializer, Shardable, ShardWriteManager, ShardReader};

use regex::Regex;

use docopt::Docopt;

const USAGE: &'static str = "
10x Genomics BAM to FASTQ converter.

Usage:
  bamtofastq [options] <bam> <output-path>
  bamtofastq (-h | --help)
  bamtofastq --version

Options:
  --sample=NAME        Sample name for FASTQ files [default=sample]
  --reads-per-fastq=N  Number of reads per FASTQ chunk [default=200000000]
  --gemcode            Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)
  --lr20               Convert a BAM produced by Longranger 2.0
  --cr11               Convert a BAM produced by Cell Ranger 1.0-1.1
  -h --help            Show this screen.
  --version            Show version.
";

#[derive(Debug, RustcDecodable)]
pub struct Args {
    arg_bam: String,
    arg_output_path: String,
    flag_reads_per_fastq: usize,
    flag_sample: String,
    flag_gemcode: bool,
    flag_lr20: bool,
    flag_cr11: bool,
}

/// A Fastq record ready to be written
#[derive(Debug, RustcEncodable, RustcDecodable, PartialEq, PartialOrd, Eq, Ord)]
struct FqRecord {
    head: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>
}

/// Which read in a pair we have
#[derive(Clone, Copy, Debug, RustcEncodable, RustcDecodable, PartialOrd, Ord, PartialEq, Eq)]
enum ReadNum {
    R1,
    R2
}

/// Internally serialized read. Used for re-uniting discordant read pairs
#[derive(Debug, RustcEncodable, RustcDecodable, PartialOrd, Ord, Eq, PartialEq)]
struct SerFq {
    rec: FqRecord,
    read_num: ReadNum,
    i1: Option<FqRecord>,
    i2: Option<FqRecord>,
}

/// Serialization of SerFq structs
#[derive(Clone)]
pub struct SerFqImpl {}

impl Serializer<SerFq> for SerFqImpl {
    fn serialize(&self, items: &Vec<SerFq>, buf: &mut Vec<u8>) {
        encode_into(items, buf, bincode::SizeLimit::Infinite).unwrap();
    }

    fn deserialize(&self, buf: &mut Vec<u8>, data: &mut Vec<SerFq>) {
        let mut buf_slice = buf.as_mut_slice();
        let r: Vec<SerFq> = decode(&mut buf_slice).unwrap();
        data.extend(r);
    }
}

/// Shard SerFq structs based on a hash of the head string
impl Shardable for SerFq {
    fn shard(&self) -> usize {
        let mut s = SipHasher::new();
        self.rec.head.hash(&mut s);
        s.finish() as usize
    }
}

/// Entry in the conversion spec from a BAM record back to a read.
/// Each read can be composed of data from a pair of tags (tag w/ sequence, tag w/ qual),
/// or a fixed-length sequence of Ns (with a default QV), or the sequence in the read.
#[derive(Clone, Debug, PartialEq)]
enum SpecEntry {
    Tags(String, String),
    Ns(usize),
    Read,
}

/// Spec for converting from a BAM record back to reads. Empty vector indicates that this read doesn't exist
/// in the output. The i1 and i2 reads should be buildable from tags in the R1 read.
#[derive(Clone)]
struct FormatBamRecords {
    r1_spec: Vec<SpecEntry>,
    r2_spec: Vec<SpecEntry>,
    i1_spec: Vec<SpecEntry>,
    i2_spec: Vec<SpecEntry>,
}

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("unrecognized"),
    }
}

/// Class that can convert a BAM record into Fastq sequences, given some conversion specs
impl FormatBamRecords {

    /// Read the conversion spec from the special @CO 10x_bam_to_fastq tags in the BAM header
    pub fn from_headers(reader: &bam::Reader) -> Option<FormatBamRecords> {

        let mut spec = Self::parse_spec(reader);
        if spec.len() == 0 {
            None
        } else {
            Some(
                FormatBamRecords {
                    r1_spec: spec.remove("R1").unwrap(),
                    r2_spec: spec.remove("R2").unwrap(),
                    i1_spec: spec.remove("I1").unwrap_or_else(|| Vec::new()),
                    i2_spec: spec.remove("I2").unwrap_or_else(|| Vec::new()),
            })
        }
    }

    /// hard-coded spec for gemcode BAM files
    pub fn gemcode() -> FormatBamRecords {

        FormatBamRecords {
            r1_spec: vec![SpecEntry::Read],
            r2_spec: vec![SpecEntry::Read],
            i1_spec: vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())],
            i2_spec: vec![SpecEntry::Tags("RX".to_string(), "QX".to_string())],
        }
    }

    // hard-coded specs for longranger 2.0 BAM files
    pub fn lr20() -> FormatBamRecords {

        FormatBamRecords {
            r1_spec: vec![SpecEntry::Tags("RX".to_string(), "QX".to_string()), SpecEntry::Ns(7), SpecEntry::Read],
            r2_spec: vec![SpecEntry::Read],
            i1_spec: vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())],
            i2_spec: vec![],
        }
    } 

    // hard-coded specs for cellranger 1.0-1.1 BAM files
    pub fn cr11() -> FormatBamRecords {

        FormatBamRecords {
            r1_spec: vec![SpecEntry::Read],
            r2_spec: vec![SpecEntry::Tags("UR".to_string(), "UQ".to_string())],
            i1_spec: vec![SpecEntry::Tags("CR".to_string(), "CQ".to_string())],
            i2_spec: vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())],
        }
    }

    /// Parse the specs from BAM headers if available
    fn parse_spec(reader: &bam::Reader) -> HashMap<String, Vec<SpecEntry>> {

        // Example header line:
        // @CO	10x_bam_to_fastq:R1(RX:QX,TR:TQ,SEQ:QUAL)
        let re = Regex::new(r"@CO\t10x_bam_to_fastq:(\S+)\((\S+)\)").unwrap();
        let text = String::from_utf8(reader.header.text()).unwrap();
        let mut spec = HashMap::new();

        for l in text.lines() {
            match re.captures(l) {
                Some(c) => {
                    let mut read_spec = Vec::new();

                    let read = c.at(1).unwrap().to_string();
                    let tag_list = c.at(2).unwrap();
                    let spec_elems = tag_list.split(',');
                    for el in spec_elems {
                        if el == "SEQ:QUAL" {
                            read_spec.push(SpecEntry::Read)
                        } else {
                            let mut parts = el.split(':');
                            let rtag = parts.next().unwrap().to_string();
                            let qtag = parts.next().unwrap().to_string();
                            read_spec.push(SpecEntry::Tags(rtag, qtag));
                        }
                    }

                    spec.insert(read, read_spec);
                }
                None => ()
            }
        }

        println!("spec: {:?}", spec);
        spec
    }

    /// Convert a BAM record to a Fq record, for internal caching
    pub fn bam_rec_to_ser(&self, rec: &Record) -> SerFq {
        match (rec.is_first_in_template(), rec.is_last_in_template()) {
            (true, false) => {
                SerFq {
                    read_num: ReadNum::R1,
                    rec: self.bam_rec_to_fq(rec, &self.r1_spec).unwrap(),
                    i1: if self.i1_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i1_spec).unwrap()) }  else { None },
                    i2: if self.i2_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i2_spec).unwrap()) }  else { None },
                    
                }
            },
            (false, true) => {
                SerFq {
                    read_num: ReadNum::R2,
                    rec: self.bam_rec_to_fq(rec, &self.r2_spec).unwrap(),
                    i1: if self.i1_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i1_spec).unwrap()) }  else { None },
                    i2: if self.i2_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i2_spec).unwrap()) }  else { None },
                }
            },
            _ => panic!("Not a valid read pair"),
        }
    }

    /// Convert a BAM record to Fq record ready to be written
    pub fn bam_rec_to_fq(&self, rec: &Record, spec: &Vec<SpecEntry>) -> Result<FqRecord> {

        let mut head = Vec::new();
        head.extend_from_slice(rec.qname());

        // Reconstitute read and QVs
        let mut r = Vec::new();
        let mut q = Vec::new();

        for (idx, item) in spec.into_iter().enumerate() {

            // It OK for the final tag in the spec to be missing from the read
            let last_item = idx == spec.len() - 1;

            match item {
                // Data from a tag
                &SpecEntry::Tags(ref read_tag, ref qv_tag) => {  

                    let rx = rec.aux(read_tag.as_bytes());
                    if rx.is_none() && last_item {
                        continue;
                    } else if rx.is_none() {
                        panic!(format!("read: {:?} missing: {:?}", rec.qname(), read_tag));
                    } else {
                        let rx = rx.unwrap().string();
                        r.extend_from_slice(rx);
                    }

                    let qx = rec.aux(qv_tag.as_bytes());
                    if qx.is_none() && last_item {
                        continue;
                    } else if qx.is_none() {
                        panic!(format!("read: {:?} missing: {:?}", rec.qname(), qv_tag));
                    } else {
                        let qx = qx.unwrap().string();
                        q.extend_from_slice(qx);   
                    }
                },

                // Just hardcode some Ns -- for cases where we didn't retain the required data
                &SpecEntry::Ns(len) => {
                    for _ in 0 .. len {
                        r.push(b'N');
                        q.push(b'J');
                    }
                }

                &SpecEntry::Read => {
                    // The underlying read
                    let mut seq = rec.seq().as_bytes();
                    let mut qual: Vec<u8> = rec.qual().iter().map(|x| x + 33).collect();

                    if rec.is_reverse() {
                        seq.reverse();
                        for b in seq.iter_mut() {
                            *b = complement(*b);
                        }

                        qual.reverse();
                    }

                    r.extend(seq);
                    q.extend(qual);
                }
            }
        }

        let fq_rec = FqRecord {
                head: head.clone(),
                seq: r,
                qual: q,
        };

        Ok(fq_rec)
    }

    pub fn format_read_pair(&self, r1_rec: &Record, r2_rec: &Record) -> Result<(FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)> {
        let r1 = self.bam_rec_to_fq(r1_rec, &self.r1_spec).unwrap();
        let r2 = self.bam_rec_to_fq(r2_rec, &self.r2_spec).unwrap();

        let i1 = if self.i1_spec.len() > 0 {
             Some(self.bam_rec_to_fq(r1_rec, &self.i1_spec).unwrap())
        } else {
            None
        };

        let i2 = if self.i2_spec.len() > 0 {
             Some(self.bam_rec_to_fq(r1_rec, &self.i2_spec).unwrap())
        } else {
            None
        };

        Ok((r1, r2, i1, i2))
    }


    pub fn format_read(&self, rec: &Record) -> Result<(FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)> {
        let r1 = self.bam_rec_to_fq(rec, &self.r1_spec).unwrap();
        let r2 = self.bam_rec_to_fq(rec, &self.r2_spec).unwrap();

        let i1 = if self.i1_spec.len() > 0 {
             Some(self.bam_rec_to_fq(rec, &self.i1_spec).unwrap())
        } else {
            None
        };

        let i2 = if self.i2_spec.len() > 0 {
             Some(self.bam_rec_to_fq(rec, &self.i2_spec).unwrap())
        } else {
            None
        };

        Ok((r1, r2, i1, i2))
    }

    /// A spec implies double-ended reads if both the R1 and R2 reads generate different BAM records.
    /// If not the R1 and R2 sequences can be derived from a single BAM entry.
    pub fn is_double_ended(&self) -> bool {
        self.r1_spec.contains(&SpecEntry::Read) && self.r2_spec.contains(&SpecEntry::Read)
    }
}

type BGW = BufWriter<GzEncoder<File>>;

/// Open Fastq files being written to
struct FastqWriter {
    formatter: FormatBamRecords,
    out_path: PathBuf,
    sample_name: String,

    r1: BGW,
    r2: BGW,
    i1: Option<BGW>,
    i2: Option<BGW>,

    chunk_written: usize,
    total_written: usize,
    n_chunks: usize,
    reads_per_fastq: usize,
}

/// Write sets of Fastq records to the open Fastq files
impl FastqWriter {

    pub fn new(out_path: &Path, formatter: FormatBamRecords, sample_name: String, reads_per_fastq: usize) -> FastqWriter {

        // open output files
        let (r1_path, r2_path, i1_path, i2_path) = Self::get_paths(&out_path, &sample_name, 0);
        let r1 = Self::open_gzip_writer(r1_path);
        let r2 = Self::open_gzip_writer(r2_path);
        let i1 = if formatter.i1_spec.len() > 0 { Some(Self::open_gzip_writer(i1_path)) } else { None };
        let i2 = if formatter.i2_spec.len() > 0 { Some(Self::open_gzip_writer(i2_path)) } else { None };

        FastqWriter {
            formatter: formatter,
            out_path: out_path.to_path_buf(),
            sample_name: sample_name,
            r1: r1,
            r2: r2,
            i1: i1,
            i2: i2,
            n_chunks: 1,
            total_written: 0,
            chunk_written: 0,
            reads_per_fastq: reads_per_fastq,
        }
    }

    fn get_paths(out_path: &Path, sample_name: &str, n_files: usize) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
        let r1 = out_path.join(format!("{}_S1_L001_R1_{:03}.fastq.gz", sample_name, n_files+1));
        let r2 = out_path.join(format!("{}_S1_L001_R2_{:03}.fastq.gz", sample_name, n_files+1));
        let i1 = out_path.join(format!("{}_S1_L001_I1_{:03}.fastq.gz", sample_name, n_files+1));
        let i2 = out_path.join(format!("{}_S1_L001_I2_{:03}.fastq.gz", sample_name, n_files+1));
        (r1, r2, i1, i2)
    }

    pub fn write_rec(w: &mut BGW, rec: &FqRecord)  {
        w.write(b"@").unwrap();
        w.write(&rec.head).unwrap();
        w.write(b"\n").unwrap();

        w.write(&rec.seq).unwrap();
        w.write(b"\n+\n").unwrap();
        w.write(&rec.qual).unwrap();
        w.write(b"\n").unwrap();
    }

    pub fn try_write_rec(w: &mut Option<BGW>, rec: &Option<FqRecord>) {
        match w {
            &mut Some(ref mut w) => match rec { &Some(ref r) => FastqWriter::write_rec(w, r), &None => panic!("setup error") },
            &mut None => ()
        }
    }

    /// Write a set of fastq records
    pub fn write(&mut self, r1: &FqRecord, r2: &FqRecord, i1: &Option<FqRecord>, i2: &Option<FqRecord>) {
        FastqWriter::write_rec(&mut self.r1, r1);
        FastqWriter::write_rec(&mut self.r2, r2);
        FastqWriter::try_write_rec(&mut self.i1, i1);
        FastqWriter::try_write_rec(&mut self.i2, i2);
        self.total_written += 1;
        self.chunk_written += 1;

        if self.chunk_written >= self.reads_per_fastq {
            self.cycle_writers()
        }
    }

    /// Open up a fresh output chunk
    fn cycle_writers(&mut self) {
        let (r1_path, r2_path, i1_path, i2_path) = Self::get_paths(&self.out_path, &self.sample_name, self.n_chunks);
        self.r1 = Self::open_gzip_writer(r1_path);
        self.r2 = Self::open_gzip_writer(r2_path);
        self.i1 = if self.formatter.i1_spec.len() > 0 { Some(Self::open_gzip_writer(i1_path)) } else { None };
        self.i2 = if self.formatter.i2_spec.len() > 0 { Some(Self::open_gzip_writer(i2_path)) } else { None };
        self.n_chunks += 1;
        self.chunk_written = 0;
    }

    fn open_gzip_writer<P: AsRef<Path>>(path: P) -> BufWriter<GzEncoder<File>> {
        let f = File::create(path).unwrap();
        let gz = GzEncoder::new(f, Compression::Fast);
        BufWriter::new(gz)
    }
}

/// Read-pair cache. Let's us stream through the BAM and find nearby mates so we can write them out immediately
struct RpCache {
    cache: HashMap<Vec<u8>, Record>
}

impl RpCache {

    pub fn new() -> RpCache {
        RpCache { cache: HashMap::new() }
    }

    pub fn cache_rec(&mut self, rec: Record) -> Option<(Record, Record)> {
        // If cache already has entry, we have a pair! Return both
        match self.cache.remove(rec.qname()) {
            Some(old_rec) => {
                if rec.is_first_in_template() && old_rec.is_last_in_template() {
                    Some((rec, old_rec))
                } else if old_rec.is_first_in_template() && rec.is_last_in_template() {
                    Some((old_rec, rec))
                } else {
                    panic!("invalid pair")
                }
            },
            None => {
                self.cache.insert(Vec::from(rec.qname()), rec);
                None
            }
        }
    }

    pub fn clear_orphans(&mut self, current_tid: i32, current_pos: i32) -> Vec<Record> {
        let mut orphans = Vec::new();
        let mut new_cache = HashMap::new();

        for (key, rec) in self.cache.drain() {
            if rec.pos() - current_pos > 5000 || rec.tid() != current_tid {
                orphans.push(rec);
            } else {
                new_cache.insert(key, rec);
            }
        }

        self.cache = new_cache;
        orphans
    }

    pub fn len(&self) -> usize 
    {
        self.cache.len()
    }
}

fn main() {
    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.decode())
                         .unwrap_or_else(|e| e.exit());
    go(args, None);
}                


pub fn go(args: Args, cache_size: Option<usize>) {

    let cache_size = cache_size.unwrap_or(100000);
    let bam = bam::Reader::new(&args.arg_bam).ok().expect("Error opening BAM file");

    let formatter = {
        let header_fmt = FormatBamRecords::from_headers(&bam);
        match header_fmt {
            Some(f) => f,
            None => {
                if args.flag_gemcode {
                    FormatBamRecords::gemcode()
                } else if args.flag_lr20 {
                    FormatBamRecords::lr20()
                } else if args.flag_cr11 {
                    FormatBamRecords::cr11()        
                } else {
                    println!("Unrecognized 10x BAM file. For BAM files produced by older pipelines, use one of the following flags:");
                    println!("--gemcode   BAM files created with GemCode data using Longranger 1.0 - 1.3");
                    println!("--lr20      BAM files created with Longranger 2.0 using Chromium Genome data");
                    println!("--cr11      BAM files created with Cell Ranger 1.0-1.1 using Single Cell 3' v1 data");
                    return
                }
            }
        }
    };

    let out_path = Path::new(&args.arg_output_path);

    match create_dir(&out_path) {
        Err(msg) => { println!("Couldn't create output directory: {:?}.  Error: {}", out_path, msg); return},
        Ok(_) => (),
    }
    
    // prep output files
    let fq = FastqWriter::new(out_path, formatter.clone(), args.flag_sample, args.flag_reads_per_fastq);
 
    if formatter.is_double_ended() {
        proc_double_ended(bam, formatter, fq, cache_size);
    } else {
        proc_single_ended(bam, formatter, fq);
    }
}


fn proc_double_ended(bam: bam::Reader, formatter: FormatBamRecords, mut fq: FastqWriter, cache_size: usize) {
    // Temp file for hold unpaired reads. Will be cleaned up automatically.
    let tmp_file = NamedTempFile::new().unwrap();

    let total_read_pairs = {
        // Cache for efficiently finding local read pairs
        let mut rp_cache = RpCache::new();

        // For chimeric read piars that are showing up in different places, we will write these to disk for later use
        let w: ShardWriteManager<SerFq, SerFqImpl> = ShardWriteManager::new(tmp_file.path(), 256, 2, SerFqImpl{});
        let mut sender = w.get_sender();

        // Count total R1s observed, so we can make sure we've preserved all read pairs
        let mut total_read_pairs = 0;

        for _rec in bam.records() {
            let rec = _rec.unwrap();

            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            if rec.is_first_in_template() {
                total_read_pairs += 1;
            }

            // Save our current location
            let tid = rec.tid();
            let pos = rec.pos();

            match rp_cache.cache_rec(rec) {
                Some((r1,r2)) => {
                    let (fq1, fq2, fq_i1, fq_i2) = formatter.format_read_pair(&r1, &r2).unwrap();
                    fq.write(&fq1, &fq2, &fq_i1, &fq_i2);
                },
                None => ()
            }

            // If cache gets too big, clear out stragglers & serialize for later
            if rp_cache.len() > cache_size {
                for orphan in rp_cache.clear_orphans(tid, pos) {
                    let ser = formatter.bam_rec_to_ser(&orphan);
                    sender.send(ser);
                }
            }
        }

        for (_, orphan) in rp_cache.cache.drain() {
            let ser = formatter.bam_rec_to_ser(&orphan);
            sender.send(ser);
        }

        total_read_pairs
    };

    // Read back the shards, sort to find pairs, and write.
    let reader = ShardReader::open(tmp_file.path(), SerFqImpl{});

    let mut ncached = 0;
    for s in 0..reader.num_shards() {
        let mut data = reader.read_shard(s);
        data.sort();

        for (_, items) in &data.iter().group_by(|x| &x.rec.head) {
            // write out items
            let mut item_vec: Vec<_> = items.collect();
            if item_vec.len() != 2 {
                panic!("didn't get both reads!: {:?}", item_vec);
            }

            item_vec.sort_by_key(|x| x.read_num);
            let r1 = item_vec.swap_remove(0);
            let r2 = item_vec.swap_remove(0);
            fq.write(&r1.rec, &r2.rec, &r1.i1, &r1.i2);
            ncached += 1;
        }
    }

    // make sure we have the right number of output reads
    println!("Writing finished.  Observed {} read pairs. Wrote {} read pairs ({} cached)", total_read_pairs, fq.total_written, ncached);
}

fn proc_single_ended(bam: bam::Reader, formatter: FormatBamRecords, mut fq: FastqWriter) {

    let total_reads = {
        // Count total R1s observed, so we can make sure we've preserved all read pairs
        let mut total_reads = 0;

        for _rec in bam.records() {
            let rec = _rec.unwrap();

            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            total_reads += 1;

            let (fq1, fq2, fq_i1, fq_i2) = formatter.format_read(&rec).unwrap();
            fq.write(&fq1, &fq2, &fq_i1, &fq_i2);
        }

        total_reads
    };

    // make sure we have the right number of output reads
    println!("Writing finished.  Observed {} read pairs. Wrote {} read pairs", total_reads, fq.total_written);
}

#[cfg(test)]
mod tests {
    use tempdir;
    use std::path::Path;
    use super::*;
    use std::fs::{File};
    use std::io::{Read, BufRead, BufReader, Lines};
    use std::collections::HashMap;
    use std::boxed::Box;
    use flate2::read::GzDecoder;

    /// Head, Seq, Qual from single FASTQ
    type FqRec = (Vec<u8>, Vec<u8>, Vec<u8>);

    /// R1, R2, optional SI
    type RawReadSet = (FqRec, FqRec, Option<FqRec>);

    /// Read a single FqRec from a line 
    pub fn get_fastq_item<R: BufRead>(lines: &mut Lines<R>) -> Option<FqRec>{
        match lines.next() {
            Some(head) => {
                let r1 = lines.next().unwrap().unwrap().into_bytes();
                let _  = lines.next().unwrap().unwrap().into_bytes();
                let q1 = lines.next().unwrap().unwrap().into_bytes();

                // trim after first space of head
                let head_str = head.unwrap();
                let mut split = head_str.split_whitespace();
                let trim_head = split.next().unwrap().to_string().into_bytes();

                Some((trim_head, r1, q1))
            },
            None => None
        }
    }

    pub struct FastqPairIter {
        r1: Lines<BufReader<Box<Read>>>,
        r2: Lines<BufReader<Box<Read>>>,
        si: Option<Lines<BufReader<Box<Read>>>>,
    }

    pub fn open_w_gz<P: AsRef<Path>>(p: P) -> Box<Read> {
        let r = File::open(p.as_ref()).unwrap();

        if p.as_ref().extension().unwrap() == "gz" {
            Box::new(GzDecoder::new(r).unwrap())
        } else {
            Box::new(r)
        }
    }


    pub fn open_fastq_pair_iter<P: AsRef<Path>>(r1: P, r2: P, si: Option<P>) -> Box<Iterator<Item=RawReadSet>> {
        Box::new(
        FastqPairIter::init(
            open_w_gz(r1),
            open_w_gz(r2),
            si.map(|si| open_w_gz(si)),
        ))
    }

    impl FastqPairIter {
        pub fn init(r1: Box<Read>, r2: Box<Read>, si: Option<Box<Read>>) -> FastqPairIter {

            FastqPairIter {
                r1: BufReader::new(r1).lines(),
                r2: BufReader::new(r2).lines(),
                si: si.map(|x| BufReader::new(x).lines()),
            }
        }
    }

    impl Iterator for FastqPairIter {
        type Item = RawReadSet;

        fn next(&mut self) -> Option<RawReadSet> {

            match get_fastq_item(&mut self.r1) {
                Some(f_r1) => {
                    let f_r2 = get_fastq_item(&mut self.r2).unwrap();
                    
                    let f_si = self.si.as_mut().map(|_si| get_fastq_item(_si).unwrap());
                    Some((f_r1, f_r2, f_si))
                }
                None => None,
            }
        }
    }

    pub struct InterleavedFastqPairIter {
        ra: Lines<BufReader<Box<Read>>>,
        si: Option<Lines<BufReader<Box<Read>>>>,
    }

    pub fn open_interleaved_fastq_pair_iter<P: AsRef<Path>>(ra: P, si: Option<P>) -> Box<Iterator<Item=RawReadSet>> {
        Box::new(
        InterleavedFastqPairIter::init(
            open_w_gz(ra),
            si.map(|si| open_w_gz(si)),
        ))
    }

    impl InterleavedFastqPairIter {

        pub fn init(ra: Box<Read>, si: Option<Box<Read>>) -> InterleavedFastqPairIter {
            InterleavedFastqPairIter {
                ra: BufReader::new(ra).lines(),
                si: si.map(|x| BufReader::new(x).lines()),
            }
        }
    }

    impl Iterator for InterleavedFastqPairIter {
        type Item = RawReadSet;

        fn next(&mut self) -> Option<RawReadSet> {

            match get_fastq_item(&mut self.ra) {
                Some(f_r1) => {
                    let f_r2 = get_fastq_item(&mut self.ra).unwrap();
                    let f_si = self.si.as_mut().map(|_si| get_fastq_item(_si).unwrap());
                    Some((f_r1, f_r2, f_si))
                }
                None => None,
            }
        }
    }

    type ReadSet = HashMap<Vec<u8>, RawReadSet>;

    pub fn load_fastq_set<I: Iterator<Item=RawReadSet>>(iter: I) -> ReadSet {
        let mut reads = ReadSet::new();

        for r in iter {
            reads.insert((r.0).0.clone(), r);
        }

        reads
    }

    pub fn compare_read_sets(orig_set: ReadSet, new_set: ReadSet) {
        assert_eq!(orig_set.len(), new_set.len());

        let mut keys1: Vec<Vec<u8>> = orig_set.keys().cloned().collect();
        keys1.sort();

        let mut keys2: Vec<Vec<u8>> = new_set.keys().cloned().collect();
        keys2.sort();

        for (k1, k2) in keys1.iter().zip(keys2.iter()) {
            assert_eq!(k1, k2);
            assert_eq!(orig_set.get(k1), new_set.get(k2));
        }

        assert_eq!(orig_set, new_set);
    }

    #[test]
    fn test_lr21() {
        let tempdir = tempdir::TempDir::new("bam_to_fq_test").expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");

        let args = Args {
            arg_bam: "test/lr21.bam".to_string(),
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            flag_gemcode: false,
            flag_lr20: false,
            flag_cr11: false,
            flag_sample: "sample".to_string(),
            flag_reads_per_fastq: 100000,
        };

        super::go(args, Some(2));

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/crg-tiny-fastq-2.0.0/read-RA_si-GTTGCAGC_lane-001-chunk-001.fastq.gz", 
            Some("test/crg-tiny-fastq-2.0.0/read-I1_si-GTTGCAGC_lane-001-chunk-001.fastq.gz"));

        let orig_set = load_fastq_set(true_fastq_read);

        let new_set = load_fastq_set(
            open_fastq_pair_iter(
                tmp_path.join("sample_S1_L001_R1_001.fastq.gz"), 
                tmp_path.join("sample_S1_L001_R2_001.fastq.gz"), 
                Some(tmp_path.join("sample_S1_L001_I1_001.fastq.gz"))));
        
        compare_read_sets(orig_set, new_set);
    }


    #[test]
    fn test_cr12() {
        let tempdir = tempdir::TempDir::new("bam_to_fq_test").expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");

        let args = Args {
            arg_bam: "test/cr12.bam".to_string(),
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            flag_gemcode: false,
            flag_lr20: false,
            flag_cr11: false,
            flag_sample: "sample".to_string(),
            flag_reads_per_fastq: 100000,
        };

        super::go(args, Some(2));

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/cellranger-tiny-fastq-1.2.0/read-RA_si-TTTCATGA_lane-008-chunk-001.fastq.gz",
            Some("test/cellranger-tiny-fastq-1.2.0/read-I1_si-TTTCATGA_lane-008-chunk-001.fastq.gz"));

        let orig_set = load_fastq_set(true_fastq_read);

        let new_set = load_fastq_set(
            open_fastq_pair_iter(
                tmp_path.join("sample_S1_L001_R1_001.fastq.gz"), 
                tmp_path.join("sample_S1_L001_R2_001.fastq.gz"), 
                Some(tmp_path.join("sample_S1_L001_I1_001.fastq.gz"))));
        
        compare_read_sets(orig_set, new_set);
    }
}
