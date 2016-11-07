extern crate docopt;
extern crate rust_htslib;
extern crate flate2;
extern crate shardio;
extern crate bincode;
extern crate itertools;
extern crate rustc_serialize;
extern crate regex;
extern crate tempfile;

use std::io::{Write, BufWriter, Result};
use std::fs::File;
use std::fs::create_dir;
use std::path::Path;
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
  --gemcode            Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)
  --lr20               Convert a BAM produced by Longranger 2.0
  --cr11               Convert a BAM produced by Cell Ranger 1.0-1.1
  -h --help            Show this screen.
  --version            Show version.
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_bam: String,
    arg_output_path: String,
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
#[derive(Debug)]
enum SpecEntry {
    Tags(String, String),
    Ns(usize),
    Read,
}

/// Spec for converting from a BAM record back to reads. Empty vector indicates that this read doesn't exist
/// in the output. The i1 and i2 reads should be buildable from tags in the R1 read.
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
        let text = reader.header.text();
        let text = text.unwrap();
        let mut spec = HashMap::new();

        for l in text.lines() {
            match re.captures(l) {
                Some(c) => {
                    let mut read_spec = Vec::new();

                    println!("got de-bam header: {:?}", c);
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

        for item in spec {
            match item {
                // Data from a tag
                &SpecEntry::Tags(ref read_tag, ref qv_tag) => {                
                    let rx = rec.aux(read_tag.as_bytes()).unwrap().string();
                    r.extend_from_slice(rx);

                    let qx = rec.aux(qv_tag.as_bytes()).unwrap().string();
                    q.extend_from_slice(qx);   
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
}

/// Open Fastq files being written to
struct FastqWriter<W: Write> {
    r1: W,
    r2: W,
    i1: Option<W>,
    i2: Option<W>,
    total_written: usize,
}

/// Write sets of Fastq records to the open Fastq files
impl<W: Write> FastqWriter<W> {

    pub fn write_rec(w: &mut W, rec: &FqRecord)  {
        w.write(b"@").unwrap();
        w.write(&rec.head).unwrap();
        w.write(b"\n").unwrap();

        w.write(&rec.seq).unwrap();
        w.write(b"\n+\n").unwrap();
        w.write(&rec.qual).unwrap();
        w.write(b"\n").unwrap();
    }

    pub fn try_write_rec(w: &mut Option<W>, rec: &Option<FqRecord>) {
        match w {
            &mut Some(ref mut w) => match rec { &Some(ref r) => FastqWriter::write_rec(w, r), &None => panic!("setup error") },
            &mut None => ()
        }
    }

    pub fn write(&mut self, r1: &FqRecord, r2: &FqRecord, i1: &Option<FqRecord>, i2: &Option<FqRecord>) {
        self.total_written += 1;

        FastqWriter::write_rec(&mut self.r1, r1);
        FastqWriter::write_rec(&mut self.r2, r2);

        FastqWriter::try_write_rec(&mut self.i1, i1);
        FastqWriter::try_write_rec(&mut self.i2, i2);
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

fn make_writer<P: AsRef<Path>>(path: P) -> BufWriter<GzEncoder<File>> {
    let f = File::create(path).unwrap();
    let gz = GzEncoder::new(f, Compression::Fast);
    BufWriter::new(gz)
}


fn main() {
    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.decode())
                         .unwrap_or_else(|e| e.exit());


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
        Err(msg) => println!("Couldn't create output directory: {:?}.  Error: {}", out_path, msg),
        Ok(_) => (),
    }
    
    // open output files
    let r1_path = out_path.join(Path::new("r1.fastq.gz"));
    let r2_path = out_path.join(Path::new("r2.fastq.gz"));
    let i1_path = out_path.join(Path::new("i1.fastq.gz"));
    let i2_path = out_path.join(Path::new("i2.fastq.gz"));

    let mut fq = FastqWriter {
        r1: make_writer(r1_path),
        r2: make_writer(r2_path),
        i1: if formatter.i1_spec.len() > 0 { Some(make_writer(i1_path)) } else { None },
        i2: if formatter.i2_spec.len() > 0 { Some(make_writer(i2_path)) } else { None },
        total_written: 0,
    };


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
            if rp_cache.len() > 1000000 {
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
        }
    }

    // make sure we have the right number of output reads
    println!("Writing finished.  Observed {} read pairs. Wrote {} read pairs", total_read_pairs, fq.total_written);
}
