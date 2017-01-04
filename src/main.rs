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

#[cfg(test)]
extern crate fastq_10x;

use std::io::{Write, BufWriter, Result};
use std::fs::File;
use std::fs::create_dir;
use std::path::{PathBuf, Path};
use std::hash::{Hash, SipHasher, Hasher};
use std::str::FromStr;

use std::collections::HashMap;
use itertools::Itertools;

use tempfile::NamedTempFile;

use flate2::write::GzEncoder;
use flate2::Compression;

use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::{Aux, Record};

use bincode::rustc_serialize::{encode_into, decode};
use shardio::shard::{Serializer, Shardable, ShardWriteManager, ShardReader};

use regex::Regex;
use docopt::Docopt;

mod locus;

const USAGE: &'static str = "
10x Genomics BAM to FASTQ converter.

    Tool for converting 10x BAMs produced by Cell Ranger or Long Ranger back to
    FASTQ files that can be used as inputs to re-run analysis. The FASTQ files
    emitted by the tool should contain the same set of sequences that were
    input to the original pipeline run, although the order will not be 
    preserved.  The FASTQs will be emitted into a directory structure that is
    compatible with the directories created by the 'mkfastq' tool.

    10x BAMs produced by Long Ranger v2.1+ and Cell Ranger v1.2+ contain header
    fields that permit automatic conversion to the correct FASTQ sequences.

    Older 10x pipelines require one of the arguments listed below to indicate 
    which pipeline created the BAM.


Usage:
  bamtofastq [options] <bam> <output-path>
  bamtofastq (-h | --help)

Options:
  --locus=<locus>      Optional. Only include read pairs mapping to locus. Use chrom:start-end format.
  --reads-per-fastq=N  Number of reads per FASTQ chunk [default: 200000000]
  --gemcode            Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)
  --lr20               Convert a BAM produced by Longranger 2.0
  --cr11               Convert a BAM produced by Cell Ranger 1.0-1.1
  -h --help            Show this screen.
";

#[derive(Debug, RustcDecodable, Clone)]
pub struct Args {
    arg_bam: String,
    arg_output_path: String,
    flag_locus: Option<String>,
    flag_reads_per_fastq: usize,
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
    read_group: Option<Rg>,
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

type Rg = (String, u32);

/// Spec for converting from a BAM record back to reads. Empty vector indicates that this read doesn't exist
/// in the output. The i1 and i2 reads should be buildable from tags in the R1 read.
#[derive(Clone)]
struct FormatBamRecords {
    rg_spec: Option<HashMap<String, Rg>>,
    r1_spec: Vec<SpecEntry>,
    r2_spec: Vec<SpecEntry>,
    i1_spec: Vec<SpecEntry>,
    i2_spec: Vec<SpecEntry>,
    rename: Option<Vec<String>>,
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
    pub fn from_headers<R: bam::Read>(reader: &R) -> Option<FormatBamRecords> {

        let mut spec = Self::parse_spec(reader);

        if spec.len() == 0 {
            None
        } else {
            Some(
                FormatBamRecords {
                    rg_spec: Self::parse_rgs(reader),
                    r1_spec: spec.remove("R1").unwrap(),
                    r2_spec: spec.remove("R2").unwrap(),
                    i1_spec: spec.remove("I1").unwrap_or_else(|| Vec::new()),
                    i2_spec: spec.remove("I2").unwrap_or_else(|| Vec::new()),
                    rename: None,
            })
        }
    }

    /// hard-coded spec for gemcode BAM files
    pub fn gemcode<R: bam::Read>(reader: &R) -> FormatBamRecords {

        FormatBamRecords {
            rg_spec: Self::parse_rgs(reader),
            r1_spec: vec![SpecEntry::Read],
            r2_spec: vec![SpecEntry::Read],
            i1_spec: vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())],
            i2_spec: vec![SpecEntry::Tags("RX".to_string(), "QX".to_string())],
            rename: None,
        }
    }

    // hard-coded specs for longranger 2.0 BAM files
    pub fn lr20<R: bam::Read>(reader: &R) -> FormatBamRecords {

        FormatBamRecords {
            rg_spec: Self::parse_rgs(reader),
            r1_spec: vec![SpecEntry::Tags("RX".to_string(), "QX".to_string()), SpecEntry::Ns(7), SpecEntry::Read],
            r2_spec: vec![SpecEntry::Read],
            i1_spec: vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())],
            i2_spec: vec![],
            rename: None,
        }
    } 

    // hard-coded specs for cellranger 1.0-1.1 BAM files
    pub fn cr11<R: bam::Read>(reader: &R) -> FormatBamRecords {

        FormatBamRecords {
            rg_spec: Self::parse_rgs(reader),
            r1_spec: vec![SpecEntry::Read],
            r2_spec: vec![SpecEntry::Tags("UR".to_string(), "UQ".to_string())],
            i1_spec: vec![SpecEntry::Tags("CR".to_string(), "CQ".to_string())],
            i2_spec: vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())],
            rename: Some(vec!["R1".to_string(), "R3".to_string(), "R2".to_string(), "I1".to_string()]),
        }
    }

    fn parse_rgs<R: bam::Read>(reader: &R) -> Option<HashMap<String, Rg>> {
        let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();
        let mut rg_items = HashMap::new();

        for l in text.lines() {
            if l.starts_with("@RG") {
                let r = Self::parse_rg_line(l);
                match r {
                    Some((id, rg, lane)) => {rg_items.insert(id, (rg, lane));},
                    None => ()
                };
            }
        }

        if rg_items.len() > 0 { Some(rg_items) } else { None }
    }

    fn parse_rg_line(line: &str) -> Option<(String, String, u32)> {
        
        let mut entries = line.split('\t');
        let _ = entries.next(); // consume @RG entry

        let mut tags = HashMap::new();
        for entry in entries {
            let mut parts = entry.splitn(2, ':');
            let tag = parts.next().unwrap();
            let val = parts.next().unwrap();
            tags.insert(tag.to_string(), val.to_string());
        }

        if tags.contains_key("ID") {
            let v = tags.remove("ID").unwrap();

            let vv = &v;
            let mut parts = vv.rsplitn(2, ':');
            let lane = parts.next().unwrap();
            let rg = parts.next().unwrap().to_string();

            Some((v.clone(), rg, u32::from_str(lane).unwrap()))
        } else {
            None
        }
    }


    /// Parse the specs from BAM headers if available
    fn parse_spec<R: bam::Read>(reader: &R) -> HashMap<String, Vec<SpecEntry>> {

        // Example header line:
        // @CO	10x_bam_to_fastq:R1(RX:QX,TR:TQ,SEQ:QUAL)
        let re = Regex::new(r"@CO\t10x_bam_to_fastq:(\S+)\((\S+)\)").unwrap();
        let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();
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

    pub fn find_rg(&self, rec: &Record) -> Option<Rg> {
        match self.rg_spec {
            Some(ref spec) => {
                let rg = rec.aux(b"RG");
                match rg {
                    Some(Aux::String(s)) => { 
                        let key = String::from_utf8(Vec::from(s)).unwrap();
                        spec.get(&key).map(|x| x.clone())
                    },
                    _ => None,
                }
            },
            None => None,
        }
    }

    /// Convert a BAM record to a Fq record, for internal caching
    pub fn bam_rec_to_ser(&self, rec: &Record) -> SerFq {
        match (rec.is_first_in_template(), rec.is_last_in_template()) {
            (true, false) => {
                SerFq {
                    read_group: self.find_rg(rec),
                    read_num: ReadNum::R1,
                    rec: self.bam_rec_to_fq(rec, &self.r1_spec).unwrap(),
                    i1: if self.i1_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i1_spec).unwrap()) }  else { None },
                    i2: if self.i2_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i2_spec).unwrap()) }  else { None },
                    
                }
            },
            (false, true) => {
                SerFq {
                    read_group: self.find_rg(rec),
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

    pub fn format_read_pair(&self, r1_rec: &Record, r2_rec: &Record) -> Result<(Option<Rg>, FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)> {
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

        let rg = self.find_rg(r1_rec);
        Ok((rg, r1, r2, i1, i2))
    }


    pub fn format_read(&self, rec: &Record) -> Result<(Option<Rg>, FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)> {
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

        let rg = self.find_rg(rec);
        Ok((rg, r1, r2, i1, i2))
    }

    /// A spec implies double-ended reads if both the R1 and R2 reads generate different BAM records.
    /// If not the R1 and R2 sequences can be derived from a single BAM entry.
    pub fn is_double_ended(&self) -> bool {
        self.r1_spec.contains(&SpecEntry::Read) && self.r2_spec.contains(&SpecEntry::Read)
    }
}

type BGW = BufWriter<GzEncoder<File>>;

struct FastqManager {
    writers: HashMap<Rg, FastqWriter>,
}

impl FastqManager {
    pub fn new(out_path: &Path, formatter: FormatBamRecords, sample_name: String, reads_per_fastq: usize) -> FastqManager {

        match formatter.rg_spec {
            Some(ref spec) => {
                // Take the read groups and generate read group paths
                let mut sample_def_paths = HashMap::new();
                let mut writers = HashMap::new();

                for (_, &(ref _samp, lane)) in spec.iter() {
                    let samp = _samp.clone();
                    let path = sample_def_paths.entry(samp).or_insert_with(|| {
                        let suffix = _samp.replace(":", "_");
                        let samp_path = out_path.join(suffix);
                        create_dir(&samp_path);
                        samp_path
                    });
                    
                    let writer = FastqWriter::new(path, formatter.clone(), "bamtofastq".to_string(), lane, reads_per_fastq);
                    writers.insert((_samp.clone(), lane), writer);
                }

                FastqManager { writers: writers }
            },
            None => {
                let mut w = HashMap::new();
                w.insert(("def".to_string(), 0), FastqWriter::new(out_path, formatter.clone(), "bamtofastq".to_string(), 1, reads_per_fastq));
                FastqManager { writers: w }
            }
        }
    }

    pub fn write(&mut self, rg: &Option<Rg>, r1: &FqRecord, r2: &FqRecord, i1: &Option<FqRecord>, i2: &Option<FqRecord>) {
        let w = self.writers.get_mut(rg.as_ref().unwrap_or(&("def".to_string(), 0))).unwrap();
        w.write(r1, r2, i1, i2);
    }

    pub fn total_written(&self) -> usize {
        self.writers.iter().map(|(_, ref w)| w.total_written).sum()
    }

    pub fn paths(&self) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {
        let mut r = Vec::new();

        for (_, w) in self.writers.iter() {
            r.extend(w.path_sets.clone());
        }

        r
    }
}


/// Open Fastq files being written to
struct FastqWriter {
    formatter: FormatBamRecords,
    out_path: PathBuf,
    sample_name: String,
    lane: u32,

    r1: BGW,
    r2: BGW,
    i1: Option<BGW>,
    i2: Option<BGW>,

    chunk_written: usize,
    total_written: usize,
    n_chunks: usize,
    reads_per_fastq: usize,
    path_sets: Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>,
}

/// Write sets of Fastq records to the open Fastq files
impl FastqWriter {

    pub fn new(out_path: &Path, formatter: FormatBamRecords, sample_name: String, lane: u32, reads_per_fastq: usize) -> FastqWriter {

        // open output files
        let paths = Self::get_paths(&out_path, &sample_name, lane, 0, &formatter);
        //let r1 = 
        let r2 = Self::open_gzip_writer(&paths.1);
        let i1 = paths.2.as_ref().map(|p| Self::open_gzip_writer(p));
        let i2 = paths.3.as_ref().map(|p| Self::open_gzip_writer(p));

        FastqWriter {
            formatter: formatter,
            out_path: out_path.to_path_buf(),
            sample_name: sample_name,
            lane: lane,
            r1: Self::open_gzip_writer(&paths.0),
            r2: r2,
            i1: i1,
            i2: i2,
            n_chunks: 1,
            total_written: 0,
            chunk_written: 0,
            reads_per_fastq: reads_per_fastq,
            path_sets: vec![paths],
        }
    }

    fn get_paths(out_path: &Path, sample_name: &str, lane: u32, n_files: usize, formatter: &FormatBamRecords) -> (PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>) {

        if formatter.rename.is_none() {
            let r1 = out_path.join(format!("{}_S1_L{:03}_R1_{:03}.fastq.gz", sample_name, lane, n_files+1));
            let r2 = out_path.join(format!("{}_S1_L{:03}_R2_{:03}.fastq.gz", sample_name, lane, n_files+1));
            let i1 = out_path.join(format!("{}_S1_L{:03}_I1_{:03}.fastq.gz", sample_name, lane, n_files+1));
            let i2 = out_path.join(format!("{}_S1_L{:03}_I2_{:03}.fastq.gz", sample_name, lane, n_files+1));

            (r1, r2, 
            if formatter.i1_spec.len() > 0 { Some(i1) } else { None }, 
            if formatter.i2_spec.len() > 0 { Some(i2) } else { None })
        } else {
            let new_read_names = formatter.rename.as_ref().unwrap();

            let r1 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, n_files+1, new_read_names[0]));
            let r2 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, n_files+1, new_read_names[1]));
            let i1 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, n_files+1, new_read_names[2]));
            let i2 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, n_files+1, new_read_names[3]));

            (r1, r2, 
            if formatter.i1_spec.len() > 0 { Some(i1) } else { None }, 
            if formatter.i2_spec.len() > 0 { Some(i2) } else { None })

        }
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
        let paths = Self::get_paths(&self.out_path, &self.sample_name, self.lane, self.n_chunks, &self.formatter);
        self.r1 = Self::open_gzip_writer(&paths.0);
        self.r2 = Self::open_gzip_writer(&paths.1);
        self.i1 = paths.2.as_ref().map(|p| Self::open_gzip_writer(p));
        self.i2 = paths.3.as_ref().map(|p| Self::open_gzip_writer(p));
        self.n_chunks += 1;
        self.chunk_written = 0;
        self.path_sets.push(paths);
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


pub fn go(args: Args, cache_size: Option<usize>) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {

    let cache_size = cache_size.unwrap_or(100000);
 
    match args.flag_locus {
        Some(ref locus) => {
            let loc = locus::Locus::from_str(locus).unwrap();
            let mut bam = bam::IndexedReader::from_path(&args.arg_bam).ok().expect("Error opening BAM file");
            let tid = bam.header().tid(loc.chrom.as_bytes()).expect("Unrecognized chromosome in locus");
            bam.seek(tid, loc.start, loc.end);
            inner(args.clone(), cache_size, bam)
        },
        None => {
            let bam = bam::Reader::from_path(&args.arg_bam).ok().expect("Error opening BAM file");
            inner(args, cache_size, bam)
        }
    }
}

pub fn inner<R: bam::Read>(args: Args, cache_size: usize, bam: R) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {

    let formatter = {
        let header_fmt = FormatBamRecords::from_headers(&bam);
        match header_fmt {
            Some(f) => f,
            None => {
                if args.flag_gemcode {
                    FormatBamRecords::gemcode(&bam)
                } else if args.flag_lr20 {
                    FormatBamRecords::lr20(&bam)
                } else if args.flag_cr11 {
                    FormatBamRecords::cr11(&bam)     
                } else {
                    println!("Unrecognized 10x BAM file. For BAM files produced by older pipelines, use one of the following flags:");
                    println!("--gemcode   BAM files created with GemCode data using Longranger 1.0 - 1.3");
                    println!("--lr20      BAM files created with Longranger 2.0 using Chromium Genome data");
                    println!("--cr11      BAM files created with Cell Ranger 1.0-1.1 using Single Cell 3' v1 data");
                    return vec![];
                }
            }
        }
    };

    let out_path = Path::new(&args.arg_output_path);

    match create_dir(&out_path) {
        Err(msg) => { println!("Couldn't create output directory: {:?}.  Error: {}", out_path, msg); return vec![]},
        Ok(_) => (),
    }
    
    // prep output files
    println!("{:?}", args);
    let fq = FastqManager::new(out_path, formatter.clone(), "bamtofastq".to_string(), args.flag_reads_per_fastq);
 
    if formatter.is_double_ended() {
        proc_double_ended(bam, formatter, fq, cache_size)
    } else {
        proc_single_ended(bam, formatter, fq)
    }
}


fn proc_double_ended<R: bam::Read>(bam: R, formatter: FormatBamRecords, mut fq: FastqManager, cache_size: usize) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {
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
                    let (rg, fq1, fq2, fq_i1, fq_i2) = formatter.format_read_pair(&r1, &r2).unwrap();
                    fq.write(&rg, &fq1, &fq2, &fq_i1, &fq_i2);
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
            fq.write(&r1.read_group, &r1.rec, &r2.rec, &r1.i1, &r1.i2);
            ncached += 1;
        }
    }

    // make sure we have the right number of output reads
    println!("Writing finished.  Observed {} read pairs. Wrote {} read pairs ({} cached)", total_read_pairs, fq.total_written(), ncached);
    fq.paths()
}

fn proc_single_ended<R: bam::Read>(bam: R, formatter: FormatBamRecords, mut fq: FastqManager) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {

    let total_reads = {
        // Count total R1s observed, so we can make sure we've preserved all read pairs
        let mut total_reads = 0;

        for _rec in bam.records() {
            let rec = _rec.unwrap();

            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            total_reads += 1;

            let (rg, fq1, fq2, fq_i1, fq_i2) = formatter.format_read(&rec).unwrap();
            fq.write(&rg, &fq1, &fq2, &fq_i1, &fq_i2);
        }

        total_reads
    };

    // make sure we have the right number of output reads
    println!("Writing finished.  Observed {} read pairs. Wrote {} read pairs", total_reads, fq.total_written());
    fq.paths()
}

#[cfg(test)]
mod tests {
    use tempdir;
    use super::*;
    use fastq_10x::*;
    use std::collections::HashMap;

    type ReadSet = HashMap<Vec<u8>, RawReadSet>;

    pub fn load_fastq_set<I: Iterator<Item=RawReadSet>>(reads: &mut ReadSet, iter: I) {
        for r in iter {
            reads.insert((r.0).0.clone(), r);
        }
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

    pub fn compare_read_sets_ignore_n(orig_set: ReadSet, new_set: ReadSet) {
        assert_eq!(orig_set.len(), new_set.len());

        let mut keys1: Vec<Vec<u8>> = orig_set.keys().cloned().collect();
        keys1.sort();

        let mut keys2: Vec<Vec<u8>> = new_set.keys().cloned().collect();
        keys2.sort();

        assert_eq!(keys1, keys2);


        for (k1, k2) in keys1.iter().zip(keys2.iter()) {
            assert_eq!(k1, k2);
            compare_raw_read_sets_ignore_n(orig_set.get(k1).unwrap(), new_set.get(k2).unwrap());
        }
    }

    // Relax the comparison for R1 -- if the v2 R1 read has 'N' or the v2 R1 qual has 'J', allow it through
    // this handles the case where the 7 trimmed bases after the BC are were not retained in Long Ranger 2.0
    // also ignore mismatches in the first 16bp, which are caused by the bug in LR 2.0 that caused the RX
    // tag to have the corrected sequence rather than the raw sequence
    pub fn compare_raw_read_sets_ignore_n(v1: &RawReadSet, v2: &RawReadSet) {
        assert_eq!(&(v1.0).0, &(v2.0).0);
        compare_bytes_ignore_n(&(v1.0).1, &(v2.0).1);
        compare_bytes_ignore_n(&(v1.0).2, &(v2.0).2);

        assert_eq!(v1.1, v2.1);
        assert_eq!(v1.2, v2.2)
    }


    pub fn compare_bytes_ignore_n(v1: &Vec<u8>, v2: &Vec<u8>) {
        assert_eq!(v1.len(), v2.len());
        for (idx, (b1, b2)) in v1.iter().zip(v2).enumerate() {
            if idx >= 16 && b1 != b2 && *b2 != b'N' && *b2 != b'J' {
                println!("got mismatch at pos: {}", idx);
                assert_eq!(v1, v2)
            }
        }
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
            flag_reads_per_fastq: 100000,
            flag_locus: None,
        };

        let out_path_sets = super::go(args, Some(2));

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/crg-tiny-fastq-2.0.0/read-RA_si-GTTGCAGC_lane-001-chunk-001.fastq.gz", 
            Some("test/crg-tiny-fastq-2.0.0/read-I1_si-GTTGCAGC_lane-001-chunk-001.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, true_fastq_read);

        let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1, r2, i1));
        }
        
        compare_read_sets(orig_reads, output_reads);
    }

    #[test]
    fn test_lr20() {
        let tempdir = tempdir::TempDir::new("bam_to_fq_test").expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");

        let args = Args {
            arg_bam: "test/lr20.bam".to_string(),
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            flag_gemcode: false,
            flag_lr20: true,
            flag_cr11: false,
            flag_reads_per_fastq: 100000,
        };

        let out_path_sets = super::go(args, Some(2));

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/crg-tiny-fastq-2.0.0/read-RA_si-GTTGCAGC_lane-001-chunk-001.fastq.gz", 
            Some("test/crg-tiny-fastq-2.0.0/read-I1_si-GTTGCAGC_lane-001-chunk-001.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, true_fastq_read);

        let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1, r2, i1));
        }
        
        // use special comparison method that ignores N's in R1
        // accounts for missing trimmed bases
        compare_read_sets_ignore_n(orig_reads, output_reads);
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
            flag_reads_per_fastq: 100000,
            flag_locus: None,
        };

        let out_path_sets = super::go(args, Some(2));

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/cellranger-tiny-fastq-1.2.0/read-RA_si-TTTCATGA_lane-008-chunk-001.fastq.gz",
            Some("test/cellranger-tiny-fastq-1.2.0/read-I1_si-TTTCATGA_lane-008-chunk-001.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, true_fastq_read);

       let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1, r2, i1));
        }
        
        compare_read_sets(orig_reads, output_reads);
    }

    #[test]
    fn test_cr12_v1() {
        let tempdir = tempdir::TempDir::new("bam_to_fq_test").expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");

        let args = Args {
            arg_bam: "test/cr12-v1.bam".to_string(),
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            flag_gemcode: false,
            flag_lr20: false,
            flag_cr11: false,
            flag_reads_per_fastq: 100000,
        };

        let out_path_sets = super::go(args, Some(2));

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/cellranger-3p-v1/read-RA_si-ACCAGTCC_lane-001-chunk-000.fastq.gz",
            Some("test/cellranger-3p-v1/read-I1_si-ACCAGTCC_lane-001-chunk-000.fastq.gz"),
        );

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, true_fastq_read);

       let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets.clone() {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1, r2, i1));
        }
        
        compare_read_sets(orig_reads, output_reads);

        // Separately test I1 & I2 as if they were the main reads.
        let true_index_reads = open_fastq_pair_iter(
            "test/cellranger-3p-v1/read-I1_si-ACCAGTCC_lane-001-chunk-000.fastq.gz",
            "test/cellranger-3p-v1/read-I2_si-ACCAGTCC_lane-001-chunk-000.fastq.gz",
            None);
        let mut orig_index_reads = ReadSet::new();
        load_fastq_set(&mut orig_index_reads, true_index_reads);


        let mut output_index_reads = ReadSet::new();
        for (_, _, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_index_reads, open_fastq_pair_iter(i1.unwrap(), i2.unwrap(), None));
        }

        compare_read_sets(orig_index_reads, output_index_reads);
    }
}
