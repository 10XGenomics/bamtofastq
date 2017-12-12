// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.

// `error_chain!` can recurse deeply
#![recursion_limit = "1024"]

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate serde_derive;

extern crate docopt;
extern crate rust_htslib;
extern crate flate2;
extern crate shardio;
extern crate bincode;
extern crate itertools;
extern crate regex;
extern crate tempfile;
extern crate tempdir;
extern crate serde;

use std::io::{Write, BufWriter};
use std::fs::File;
use std::fs::create_dir;
use std::path::{PathBuf, Path};
use std::hash::{Hash, SipHasher, Hasher};
use std::str::FromStr;
use std::result;

use std::collections::HashMap;

use tempfile::NamedTempFile;
use itertools::Itertools;

use flate2::write::GzEncoder;
use flate2::Compression;

use rust_htslib::bam::{self, Read, ReadError};
use rust_htslib::bam::record::{Aux, Record};

use shardio::{ShardDef, ShardWriteManager, ShardReader};
use shardio::helper::ThreadProxyWriter;

use regex::Regex;
use docopt::Docopt;

mod locus;
mod rpcache;
mod bx_index;
mod fastq_reader;
use rpcache::RpCache;
use bx_index::BxListIter;


// We'll put our errors in an `errors` module, and other modules in
// this crate will `use errors::*;` to get access to everything
// `error_chain!` creates.
mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! { 
        foreign_links {
            Io(::std::io::Error) #[doc = "Link to a `std::error::Error` type."];
            Bam(::rust_htslib::bam::FetchError);
        }
    }
}

use errors::*;



const VERSION: &'static str = env!("CARGO_PKG_VERSION");

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

    NOTE: BAMs created by non-10x pipelines are unlikely to work correctly,
    unless all the relevant tags have been recreated.

    NOTE: BAM produced by the BASIC and ALIGNER pipeline from Long Ranger 2.1.2 and earlier
    are not compatible with bamtofastq

    NOTE: BAM files created by CR < 1.3 do not have @RG headers, so bamtofastq will use the GEM group
    annotations attached to the CB (cell barcode) tag to split data from multiple input libraries.
    Reads without a valid barcode do not carry the CB tag and will be dropped. These reads would 
    not be included in any valid cell.

Usage:
  bamtofastq [options] <bam> <output-path>
  bamtofastq (-h | --help)

Options:
  --locus=<locus>      Optional. Only include read pairs mapping to locus. Use chrom:start-end format.
  --reads-per-fastq=N  Number of reads per FASTQ chunk [default: 50000000]
  --gemcode            Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)
  --lr20               Convert a BAM produced by Longranger 2.0
  --cr11               Convert a BAM produced by Cell Ranger 1.0-1.1
  --bx-list=L          Only include BX values listed in text file L. Requires BX-sorted and index BAM file (see Long Ranger support for details).
  -h --help            Show this screen.
";


/*
== Dev Notes ==
This code has bunch of special cases to workaround deficiences of the BAM files produced by older pipelines,
before we started enforcing the notion that BAM files should be easily convertible back to the original 
FASTQ data that was input.  Once these older pipelines & chemistries are out of service, the code could
be made considerably simpler.

1) Workaround for CR < 1.3:  there are no RG headers or RG tags, so we can't accurately get back to 
per-Gem Group FASTQs, which is important because multi-gem-group experiments are common.  If we don't
have RG headers, we will set up files for 20 gem groups, and use the gem-group suffix on the CB tag to
determine the Gem group.  Reads without a CB tag will get dropped.


*/




#[derive(Debug, Deserialize, Clone)]
pub struct Args {
    arg_bam: String,
    arg_output_path: String,
    flag_locus: Option<String>,
    flag_bx_list: Option<String>,
    flag_reads_per_fastq: usize,
    flag_gemcode: bool,
    flag_lr20: bool,
    flag_cr11: bool,
}

/// A Fastq record ready to be written
#[derive(Debug, Serialize, Deserialize, PartialEq, PartialOrd, Eq, Ord)]
struct FqRecord {
    head: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>
}

/// Which read in a pair we have
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialOrd, Ord, PartialEq, Eq)]
enum ReadNum {
    R1,
    R2
}

/// Internally serialized read. Used for re-uniting discordant read pairs
#[derive(Debug, Serialize, Deserialize, PartialOrd, Ord, Eq, PartialEq)]
struct SerFq {
    read_group: Option<Rg>,
    header_key: Vec<u8>,
    rec: FqRecord,
    read_num: ReadNum,
    i1: Option<FqRecord>,
    i2: Option<FqRecord>,
}


struct SerFqShard;

/// Shard SerFq structs based on a hash of the head string
impl ShardDef<SerFq> for SerFqShard {
    fn get_shard(v: &SerFq) -> usize {
        let mut s = SipHasher::new();
        v.header_key.hash(&mut s);
        s.finish() as usize
    }
}

/// Entry in the conversion spec from a BAM record back to a read.
/// Each read can be composed of data from a pair of tags (tag w/ sequence, tag w/ qual),
/// or a fixed-length sequence of Ns (with a default QV), or the sequence in the read.
#[derive(Clone, Debug, PartialEq, Eq)]
enum SpecEntry {
    Tags(String, String),
    Ns(usize),
    Read,
}

type Rg = (String, u32);

/// Spec for converting from a BAM record back to reads. Empty vector indicates that this read doesn't exist
/// in the output. The i1 and i2 reads should be buildable from tags in the R1 read.
#[derive(Clone, Debug)]
struct FormatBamRecords {
    rg_spec: HashMap<String, Rg>,
    r1_spec: Vec<SpecEntry>,
    r2_spec: Vec<SpecEntry>,
    i1_spec: Vec<SpecEntry>,
    i2_spec: Vec<SpecEntry>,
    rename: Option<Vec<String>>,
    order: [u32; 4],
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
        let rgs = Self::parse_rgs(reader);

        if spec.len() == 0 {
            None
        } else {
            Some(
                FormatBamRecords {
                    rg_spec: rgs,
                    r1_spec: spec.remove("R1").unwrap(),
                    r2_spec: spec.remove("R2").unwrap(),
                    i1_spec: spec.remove("I1").unwrap_or_else(|| Vec::new()),
                    i2_spec: spec.remove("I2").unwrap_or_else(|| Vec::new()),
                    rename: None,
                    order: [1,3,2,4],
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
            rename: Some(vec!["R1".to_string(), "R3".to_string(), "I1".to_string(), "R2".to_string()]),
            order: [1,4,2,3]
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
            order: [1,3,2,0]
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
            order: [1,4,3,2],
        }
    }

    fn parse_rgs<R: bam::Read>(reader: &R) -> HashMap<String, Rg> {
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

        if rg_items.len() == 0 {
            println!("WARNING: no @RG (read group) headers found in BAM file. Splitting data by the GEM group marked in the corrected barcode tag.");
            println!("Reads without a corrected barcode will not appear in output FASTQs");
            // No RG items in header -- invent a set fixed set of RGs
            // each observed Gem group in the BAM file will get mapped to these.
            for i in 1 .. 100 {
                let name = format!("gemgroup{:03}", i);
                rg_items.insert(name.clone(), (name, 0));
            }
        }

        rg_items
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
            let rg = match parts.next() {
                Some(v) => v.to_string(),
                None => return None,
            };

            match u32::from_str(lane) {
                Ok(n) => Some((v.clone(), rg, n)),
                Err(_) => {
                    // Handle case in ALIGNER pipeline prior to 2.1.3 -- samtools merge would append a unique identifier to each RG ID tags
                    // Detect this condition and remove from lane
                    let re = Regex::new(r"^([0-9]+)-[0-9A-F]+$").unwrap();
                    let cap = re.captures(lane);


                    if cap.is_none() {
                        None
                    } else {
                        let lane_u32 = u32::from_str(cap.unwrap().get(1).unwrap().as_str()).unwrap();
                        Some((v.clone(), rg, lane_u32))
                    }
                }
            }
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

                    let read = c.get(1).unwrap().as_str().to_string();
                    let tag_list = c.get(2).unwrap().as_str();
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

        spec
    }

    fn try_get_rg(&self, rec: &Record) -> Option<Rg> {
        let rg = rec.aux(b"RG");
        match rg {
            Some(Aux::String(s)) => { 
                let key = String::from_utf8(Vec::from(s)).unwrap();
                self.rg_spec.get(&key).map(|x| x.clone())
            },
            Some(..) => panic!("invalid type of RG header. record: {}", std::str::from_utf8(rec.qname()).unwrap()),
            None => None,
        }
    }

    pub fn find_rg(&self, rec: &Record) -> Option<Rg> {

        let main_rg_tag = self.try_get_rg(rec);

        if main_rg_tag.is_some() {
            return main_rg_tag
        } else {
            let emit = |tag| {
                let corrected_bc = String::from_utf8(Vec::from(tag)).unwrap();
                let mut parts = (&corrected_bc).split("-");
                let _ = parts.next();
                match parts.next() {
                    Some(v) => {
                        match u32::from_str(v) {
                            Ok(v) => {
                                //println!("got gg: {}", v);
                                let name = format!("gemgroup{:03}", v);
                                self.rg_spec.get(&name).map(|x| x.clone())
                            },
                            _ => None
                        }
                    },
                    _ => None,
                }
            };

            // Workaround for early CR 1.1 and 1.2 data
            // Attempt to extract the gem group out of the corrected barcode tag (CB)
            match rec.aux(b"CB") {
                Some(Aux::String(s)) => return emit(s),
                _ => (),
            }

            // Workaround for GemCode (Long Ranger 1.3) data
            // Attempt to extract the gem group out of the corrected barcode tag (BX)
            match rec.aux(b"BX") {
                Some(Aux::String(s)) => return emit(s),
                _ => (),
            }

            None
        }
    }

    /// Convert a BAM record to a Fq record, for internal caching
    pub fn bam_rec_to_ser(&self, rec: &Record) -> SerFq {
        match (rec.is_first_in_template(), rec.is_last_in_template()) {
            (true, false) => {
                SerFq {
                    header_key: rec.qname().to_vec(),
                    read_group: self.find_rg(rec),
                    read_num: ReadNum::R1,
                    rec: self.bam_rec_to_fq(rec, &self.r1_spec, self.order[0]).unwrap(),
                    i1: if self.i1_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i1_spec, self.order[2]).unwrap()) }  else { None },
                    i2: if self.i2_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i2_spec, self.order[3]).unwrap()) }  else { None },
                    
                }
            },
            (false, true) => {
                SerFq {
                    header_key: rec.qname().to_vec(),
                    read_group: self.find_rg(rec),
                    read_num: ReadNum::R2,
                    rec: self.bam_rec_to_fq(rec, &self.r2_spec, self.order[1]).unwrap(),
                    i1: if self.i1_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i1_spec, self.order[2]).unwrap()) }  else { None },
                    i2: if self.i2_spec.len() > 0 { Some(self.bam_rec_to_fq(rec, &self.i2_spec, self.order[3]).unwrap()) }  else { None },
                }
            },
            _ => panic!("Not a valid read pair"),
        }
    }

    /// Convert a BAM record to Fq record ready to be written
    pub fn bam_rec_to_fq(&self, rec: &Record, spec: &Vec<SpecEntry>, read_number: u32) -> Result<FqRecord> {

        let mut head = Vec::new();
        head.extend_from_slice(rec.qname());
        let head_suffix = format!(" {}:N:0:0", read_number);
        head.extend(head_suffix.as_bytes());

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
                        panic!(format!("Invalid BAM record: read: {:?} is missing tag: {:?}", std::str::from_utf8(rec.qname()).unwrap(), read_tag));
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
        let r1 = self.bam_rec_to_fq(r1_rec, &self.r1_spec, self.order[0]).unwrap();
        let r2 = self.bam_rec_to_fq(r2_rec, &self.r2_spec, self.order[1]).unwrap();

        let i1 = if self.i1_spec.len() > 0 {
             Some(self.bam_rec_to_fq(r1_rec, &self.i1_spec, self.order[2]).unwrap())
        } else {
            None
        };

        let i2 = if self.i2_spec.len() > 0 {
             Some(self.bam_rec_to_fq(r1_rec, &self.i2_spec, self.order[3]).unwrap())
        } else {
            None
        };

        let rg = self.find_rg(r1_rec);
        Ok((rg, r1, r2, i1, i2))
    }


    pub fn format_read(&self, rec: &Record) -> Result<(Option<Rg>, FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)> {
        let r1 = self.bam_rec_to_fq(rec, &self.r1_spec, self.order[0]).unwrap();
        let r2 = self.bam_rec_to_fq(rec, &self.r2_spec, self.order[1]).unwrap();

        let i1 = if self.i1_spec.len() > 0 {
             Some(self.bam_rec_to_fq(rec, &self.i1_spec, self.order[2]).unwrap())
        } else {
            None
        };

        let i2 = if self.i2_spec.len() > 0 {
             Some(self.bam_rec_to_fq(rec, &self.i2_spec, self.order[3]).unwrap())
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

type BGW = ThreadProxyWriter<BufWriter<GzEncoder<File>>>;

struct FastqManager {
    writers: HashMap<Rg, FastqWriter>,
}

impl FastqManager {
    pub fn new(out_path: &Path, formatter: FormatBamRecords, sample_name: String, reads_per_fastq: usize) -> FastqManager {

        // Take the read groups and generate read group paths
        let mut sample_def_paths = HashMap::new();
        let mut writers = HashMap::new();

        for (_, &(ref _samp, lane)) in formatter.rg_spec.iter() {
            let samp = _samp.clone();
            let path = sample_def_paths.entry(samp).or_insert_with(|| {
                let suffix = _samp.replace(":", "_");
                let samp_path = out_path.join(suffix);
                //create_dir(&samp_path).expect("couldn't create output directory");
                samp_path
            });
            
            let writer = FastqWriter::new(path, formatter.clone(), "bamtofastq".to_string(), lane, reads_per_fastq);
            writers.insert((_samp.clone(), lane), writer);
        }

        FastqManager { writers: writers }
    }

    pub fn write(&mut self, rg: &Option<Rg>, r1: &FqRecord, r2: &FqRecord, i1: &Option<FqRecord>, i2: &Option<FqRecord>) {
        match rg {
            &Some(ref rg) => {
                 self.writers.get_mut(rg).map(|w| w.write(r1,r2,i1,i2));
            },
            _ => ()
        }
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

    r1: Option<BGW>,
    r2: Option<BGW>,
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
        FastqWriter {
            formatter: formatter,
            out_path: out_path.to_path_buf(),
            sample_name: sample_name,
            lane: lane,
            r1: None,
            r2: None,
            i1: None,
            i2: None,
            n_chunks: 0,
            total_written: 0,
            chunk_written: 0,
            reads_per_fastq: reads_per_fastq,
            path_sets: vec![],
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

            let r1 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, new_read_names[0], n_files+1));
            let r2 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, new_read_names[1], n_files+1));
            let i1 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, new_read_names[2], n_files+1));
            let i2 = out_path.join(format!("{}_S1_L{:03}_{}_{:03}.fastq.gz", sample_name, lane, new_read_names[3], n_files+1));

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

    pub fn try_write_rec2(w: &mut Option<BGW>, rec: &FqRecord) {
        match w {
            &mut Some(ref mut w) => FastqWriter::write_rec(w, rec),
            &mut None => ()
        }
    }

    /// Write a set of fastq records
    pub fn write(&mut self, r1: &FqRecord, r2: &FqRecord, i1: &Option<FqRecord>, i2: &Option<FqRecord>) {

        if self.total_written == 0 {
            // Create the output dir if needed:
            let _ = create_dir(&self.out_path);

            self.cycle_writers()
        }

        FastqWriter::try_write_rec2(&mut self.r1, r1);
        FastqWriter::try_write_rec2(&mut self.r2, r2);
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
        self.r1 = Some(Self::open_gzip_writer(&paths.0));
        self.r2 = Some(Self::open_gzip_writer(&paths.1));
        self.i1 = paths.2.as_ref().map(|p| Self::open_gzip_writer(p));
        self.i2 = paths.3.as_ref().map(|p| Self::open_gzip_writer(p));
        self.n_chunks += 1;
        self.chunk_written = 0;
        self.path_sets.push(paths);
    }

    fn open_gzip_writer<P: AsRef<Path>>(path: P) -> ThreadProxyWriter<BufWriter<GzEncoder<File>>> {
        let f = File::create(path).unwrap();
        let gz = GzEncoder::new(f, Compression::Fast);
        ThreadProxyWriter::new(BufWriter::new(gz), 4096)
    }
}


fn main() {
    if let Err(ref e) = run() {
        println!("error: {}", e);

        for e in e.iter().skip(1) {
            println!("caused by: {}", e);
        }

        // The backtrace is not always generated. Try to run this example
        // with `RUST_BACKTRACE=1`.
        /*
        if let Some(backtrace) = e.backtrace() {
            println!("------------");
            println!("If you believe this is a bug in bamtofastq, please report a bug to software@10xgenomics.com.");
            println!("{:?}", backtrace);
        }
        */

        ::std::process::exit(1);
    }
}

// Most functions will return the `Result` type, imported from the
// `errors` module. It is a typedef of the standard `Result` type
// for which the error type is always our own `Error`.
fn run() -> Result<()> {

    println!("bamtofastq v{}", VERSION);
    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.deserialize())
                         .unwrap_or_else(|e| e.exit());

    let _ = try!(go(args, None));
    Ok(())
}          


pub fn go(args: Args, cache_size: Option<usize>) -> Result<Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>> {

    let cache_size = cache_size.unwrap_or(500000);
 
    match args.flag_locus {
        Some(ref locus) => {
            let loc = try!(locus::Locus::from_str(locus).chain_err(|| "Invalid locus argument. Please use format: 'chr1:123-456'"));
            let mut bam = try!(bam::IndexedReader::from_path(&args.arg_bam).chain_err(|| "Error opening BAM file. The BAM file must be indexed when using --locus"));
            let tid = try!(bam.header().tid(loc.chrom.as_bytes()).ok_or(format!("Requested chromosome not present: {}", loc.chrom)));
            try!(bam.fetch(tid, loc.start, loc.end));
            inner(args.clone(), cache_size, bam)
        },
        None => {
            let bam = try!(bam::Reader::from_path(&args.arg_bam).chain_err(|| "Error opening BAM file"));
            inner(args, cache_size, bam)
        }
    }
}

pub fn inner<R: bam::Read>(args: Args, cache_size: usize, bam: R) -> Result<Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>> {

    let formatter = {
        let header_fmt = FormatBamRecords::from_headers(&bam);
        match header_fmt {
            Some(mut f) => {
                // HACK -- special case
                // we need a special case for SC 3' v1 data
                // because of the bc-in-index setup, it needs the rename field set,
                // even though the BAM headers support in theory tell us what to do.
                // detect this case here and set the right rename field.
                if f.r1_spec == vec![SpecEntry::Read] && 
                   f.i1_spec == vec![SpecEntry::Tags("CR".to_string(), "CY".to_string())] &&
                   f.i2_spec == vec![SpecEntry::Tags("BC".to_string(), "QT".to_string())] {

                    f.rename = Some(vec!["R1".to_string(), "R3".to_string(), "R2".to_string(), "I1".to_string()])
                }

                if args.flag_gemcode {
                    return Err("Do not use a pipeline-specific command-line flag: --gemcode. Supplied BAM file already contains bamtofastq headers.".into())
                }

                if args.flag_lr20 {
                    return Err("Do not use a pipeline-specific command-line flag: --lr20. Supplied BAM file already contains bamtofastq headers.".into())
                }

                if args.flag_cr11 {
                    return Err("Do not use a pipeline-specific command-line flag: --cr11. Supplied BAM file already contains bamtofastq headers.".into())
                }                


                f
            },
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
                    return Ok(vec![]);
                }
            }
        }
    };

    let out_path = Path::new(&args.arg_output_path);
    match create_dir(&out_path) {
        Err(e) => return Err(format!("error creating output directory: {:?}. Does it already exist?", &out_path).into()),
        _ => ()
    }
    
    // prep output files
    println!("{:?}", args);
    let fq = FastqManager::new(out_path, formatter.clone(), "bamtofastq".to_string(), args.flag_reads_per_fastq);
 
    if formatter.is_double_ended() {
        if args.flag_bx_list.is_some() {

            // BX-sorted case: get a selected set of BXs
            let bxi = try!(bx_index::BxIndex::new(args.arg_bam));
            let bx_iter = try!(BxListIter::from_path(args.flag_bx_list.unwrap(), bxi, bam));
            proc_double_ended(bx_iter, formatter, fq, cache_size, false)
        } else {

            // Standard pos-sorted case
            proc_double_ended(bam.records(), formatter, fq, cache_size, args.flag_locus.is_some())
        }
    } else {
        if args.flag_bx_list.is_some() {
            // BX-sorted case: get a selected set of BXs
            let bxi = try!(bx_index::BxIndex::new(args.arg_bam));
            let bx_iter = try!(BxListIter::from_path(args.flag_bx_list.unwrap(), bxi, bam));
            proc_double_ended(bx_iter, formatter, fq, cache_size, false)
        } else {
            proc_single_ended(bam.records(), formatter, fq)
        }
    } 
}


fn proc_double_ended<I>(records: I, formatter: FormatBamRecords, mut fq: FastqManager, cache_size: usize, restricted_locus: bool) -> 
Result<Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>> 
where I: Iterator<Item=result::Result<Record, ReadError>> {
    // Temp file for hold unpaired reads. Will be cleaned up automatically.
    let tmp_file = try!(NamedTempFile::new());

    let total_read_pairs = {
        // Cache for efficiently finding local read pairs
        let mut rp_cache = RpCache::new(cache_size);

        // For chimeric read piars that are showing up in different places, we will write these to disk for later use
        let w: ShardWriteManager<SerFq, SerFqShard> = ShardWriteManager::new(tmp_file.path(), 2048, 256, 2);
        let mut sender = w.get_sender();

        // Count total R1s observed, so we can make sure we've preserved all read pairs
        let mut total_read_pairs = 0;

        for _rec in records {
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
    let reader = ShardReader::<SerFq>::open(tmp_file.path());

    let mut ncached = 0;
    for s in 0..reader.num_shards() {
        let mut data = reader.read_shard(s);
        data.sort();

        for (_, items) in &data.iter().group_by(|x| &x.header_key) {
            // write out items
            let mut item_vec: Vec<_> = items.collect();

            // We're missing a read in the pair, and we would expect it.
            if item_vec.len() != 2 && !restricted_locus {
                panic!("didn't get both reads!: {:?}", item_vec);
            }

            // We're missing a read in the pair, and we would expect it.
            if item_vec.len() != 2 && restricted_locus {
                continue
            }

            item_vec.sort_by_key(|x| x.read_num);
            let r1 = item_vec.swap_remove(0);
            let r2 = item_vec.swap_remove(0);
            fq.write(&r1.read_group, &r1.rec, &r2.rec, &r1.i1, &r1.i2);
            ncached += 1;
        }
    }

    // make sure we have the right number of output reads
    println!("Writing finished.  Observed {} unique read ids. Wrote {} read pairs ({} cached)", total_read_pairs, fq.total_written(), ncached);
    Ok(fq.paths())
}

fn proc_single_ended<I>(records: I, formatter: FormatBamRecords, mut fq: FastqManager) ->
  Result<Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>> 
  where I: Iterator<Item=result::Result<Record, ReadError>> {

    let total_reads = {
        // Count total R1s observed, so we can make sure we've preserved all read pairs
        let mut total_reads = 0;

        for _rec in records {
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
    Ok(fq.paths())
}


#[cfg(test)]
mod tests {
    use tempdir;
    use super::*;
    use fastq_reader::*;
    use std::collections::HashMap;

    type ReadSet = HashMap<Vec<u8>, RawReadSet>;

    fn strip_extra_headers(header: &Vec<u8>) -> Vec<u8> {
        let head_str = String::from_utf8(header.clone()).unwrap();
        let mut split = head_str.split_whitespace();
        split.next().unwrap().to_string().into_bytes()
    }

    fn strip_header_fqrec(r: FqRec) -> FqRec {
        (strip_extra_headers(&(r.0)), r.1, r.2)
    }

    fn strip_header_raw_read_set(r: RawReadSet) -> RawReadSet {
        (strip_header_fqrec(r.0), strip_header_fqrec(r.1), r.2.map(|x| strip_header_fqrec(x)))
    }

    // Load fastqs, but strip extra elements of the FASTQ header beyond the first space -- they will not be in the BAM
    pub fn load_fastq_set<I: Iterator<Item=RawReadSet>>(reads: &mut ReadSet, iter: I) {
        for r in iter {
            reads.insert(strip_extra_headers(&((r.0).0)), strip_header_raw_read_set(r));
        }
    }


    pub fn strict_compare_read_sets(orig_set: ReadSet, new_set: ReadSet) {
        
        assert_eq!(orig_set.len(), new_set.len());

        let mut keys1: Vec<Vec<u8>> = orig_set.keys().cloned().collect();
        keys1.sort();

        let mut keys2: Vec<Vec<u8>> = new_set.keys().cloned().collect();
        keys2.sort();

        for (k1, k2) in keys1.iter().zip(keys2.iter()) {
            assert_eq!(k1, k2);
            assert_eq!(orig_set.get(k1), new_set.get(k2));
        }
    }

    pub fn subset_compare_read_sets(orig_set: ReadSet, new_set: ReadSet) {

        assert!(orig_set.len() > new_set.len());

        for k in new_set.keys() {
            assert_eq!(new_set.get(k), orig_set.get(k))
        }
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
            flag_bx_list: None,
        };

        let out_path_sets = super::go(args, Some(2)).unwrap();

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/crg-tiny-fastq-2.0.0/read-RA_si-GTTGCAGC_lane-001-chunk-001.fastq.gz", 
            Some("test/crg-tiny-fastq-2.0.0/read-I1_si-GTTGCAGC_lane-001-chunk-001.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, true_fastq_read);

        let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1, r2, i1));
        }
        
        strict_compare_read_sets(orig_reads, output_reads);
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
            flag_locus: None,
            flag_bx_list: None,
        };

        let out_path_sets = super::go(args, Some(2)).unwrap();

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
            flag_bx_list: None,
        };

        let out_path_sets = super::go(args, Some(2)).unwrap();

        let true_fastq_read = open_interleaved_fastq_pair_iter(
            "test/cellranger-tiny-fastq-1.2.0/read-RA_si-TTTCATGA_lane-008-chunk-001.fastq.gz",
            Some("test/cellranger-tiny-fastq-1.2.0/read-I1_si-TTTCATGA_lane-008-chunk-001.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, true_fastq_read);

       let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1, r2, i1));
        }
        
        subset_compare_read_sets(orig_reads, output_reads);
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
            flag_locus: None,
            flag_bx_list: None,
        };

        let out_path_sets = super::go(args, Some(2)).unwrap();

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
        
        subset_compare_read_sets(orig_reads, output_reads);

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

        subset_compare_read_sets(orig_index_reads, output_index_reads);
    }
}
