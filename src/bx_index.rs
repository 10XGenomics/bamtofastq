// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.

extern crate csv;

use std::path::{Path, PathBuf};
use std::fs::File;
use rust_htslib::bam::{self, Read, ReadError};
use rust_htslib::bam::record::{Aux, Record};
use std::result;

use errors::*;

#[derive(Deserialize, Ord, PartialOrd, Eq, PartialEq)]
struct BcObs {
    bx: String,
    offset: usize,
}

pub struct BxIndex {
    obs: Vec<BcObs>
}

impl BxIndex {
    pub fn new<P: AsRef<Path>>(bam_file: P) -> Result<BxIndex> {
        let bxi_fn = bam_file.as_ref().with_extension("bam.bxi");


        let f = try!(File::open(bxi_fn.clone()).
            chain_err(|| 
                format!("Couldn't find BX index: '{:?}'. You must sort you BAM file with 'samtools sort -t BX' and index with 'bxindex'", bxi_fn)
            ));
        let mut reader = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(f);
        let obs: Vec<BcObs> = reader.deserialize().map(|x| x.unwrap()).collect();

        let res = BxIndex { obs: obs };
        Ok(res)
    }

    pub fn get_voffset(&self, bx: &String) -> usize
    {
        match self.obs.binary_search_by_key(&bx, |x| &x.bx) {
            Ok(pos) => {
                self.obs[pos].offset
            },
            Err(pos) => {
                self.obs[pos-1].offset
            }
        }
    }

    pub fn get_bx_list<P: AsRef<Path>>(bx_list_file: P) -> Result<Vec<String>> {
        let f = try!(File::open(bx_list_file));
        let mut reader = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(f);
        Ok(reader.deserialize().map(|x| x.unwrap()).collect())
    }
}

pub fn get_records_for_bx<R: Read>(index: &BxIndex, reader: &mut R, bx: &String) -> Vec<Record> {

    let start = index.get_voffset(bx);
    reader.seek(start as i64);

    let mut rec_iter = reader.records();
    let mut recs = Vec::new();

    loop {
        let rec = match rec_iter.next() {
            Some(Ok(r)) => r,
            _ => break,
        };

        let bx_read = 
            match rec.aux(b"BX") {
                Some(Aux::String(s)) => {
                    String::from_utf8(Vec::from(s)).unwrap()
                },
                _ => "".to_string(),
            };

        if &bx_read == bx {
            recs.push(rec);
        } else if &bx_read < bx {
            continue;
        } else {
            break;
        }
    }

    recs
}

//pub fn go(index: &BxIndex, reader: &mut R, bx: &Vec<String>) -> 


pub struct BxListIter<R: Read> {
    index: BxIndex,
    reader: R,
    bx_list: Vec<String>,

    cur_bx: usize,
    cur_vec: Vec<Record>
}

impl<R: Read> BxListIter<R> {
    pub fn new(bx_list: Vec<String>, index: BxIndex, mut reader: R) -> BxListIter<R> {
        let cur_vec =
            if bx_list.len() > 0 { 
                let mut v = get_records_for_bx(&index, &mut reader, &bx_list[0]);
                v.reverse();
                v
            } else { 
                vec![]
            };

        BxListIter {
            index: index,
            reader: reader,
            bx_list: bx_list,
            cur_bx: 0,
            cur_vec: cur_vec
        }
    }

    pub fn from_path(bx_list: String, index: BxIndex, reader: R) -> Result<BxListIter<R>> {
        let list = try!(BxIndex::get_bx_list(PathBuf::from(bx_list)));
        Ok(BxListIter::new(list, index, reader))
    }
}

impl<R: Read> Iterator for BxListIter<R> {
    type Item = result::Result<Record, ReadError>;

    fn next(&mut self) -> Option<result::Result<Record, ReadError>> {

        while self.cur_bx <= self.bx_list.len() {

            if self.cur_vec.len() == 0 {
                self.cur_bx += 1;
                if self.cur_bx == self.bx_list.len() {
                    break;
                }

                self.cur_vec = get_records_for_bx(&self.index, &mut self.reader, &self.bx_list[self.cur_bx]);
                self.cur_vec.reverse();
            }

            match self.cur_vec.pop() {
                Some(v) => return Some(Ok(v)),
                None => (),
            }
        }

        return None
    }
}