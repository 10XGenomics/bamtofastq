// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.

use std::collections::HashMap;
use rust_htslib::bam::record::{Aux, Record};

/// Read-pair cache. Let's us stream through the BAM and find nearby mates so we can write them out immediately
/// Reads whose mate is not found promptly are written to disk and matched up later
pub struct RpCache {
    pub cache_size: usize,
    pub cache: HashMap<Vec<u8>, Record>
}

impl RpCache {

    pub fn new(cache_size: usize) -> RpCache {
        RpCache { cache: HashMap::new(), cache_size }
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
                    println!("Found invalid set of BAM record for qname: {}.", String::from_utf8_lossy(rec.qname()));
                    println!("This may be caused by inputting the same FASTQ record to Long Ranger twice");
                    panic!("invalid BAM record detected: {}", String::from_utf8_lossy(rec.qname()));
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

        let mut dist = 5000;
       

        while self.cache.len() > self.cache_size / 2 && dist > 100{
            let mut orphan_keys = Vec::new();

            for (key, rec) in self.cache.iter() {
                // Evict unmapped reads, reads on a previous chromosome, or reads that are >5kb behind the current position
                if rec.tid() == -1 || (current_pos - rec.pos()).abs() > dist || rec.tid() != current_tid {
                    orphan_keys.push(key.clone());
                }
            }

            for k in orphan_keys {
                let rec = self.cache.remove(&k).unwrap();
                orphans.push(rec);
            }

            dist = dist / 2;
        }

        // Cache got too full -- just clear it
        if dist <= 100 {
            for (_, rec) in self.cache.drain() {
                orphans.push(rec)
            }

            self.cache = HashMap::new();
        }

        orphans
    }

    pub fn len(&self) -> usize 
    {
        self.cache.len()
    }
}
