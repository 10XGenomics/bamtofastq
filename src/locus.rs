// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.

use std::str::FromStr;
use regex::Regex;
use std::fmt;

error_chain! {
    errors {
        LocusParseError {
            description("invalid locus string")
        }
    }
}


#[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Deserialize, Clone)]
pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

impl Locus {
    pub fn from_string(s: &str) -> Locus {
        FromStr::from_str(s).ok().unwrap()
    }
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"{}:{}-{}", self.chrom, self.start, self.end)
    }
}




fn remove_commas(s: &str) -> String {
    let ss = s.to_string();
    ss.replace(",", "")
}


impl FromStr for Locus {
    type Err = Error;

    fn from_str(s: &str) -> Result<Locus> {
        let re = Regex::new(r"^(.*):([0-9,]+)(-|..)([0-9,]+)$").unwrap();
        let cap = re.captures(s);

        if cap.is_none() {
            return Err(ErrorKind::LocusParseError.into())
        }

        let cap = cap.unwrap();

        let start_s = remove_commas(cap.get(2).unwrap().as_str());
        let end_s = remove_commas(cap.get(4).unwrap().as_str());

        Ok(Locus {
            chrom: cap.get(1).unwrap().as_str().to_string(),
            start: FromStr::from_str(&start_s).unwrap(),
            end: FromStr::from_str(&end_s).unwrap(),
        })
    }
}
