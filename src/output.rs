use std::io::Write;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::compute::InnerDistanceResult;

/// Write `<prefix>.inner_distance.txt` (per-pair) and
/// `<prefix>.inner_distance_freq.txt` (histogram).
pub fn write_output(result: &InnerDistanceResult, prefix: &str) -> Result<()> {
    let txt_path = format!("{prefix}.inner_distance.txt");
    let mut fo = std::fs::File::create(&txt_path).map_err(|e| {
        RsomicsError::Io(std::io::Error::other(format!("creating {txt_path}: {e}")))
    })?;
    for rec in &result.pairs {
        match rec.inner_distance {
            None => writeln!(fo, "{}\tNA\t{}", rec.read_name, rec.category.as_str()),
            Some(d) => writeln!(fo, "{}\t{}\t{}", rec.read_name, d, rec.category.as_str()),
        }
        .map_err(RsomicsError::Io)?;
    }

    let freq_path = format!("{prefix}.inner_distance_freq.txt");
    let mut fq = std::fs::File::create(&freq_path).map_err(|e| {
        RsomicsError::Io(std::io::Error::other(format!("creating {freq_path}: {e}")))
    })?;
    for &(st, end, count) in &result.histogram {
        writeln!(fq, "{st}\t{end}\t{count}").map_err(RsomicsError::Io)?;
    }

    Ok(())
}

/// JSON-serialisable summary (emitted on `--json`).
#[derive(Debug, Serialize)]
pub struct InnerDistanceSummary {
    pub pair_num: u64,
    pub output_prefix: String,
}
