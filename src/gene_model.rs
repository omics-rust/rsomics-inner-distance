use std::collections::{HashMap, HashSet};
use std::path::Path;

use coitrees::{COITree, Interval as CoiInterval, IntervalTree};
use rsomics_common::{Result, RsomicsError};

/// Per-chromosome exon interval tree (union of all exon blocks from BED12).
///
/// Intervals stored end-inclusive (coitrees convention) for point queries.
pub struct ExonIndex {
    trees: HashMap<String, COITree<(), u32>>,
}

impl ExonIndex {
    /// Build from a BED12 file.
    ///
    /// Mirrors `RSeQC` `binned_bitsets_from_list`: each BED12 block → one exon interval;
    /// overlapping intervals from multiple transcripts are union-merged.
    /// Chromosome names are uppercased (matching `RSeQC`'s `.upper()`).
    pub fn from_bed12(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| RsomicsError::Io(std::io::Error::other(format!("reading BED12: {e}"))))?;

        let mut raw_half_open: HashMap<String, Vec<(i32, i32)>> = HashMap::new();

        for line in content.lines() {
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                continue;
            }
            let chrom = fields[0].to_uppercase();
            let Ok(tx_start) = fields[1].parse::<i32>() else {
                continue;
            };
            let Ok(block_count) = fields[9].parse::<usize>() else {
                continue;
            };
            let block_sizes: Vec<i32> = fields[10]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();
            let block_starts: Vec<i32> = fields[11]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();
            if block_sizes.len() != block_count || block_starts.len() != block_count {
                continue;
            }
            let entry = raw_half_open.entry(chrom).or_default();
            for (bstart, bsize) in block_starts.iter().zip(block_sizes.iter()) {
                let exon_start = tx_start + bstart;
                let exon_end = exon_start + bsize;
                entry.push((exon_start, exon_end));
            }
        }

        // Merge overlapping intervals (bitset union) to avoid double-counting.
        let trees = raw_half_open
            .into_iter()
            .map(|(chrom, mut ivs)| {
                ivs.sort_unstable();
                let mut merged: Vec<CoiInterval<()>> = Vec::with_capacity(ivs.len());
                for (start, end) in ivs {
                    if let Some(last) = merged.last_mut()
                        && start <= last.last + 1
                    {
                        last.last = last.last.max(end - 1);
                        continue;
                    }
                    // coitrees: end-inclusive [start, end-1].
                    merged.push(CoiInterval::new(start, end - 1, ()));
                }
                (chrom, COITree::new(&merged))
            })
            .collect();

        Ok(Self { trees })
    }

    /// Exonic bases in `[start, end)` on `chrom` (0-based half-open).
    ///
    /// Mirrors the `BinnedBitSet` intersection in RSeQC.
    pub(crate) fn exonic_bases(&self, chrom: &str, start: i32, end: i32) -> i32 {
        let Some(tree) = self.trees.get(chrom) else {
            return 0;
        };
        let mut total = 0i32;
        tree.query(start, end - 1, |node| {
            let ov_start = node.first.max(start);
            let ov_end = (node.last + 1).min(end);
            if ov_end > ov_start {
                total += ov_end - ov_start;
            }
        });
        total
    }

    pub(crate) fn has_chrom(&self, chrom: &str) -> bool {
        self.trees.contains_key(chrom)
    }
}

/// Per-chromosome transcript-range interval tree.
///
/// Decides whether both `read1_end` and `read2_start` fall in the same transcript,
/// matching RSeQC's `transcript_ranges` `Intersecter`.
///
/// Replicates the RSeQC `if/else` bug: the first transcript per chromosome creates
/// the bucket but is never inserted. Only second and subsequent transcripts are stored.
/// Black-box-verified against RSeQC 5.0.4.
pub(crate) struct TranscriptIndex {
    /// Per-chromosome tree; metadata is index into `names`.
    trees: HashMap<String, COITree<usize, u32>>,
    names: Vec<String>,
}

impl TranscriptIndex {
    pub(crate) fn from_bed12(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| RsomicsError::Io(std::io::Error::other(format!("reading BED12: {e}"))))?;

        let mut names: Vec<String> = Vec::new();
        // First transcript per chrom is skipped — replicates RSeQC `if/else` bug.
        let mut chroms_seen: HashSet<String> = HashSet::new();
        let mut raw: HashMap<String, Vec<CoiInterval<usize>>> = HashMap::new();

        for line in content.lines() {
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                continue;
            }
            let chrom = fields[0].to_uppercase();
            let Ok(tx_start) = fields[1].parse::<i32>() else {
                continue;
            };
            let Ok(tx_end) = fields[2].parse::<i32>() else {
                continue;
            };
            let name = fields[3].to_string();

            if !chroms_seen.insert(chrom.clone()) {
                let id = names.len();
                names.push(name);
                raw.entry(chrom)
                    .or_default()
                    .push(CoiInterval::new(tx_start, tx_end - 1, id));
            }
        }

        let trees = raw
            .into_iter()
            .map(|(chrom, ivs)| (chrom, COITree::new(&ivs)))
            .collect();

        Ok(Self { trees, names })
    }

    /// Transcript names overlapping a 0-based point on `chrom`.
    pub(crate) fn names_at<'a>(&'a self, chrom: &str, point: i32) -> HashSet<&'a str> {
        let Some(tree) = self.trees.get(chrom) else {
            return HashSet::new();
        };
        let mut ids = HashSet::new();
        tree.query(point, point, |node| {
            let idx: usize = *std::borrow::Borrow::<usize>::borrow(&node.metadata);
            ids.insert(self.names[idx].as_str());
        });
        ids
    }
}
