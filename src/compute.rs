use std::num::NonZero;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

use crate::cigar::{cigar_intron_len, cigar_qlen, count_overlap_exon_bases};
use crate::gene_model::{ExonIndex, TranscriptIndex};

// BAM flag bits (SAM spec §1.4).
const FLAG_PAIRED: u16 = 0x0001;
const FLAG_QCFAIL: u16 = 0x0200;
const FLAG_DUPLICATE: u16 = 0x0400;
const FLAG_SECONDARY: u16 = 0x0100;
const FLAG_UNMAPPED: u16 = 0x0004;
const FLAG_MATE_UNMAPPED: u16 = 0x0008;
const FLAG_READ1: u16 = 0x0040;

/// Category string written to `inner_distance.txt` per read pair.
#[derive(Debug, Clone)]
pub enum Category {
    SameChromNo,
    SameTranscriptNo,
    SameTranscriptYesSameExonYes,
    SameTranscriptYesSameExonNo,
    SameTranscriptYesNonExonic,
    UnknownChromosome,
    ReadPairOverlap,
}

impl Category {
    pub(crate) fn as_str(&self) -> &'static str {
        match self {
            Category::SameChromNo => "sameChrom=No",
            Category::SameTranscriptNo => "sameTranscript=No,dist=genomic",
            Category::SameTranscriptYesSameExonYes => "sameTranscript=Yes,sameExon=Yes,dist=mRNA",
            Category::SameTranscriptYesSameExonNo => "sameTranscript=Yes,sameExon=No,dist=mRNA",
            Category::SameTranscriptYesNonExonic => "sameTranscript=Yes,nonExonic=Yes,dist=genomic",
            Category::UnknownChromosome => "unknownChromosome,dist=genomic",
            Category::ReadPairOverlap => "readPairOverlap",
        }
    }
}

/// One per-read-pair result row.
pub struct PairRecord {
    pub read_name: String,
    /// `None` if mates on different chromosomes.
    pub inner_distance: Option<i32>,
    pub category: Category,
}

/// Results of inner-distance computation.
pub struct InnerDistanceResult {
    pub pairs: Vec<PairRecord>,
    /// Histogram bins: `(bin_start, bin_end_exclusive, count)`.
    pub histogram: Vec<(i32, i32, u64)>,
    pub pair_num: u64,
}

/// Run the inner-distance computation.
#[allow(clippy::too_many_arguments)]
pub fn compute_inner_distance(
    bam_path: &std::path::Path,
    bed_path: &std::path::Path,
    sample_size: u64,
    mapq_cut: u8,
    low_bound: i32,
    up_bound: i32,
    step: i32,
) -> Result<InnerDistanceResult> {
    let exon_index = ExonIndex::from_bed12(bed_path)?;
    let tx_index = TranscriptIndex::from_bed12(bed_path)?;

    let mut pairs: Vec<PairRecord> = Vec::new();
    let mut pair_num: u64 = 0;

    #[allow(clippy::cast_sign_loss)]
    let window_left_bounds: Vec<i32> = (low_bound..up_bound).step_by(step as usize).collect();
    let mut distances: Vec<i32> = Vec::new();

    let workers = NonZero::new(1usize).unwrap();
    let mut reader = rsomics_bamio::open_with_workers(bam_path, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    let mut rec = RawRecord::default();
    loop {
        if pair_num >= sample_size {
            break;
        }
        let bytes_read = raw::read_record(reader.get_mut(), &mut rec)?;
        if bytes_read == 0 {
            break;
        }

        let flags = rec.flags();
        // RSeQC filter order (from source).
        if flags & FLAG_QCFAIL != 0 {
            continue;
        }
        if flags & FLAG_DUPLICATE != 0 {
            continue;
        }
        if flags & FLAG_SECONDARY != 0 {
            continue;
        }
        if flags & FLAG_UNMAPPED != 0 {
            continue;
        }
        if flags & FLAG_PAIRED == 0 {
            continue;
        }
        if flags & FLAG_MATE_UNMAPPED != 0 {
            continue;
        }
        if rec.mapping_quality() < mapq_cut {
            continue;
        }

        let read1_start = rec.alignment_start();
        let mate_start = rec.mate_alignment_start();

        // Skip mate-upstream reads (BAM sorted; mate already processed).
        if mate_start < read1_start {
            continue;
        }
        // Same position and read1 → inner_distance = 0; RSeQC skips output.
        if mate_start == read1_start && (flags & FLAG_READ1 != 0) {
            continue;
        }

        pair_num += 1;

        let read_name = String::from_utf8_lossy(rec.name())
            .trim_end_matches('\0')
            .to_string();

        let tid = rec.reference_sequence_id();
        let rnext = rec.mate_reference_sequence_id();
        if tid != rnext || tid < 0 {
            pairs.push(PairRecord {
                read_name,
                inner_distance: None,
                category: Category::SameChromNo,
            });
            continue;
        }

        #[allow(clippy::cast_sign_loss)]
        let chrom = {
            let idx = tid as usize;
            let Some(name) = ref_names.get(idx) else {
                pairs.push(PairRecord {
                    read_name,
                    inner_distance: Some(0),
                    category: Category::UnknownChromosome,
                });
                continue;
            };
            name.to_uppercase()
        };

        let cigar_vec: Vec<(u8, u32)> = rec.cigar_ops().collect();
        let read1_len = cigar_qlen(cigar_vec.iter().copied());
        let splice_intron_size = cigar_intron_len(cigar_vec.iter().copied());

        // read1_end = read1_start + qlen + intron (RSeQC sign convention).
        let read1_end = read1_start + read1_len + splice_intron_size;
        let read2_start = mate_start;

        let category: Category;
        let written_distance: i32;

        if read2_start >= read1_end {
            let geo_distance = read2_start - read1_end;

            let read1_end_genes = tx_index.names_at(&chrom, read1_end - 1);
            let read2_start_genes = tx_index.names_at(&chrom, read2_start);
            let same_tx = read1_end_genes
                .intersection(&read2_start_genes)
                .next()
                .is_some();

            if !same_tx {
                written_distance = geo_distance;
                category = Category::SameTranscriptNo;
                distances.push(geo_distance);
            } else if exon_index.has_chrom(&chrom) {
                let exon_bases = exon_index.exonic_bases(&chrom, read1_end, read2_start);
                if exon_bases == geo_distance {
                    written_distance = exon_bases;
                    category = Category::SameTranscriptYesSameExonYes;
                    distances.push(exon_bases);
                } else if exon_bases > 0 && exon_bases < geo_distance {
                    written_distance = exon_bases;
                    category = Category::SameTranscriptYesSameExonNo;
                    distances.push(exon_bases);
                } else {
                    written_distance = geo_distance;
                    category = Category::SameTranscriptYesNonExonic;
                    distances.push(geo_distance);
                }
            } else {
                written_distance = geo_distance;
                category = Category::UnknownChromosome;
                distances.push(geo_distance);
            }
        } else {
            // Overlapping mates: negative inner distance.
            let overlap_dist = -count_overlap_exon_bases(
                read1_start,
                read1_end,
                read2_start,
                cigar_vec.iter().copied(),
            );
            written_distance = overlap_dist;

            let read1_end_genes = tx_index.names_at(&chrom, read1_end - 1);
            let read2_start_genes = tx_index.names_at(&chrom, read2_start);
            let same_tx = read1_end_genes
                .intersection(&read2_start_genes)
                .next()
                .is_some();

            if same_tx {
                category = Category::ReadPairOverlap;
            } else {
                category = Category::SameTranscriptNo;
            }
            distances.push(overlap_dist);
        }

        pairs.push(PairRecord {
            read_name,
            inner_distance: Some(written_distance),
            category,
        });
    }

    // RSeQC histogram: Interval(d-1, d) per distance d; bin [st, st+step) counts d where
    // st < d <= st+step (bx-python find semantics on half-open intervals).
    let histogram: Vec<(i32, i32, u64)> = window_left_bounds
        .iter()
        .map(|&st| {
            let count = distances
                .iter()
                .filter(|&&d| d > st && d <= st + step)
                .count() as u64;
            (st, st + step, count)
        })
        .collect();

    Ok(InnerDistanceResult {
        pairs,
        histogram,
        pair_num,
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn histogram_bin_boundary() {
        // d=50 falls in bin [45,50) since condition is d > st AND d <= st+step
        let distances = [50i32];
        let st = 45i32;
        let step = 5i32;
        let count = distances
            .iter()
            .filter(|&&d| d > st && d <= st + step)
            .count();
        assert_eq!(count, 1);

        // d=50 must NOT fall in bin [50,55): 50 > 50 is false
        let st2 = 50i32;
        let count2 = distances
            .iter()
            .filter(|&&d| d > st2 && d <= st2 + step)
            .count();
        assert_eq!(count2, 0);
    }
}
