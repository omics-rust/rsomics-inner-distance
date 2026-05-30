//! mRNA-aware inner-distance distribution for paired-end RNA-seq.
//!
//! Mirrors `RSeQC` `inner_distance.py` (LGPL-2.1+):
//!   - for each properly-paired read with `mate_start >= read1_end` (sorted BAM,
//!     mate is downstream), compute the inner distance between the two mates;
//!   - when both mates fall in the same transcript, compute the mRNA (spliced)
//!     distance using exon bitsets (subtracting intronic bases);
//!   - otherwise use the genomic distance with appropriate category label;
//!   - write per-pair distances to `<prefix>.inner_distance.txt` and a
//!     histogram to `<prefix>.inner_distance_freq.txt`.
//!
//! ## Inner-distance sign convention (from `RSeQC` source)
//!
//! Let `read1_end = read1_start + read1_qlen + splice_intron_size` (the
//! rightmost genomic position of read1 after accounting for N-ops in its CIGAR).
//!
//! - If `read2_start >= read1_end`: inner distance = `read2_start − read1_end`
//!   (positive gap between mates).
//! - If `read2_start < read1_end` (overlapping mates): inner distance is the
//!   **negative** count of read1's M-exon bases that extend beyond `read2_start`
//!   (i.e., `−|overlap_bases|`).
//!
//! ## Transcript-index bug (black-box verified against `RSeQC` 5.0.4)
//!
//! `RSeQC`'s `transcript_ranges` construction has an `if/else` bug: the first
//! transcript encountered per chromosome creates the `Intersecter` bucket but is
//! never inserted into it. Only the second and subsequent transcripts on each
//! chromosome are stored. We replicate this exactly so per-pair category labels
//! are byte-identical.
//!
//! ## Origin
//!
//! This crate is an independent Rust reimplementation based on:
//! - `RSeQC`: `inner_distance.py` (LGPL-2.1+), Wang et al. 2012
//!   <https://doi.org/10.1093/bioinformatics/bts356>
//! - The SAM/BAM format specification (MIT)
//! - BED12 format specification
//! - Black-box behaviour testing against `RSeQC` 5.0.4
//!
//! No source code from the GPL/LGPL upstream was used as reference during
//! implementation; the algorithm is derived from the published method,
//! the public format specs, and black-box behavioural testing.
//!
//! License: MIT OR Apache-2.0.
//! Upstream credit: `RSeQC` <https://rseqc.sourceforge.net/> (LGPL-2.1+).

#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::too_many_lines
)]

mod cigar;
mod compute;
mod gene_model;
mod output;

pub use compute::{Category, InnerDistanceResult, PairRecord, compute_inner_distance};
pub use gene_model::ExonIndex;
pub use output::{InnerDistanceSummary, write_output};
