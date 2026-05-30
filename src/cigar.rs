// BAM CIGAR opcodes (SAM spec §1.4.6, 0-indexed).
pub(crate) const CIGAR_MATCH: u8 = 0;
pub(crate) const CIGAR_INS: u8 = 1;
pub(crate) const CIGAR_DEL: u8 = 2;
pub(crate) const CIGAR_REF_SKIP: u8 = 3;
pub(crate) const CIGAR_SOFT_CLIP: u8 = 4;
pub(crate) const CIGAR_SEQ_MATCH: u8 = 7;
pub(crate) const CIGAR_SEQ_MISMATCH: u8 = 8;

/// Query sequence length — bases that consume query (M/I/S/=/X).
///
/// Mirrors `pysam` `aligned_read.qlen`.
pub(crate) fn cigar_qlen(cigar: impl Iterator<Item = (u8, u32)>) -> i32 {
    let mut q = 0i32;
    for (op, len) in cigar {
        match op {
            CIGAR_MATCH | CIGAR_INS | CIGAR_SOFT_CLIP | CIGAR_SEQ_MATCH | CIGAR_SEQ_MISMATCH => {
                q += len as i32;
            }
            _ => {}
        }
    }
    q
}

/// Total intron (N-op) length from a CIGAR.
///
/// Mirrors `bam_cigar.fetch_intron` summing `intron[2] - intron[1]`.
pub(crate) fn cigar_intron_len(cigar: impl Iterator<Item = (u8, u32)>) -> i32 {
    cigar
        .filter(|&(op, _)| op == CIGAR_REF_SKIP)
        .map(|(_, len)| len as i32)
        .sum()
}

/// M-bases of read1 that extend beyond `read2_start` (0-based, half-open).
///
/// RSeQC counts 1-based positions from M blocks where `i > read2_start AND i <= read1_end`.
/// Translated to 0-based half-open: M-bases in `[read2_start, read1_end)`.
pub(crate) fn count_overlap_exon_bases(
    read1_start: i32,
    read1_end: i32,
    read2_start: i32,
    cigar: impl Iterator<Item = (u8, u32)>,
) -> i32 {
    let mut ref_pos = read1_start;
    let mut count = 0i32;
    for (op, len) in cigar {
        #[allow(clippy::cast_possible_wrap)]
        let l = len as i32;
        match op {
            CIGAR_MATCH | CIGAR_SEQ_MATCH | CIGAR_SEQ_MISMATCH => {
                let ov_start = ref_pos.max(read2_start);
                let ov_end = (ref_pos + l).min(read1_end);
                if ov_end > ov_start {
                    count += ov_end - ov_start;
                }
                ref_pos += l;
            }
            CIGAR_DEL | CIGAR_REF_SKIP => {
                ref_pos += l;
            }
            _ => {}
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_qlen_basic() {
        // 5M 2I 3M → qlen = 10
        let cigar = [(CIGAR_MATCH, 5u32), (CIGAR_INS, 2), (CIGAR_MATCH, 3)];
        assert_eq!(cigar_qlen(cigar.iter().copied()), 10);
    }

    #[test]
    fn cigar_intron_basic() {
        // 100M 500N 100M → intron = 500
        let cigar = [
            (CIGAR_MATCH, 100u32),
            (CIGAR_REF_SKIP, 500),
            (CIGAR_MATCH, 100),
        ];
        assert_eq!(cigar_intron_len(cigar.iter().copied()), 500);
    }

    #[test]
    fn overlap_bases_non_overlapping() {
        // read1: 1000..1100, read1_end=1100, read2_start=1200 → no overlap
        let cigar = [(CIGAR_MATCH, 100u32)];
        assert_eq!(
            count_overlap_exon_bases(1000, 1100, 1200, cigar.iter().copied()),
            0
        );
    }

    #[test]
    fn overlap_bases_full() {
        // read1: 1000..1200, read1_end=1200, read2_start=1100 → overlap [1100,1200) = 100
        let cigar = [(CIGAR_MATCH, 200u32)];
        assert_eq!(
            count_overlap_exon_bases(1000, 1200, 1100, cigar.iter().copied()),
            100
        );
    }
}
