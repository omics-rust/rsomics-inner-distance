# rsomics-inner-distance

Compute the mRNA-aware **inner distance** (insert size) distribution of
paired-end RNA-seq fragments from a BAM + BED12 gene model — a Rust port of
RSeQC's `inner_distance.py`. For pairs whose mates fall in the same transcript,
the spliced (mRNA) distance is used; otherwise the genomic distance.

Byte-exact with RSeQC 5.0.4 (per-pair table + frequency histogram) and faster
single-threaded than the Python upstream.

## Install

```sh
cargo install rsomics-inner-distance
```

## Usage

```sh
rsomics-inner-distance -i sample.bam -o prefix -r genes.bed12
```

Writes `<prefix>.inner_distance.txt` (per-pair distances + category) and
`<prefix>.inner_distance_freq.txt` (histogram).

| flag | meaning | default |
|---|---|---|
| `-i, --input` | coordinate-sorted, indexed BAM | required |
| `-o, --out-prefix` | output file prefix | required |
| `-r, --refgene` | reference gene model (BED12) | required |
| `-k, --sample-size` | read pairs sampled | 1000000 |
| `-l/-u, --lower/upper-bound` | histogram bounds (bp) | -250 / 250 |
| `-s, --step` | histogram bin width (bp) | 5 |
| `--mapq` | minimum mapping quality | 30 |

## Origin

Independent Rust reimplementation of `inner_distance.py` (RSeQC) based on the
published method, the BED12 + SAM/BAM specs, and black-box testing against
`inner_distance.py` 5.0.4. To stay byte-exact, it reproduces RSeQC's actual
behavior (including its first-transcript-per-chromosome handling). No GPL RSeQC
source was used.

License: MIT OR Apache-2.0.
Upstream credit: [RSeQC](https://rseqc.sourceforge.net/) (GNU GPL).
