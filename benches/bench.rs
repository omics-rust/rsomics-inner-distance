use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use tempfile;

fn bench_inner_distance(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-inner-distance");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/pairs.bam");
    let bed = manifest.join("tests/golden/genes.bed12");
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    c.bench_function("rsomics-inner-distance golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([
                    "-i",
                    bam.to_str().unwrap(),
                    "-r",
                    bed.to_str().unwrap(),
                    "-o",
                    prefix.to_str().unwrap(),
                ])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_inner_distance);
criterion_main!(benches);
