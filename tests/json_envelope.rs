//! `--json` must emit exactly ONE envelope document on stdout carrying a
//! populated `result` object — not the manual summary followed by a second
//! `result: null` envelope (which made `json.load` fail with "Extra data").
//!
//! Uses a COMMITTED expectation (golden fixture has 4 read pairs); runs no live
//! oracle so it is portable across CI runners.

use std::path::Path;
use std::process::Command;

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");
const BIN: &str = env!("CARGO_BIN_EXE_rsomics-inner-distance");

#[test]
fn json_is_single_doc_with_populated_result() {
    let bam = Path::new(GOLDEN).join("pairs.bam");
    let bed = Path::new(GOLDEN).join("genes.bed12");
    if !bam.exists() || !bed.exists() {
        eprintln!("SKIP: golden fixture not found");
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("json_test").to_string_lossy().into_owned();

    let out = Command::new(BIN)
        .args([
            "-i",
            bam.to_str().unwrap(),
            "-r",
            bed.to_str().unwrap(),
            "-o",
            &prefix,
            "--mapq",
            "30",
            "--json",
        ])
        .output()
        .expect("failed to run rsomics-inner-distance");

    assert!(
        out.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // A single valid JSON document — `from_str` rejects trailing "Extra data".
    let stdout = String::from_utf8(out.stdout).unwrap();
    let doc: serde_json::Value = serde_json::from_str(stdout.trim())
        .unwrap_or_else(|e| panic!("stdout not a single JSON document ({e}): {stdout:?}"));

    assert_eq!(doc["status"], "ok");
    let result = &doc["result"];
    assert!(
        !result.is_null(),
        "result must be populated, got null: {doc}"
    );
    assert_eq!(result["pair_num"], 4, "golden fixture has 4 read pairs");
    assert_eq!(result["output_prefix"], prefix);
}
