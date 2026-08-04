//! The first-stage codebook search of the LSF quantiser.
//!
//! Columns: the ten values being matched, then the entry the reference chose.

use mcelp::lpc::ORDER;
use mcelp::lsf_weight::search_stage1;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut target = [0i16; ORDER];
        target.copy_from_slice(&v[..ORDER]);
        assert_eq!(search_stage1(&target), v[ORDER], "line {i}");
    }
    text.lines().count()
}

#[test]
fn stage1_sequence() {
    let path = format!("{}/tests/data/stage1.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} searches verified", replay(&path));
}
