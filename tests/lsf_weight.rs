//! The weights the LSF quantiser measures distance with.
//!
//! Columns: the ten line spectral frequencies, then the ten weights.

use mcelp::lpc::ORDER;
use mcelp::lsf_weight::weights;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut lsf = [0i16; ORDER];
        lsf.copy_from_slice(&v[..ORDER]);
        assert_eq!(weights(&lsf).to_vec(), v[ORDER..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn weight_sequence() {
    let path = format!("{}/tests/data/lsf_weight.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
