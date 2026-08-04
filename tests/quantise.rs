//! The LSF quantiser.
//!
//! Columns: ten line spectral frequencies, ten weights, the thirty words of
//! predictor memory, then the two transport fields the reference produced.

use mcelp::lpc::ORDER;
use mcelp::lsf_weight::quantise;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut lsf = [0i16; ORDER];
        let mut weight = [0i16; ORDER];
        lsf.copy_from_slice(&v[..ORDER]);
        weight.copy_from_slice(&v[ORDER..2 * ORDER]);
        let mut history = [[0i16; ORDER]; 3];
        for (j, row) in history.iter_mut().enumerate() {
            row.copy_from_slice(&v[2 * ORDER + j * ORDER..2 * ORDER + (j + 1) * ORDER]);
        }
        let got = quantise(&lsf, &weight, &history);
        assert_eq!(got.fields(), (v[5 * ORDER], v[5 * ORDER + 1]), "line {i}");
    }
    text.lines().count()
}

#[test]
fn quantise_sequence() {
    let path = format!("{}/tests/data/quantise.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
