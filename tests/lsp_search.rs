//! Half-polynomials in, line spectral pair out.
//!
//! Columns: the six terms of each half-polynomial, then the ten line spectral
//! frequencies the reference encoder found.

use mcelp::lpc::{ORDER, line_spectrum};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut sum = [0i16; ORDER / 2 + 1];
        let mut difference = [0i16; ORDER / 2 + 1];
        sum.copy_from_slice(&v[..6]);
        difference.copy_from_slice(&v[6..12]);
        assert_eq!(
            line_spectrum(&sum, &difference).to_vec(),
            v[12..].to_vec(),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn lsp_sequence() {
    let path = format!("{}/tests/data/lsp_search.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
