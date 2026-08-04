//! Splitting the prediction filter in two.
//!
//! Columns: the eleven prediction coefficients, then the six terms of each
//! half-polynomial the reference encoder built.

use mcelp::lpc::{ORDER, split_polynomials};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut a = [0i16; ORDER + 1];
        a.copy_from_slice(&v[..ORDER + 1]);
        let (sum, difference) = split_polynomials(&a);
        assert_eq!(
            sum.to_vec(),
            v[ORDER + 1..ORDER + 7].to_vec(),
            "sum, line {i}"
        );
        assert_eq!(
            difference.to_vec(),
            v[ORDER + 7..].to_vec(),
            "difference, line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn polynomial_sequence() {
    let path = format!("{}/tests/data/polynomials.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
