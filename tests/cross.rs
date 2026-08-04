//! The cross-correlation and its scaled pair.
//!
//! Columns: 80 samples of each stretch, then the rounded correlation, the
//! mantissa and the exponent.

use mcelp::convolve::SUBFRAME;
use mcelp::pitch_search::cross;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let (rounded, pair, _against) = cross(&v[..SUBFRAME], &v[SUBFRAME..2 * SUBFRAME]);
        let base = 2 * SUBFRAME;
        assert_eq!(
            (rounded, pair.mantissa, pair.exponent),
            (v[base], v[base + 1], v[base + 2]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn cross_sequence() {
    let path = format!("{}/tests/data/cross.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} subframes verified", replay(&path));
}
