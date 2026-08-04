//! The adaptive-codebook gain.
//!
//! Columns: 80 samples, the correlation and its exponent, then the gain and the
//! energy's mantissa and exponent.

use mcelp::convolve::SUBFRAME;
use mcelp::pitch_search::adaptive_gain;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let (gain, pair) = adaptive_gain(&v[..SUBFRAME], v[SUBFRAME], v[SUBFRAME + 1]);
        let base = SUBFRAME + 2;
        assert_eq!(
            (gain, pair.mantissa, pair.exponent),
            (v[base], v[base + 1], v[base + 2]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn gain_excerpt() {
    assert_eq!(replay("tests/data/adaptive_gain.txt"), 120);
}
