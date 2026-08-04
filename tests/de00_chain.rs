//! The adaptive-codebook measures composed: cross-correlation then the gain.
//!
//! Columns: the target, the filtered adaptive vector, then the gain and the
//! energy's mantissa and exponent.  The second block lines its energy up
//! against the scaling the correlation was taken at, not the pair's exponent.

use mcelp::convolve::SUBFRAME;
use mcelp::pitch_search::{adaptive_gain, cross};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let (rounded, _, against) = cross(&v[..SUBFRAME], &v[SUBFRAME..2 * SUBFRAME]);
        let (gain, energy) = adaptive_gain(&v[SUBFRAME..2 * SUBFRAME], rounded, against);
        let at = 2 * SUBFRAME;
        assert_eq!(
            (gain, energy.mantissa, energy.exponent),
            (v[at], v[at + 1], v[at + 2]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn de00_chain_excerpt() {
    assert_eq!(replay("tests/data/de00_chain.txt"), 120);
}
