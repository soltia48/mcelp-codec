//! Energy and two correlations, normalised.
//!
//! Columns: 80 target samples, 80 samples before, 80 after, then three
//! mantissas and three exponents.

use mcelp::convolve::SUBFRAME;
use mcelp::pitch_search::measures;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let got = measures(
            &v[..SUBFRAME],
            &v[SUBFRAME..2 * SUBFRAME],
            &v[2 * SUBFRAME..3 * SUBFRAME],
        );
        let base = 3 * SUBFRAME;
        for k in 0..3 {
            assert_eq!(
                (got[k].mantissa, got[k].exponent),
                (v[base + k], v[base + 3 + k]),
                "measure {k}, line {i}"
            );
        }
    }
    text.lines().count()
}

#[test]
fn measure_sequence() {
    let path = format!("{}/tests/data/measures.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} subframes verified", replay(&path));
}
