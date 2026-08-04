//! Best lag in a range and its normalised correlation.
//!
//! Columns: the longest lag, the count, 343 samples, then the lag and value the
//! reference produced.

use mcelp::pitch_search::{SEARCH_SPAN, normalised_peak};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let signal: Vec<i16> = v[2..2 + SEARCH_SPAN].iter().map(|&x| x as i16).collect();
        let (lag, value) = normalised_peak(&signal, v[0] as i16, v[1] as i16);
        assert_eq!(
            (lag as i64, value),
            (v[2 + SEARCH_SPAN], v[3 + SEARCH_SPAN]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn peak_sequence() {
    let path = format!("{}/tests/data/peak.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} peaks verified", replay(&path));
}
