//! The lag with the largest correlation.
//!
//! Columns: the longest lag tried, how many below it, 343 normalised samples,
//! then the lag the reference picked.

use mcelp::pitch_search::{SEARCH_SPAN, best_lag};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let signal = &v[2..2 + SEARCH_SPAN];
        assert_eq!(best_lag(signal, v[0], v[1]), v[2 + SEARCH_SPAN], "line {i}");
    }
    text.lines().count()
}

#[test]
fn lag_sequence() {
    let path = format!("{}/tests/data/best_lag.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} searches verified", replay(&path));
}
