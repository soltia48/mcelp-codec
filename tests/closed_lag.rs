//! The best integer lag in the window.
//!
//! Columns: the low and high bounds, 32 correlations, then the lag chosen.

use mcelp::pitch_search::best_closed_lag;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        assert_eq!(best_closed_lag(&v[2..34], v[0], v[1]), v[34], "line {i}");
    }
    text.lines().count()
}

#[test]
fn closed_sequence() {
    let path = format!("{}/tests/data/closed_lag.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} lags verified", replay(&path));
}
