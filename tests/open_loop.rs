//! The open-loop pitch lag.
//!
//! Columns: 343 samples, then the lag the reference settled on.

use mcelp::pitch_search::{SEARCH_SPAN, open_loop_lag};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut signal = [0i16; SEARCH_SPAN];
        signal.copy_from_slice(&v[..SEARCH_SPAN]);
        assert_eq!(open_loop_lag(&signal), v[SEARCH_SPAN], "line {i}");
    }
    text.lines().count()
}

#[test]
fn open_loop_sequence() {
    let path = format!("{}/tests/data/open_loop.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
