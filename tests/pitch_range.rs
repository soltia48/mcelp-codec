//! The lag window each subframe searches in.
//!
//! Columns: which subframe (0 or 1), the current low and high bounds, the
//! previous lag, the base offset, then the bounds and origin the reference set.

use mcelp::pitch_search::{first, second};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let (which, low, _high, lag, base) = (v[0], v[1], v[2], v[3], v[4]);
        if which == 0 {
            let w = first(lag, base);
            assert_eq!((w.low, w.high, w.origin), (v[5], v[6], v[7]), "line {i}");
        } else {
            assert_eq!(second(low, lag, base), v[7], "line {i}");
        }
    }
    text.lines().count()
}

#[test]
fn range_sequence() {
    let path = format!("{}/tests/data/pitch_range.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} windows verified", replay(&path));
}
