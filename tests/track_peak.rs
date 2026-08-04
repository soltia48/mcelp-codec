//! The running record of recent correlation peaks.
//!
//! Columns: the lag, the gain, 64 stored peaks, then the first four after the
//! update.

use mcelp::pitch_search::track_peak;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut peaks: Vec<i64> = v[2..66].to_vec();
        track_peak(&mut peaks, v[0] as i16, v[1] as i16);
        assert_eq!(peaks[..4].to_vec(), v[66..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn track_sequence() {
    let path = format!("{}/tests/data/track_peak.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} updates verified", replay(&path));
}
