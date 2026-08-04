//! Scaling the signal for the open-loop search.
//!
//! Columns: 343 samples in, 343 out.

use mcelp::pitch_search::{SEARCH_SPAN, normalise};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut signal = [0i16; SEARCH_SPAN];
        signal.copy_from_slice(&v[..SEARCH_SPAN]);
        assert_eq!(
            normalise(&signal).to_vec(),
            v[SEARCH_SPAN..].to_vec(),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn normalise_sequence() {
    let path = format!("{}/tests/data/normalise.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} blocks verified", replay(&path));
}
