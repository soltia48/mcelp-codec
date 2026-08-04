//! One filtered response per track.
//!
//! Columns: the mode, the lag, the fractional part, the gain, the eighty
//! impulse response taps, then the responses, eighty samples per track.

use mcelp::pulses::track_responses;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let got = track_responses(v[0] as usize, &v[4..4 + SPAN], v[1], v[2], v[3]);
        assert_eq!(got.as_slice(), &v[4 + SPAN..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn track_responses_excerpt() {
    assert_eq!(replay("tests/data/track_responses.txt"), 60);
}
