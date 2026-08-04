//! Innovation building as the encoder calls it: the chosen pulses laid out.
//!
//! Columns: the mode, three chosen slots, three sign flags, the track
//! responses, then the eighty-sample innovation.

use mcelp::pulses::{place_pulses, tracks};

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let n = SPAN * tracks(mode);
        let got = place_pulses(mode, &v[1..4], &v[4..7], &v[7..7 + n]);
        assert_eq!(got.as_slice(), &v[7 + n..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn place_pulses_excerpt() {
    assert_eq!(replay("tests/data/place_pulses.txt"), 60);
}
