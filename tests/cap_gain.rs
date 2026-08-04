//! The adaptive gain capped on a periodic half-frame.
//!
//! Columns: the periodicity flag, the gain going in, the gain coming out and
//! the flag as it was stored.

use mcelp::excitation::cap_gain;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut gain = v[1] as i16;
        cap_gain(v[0] as i16, &mut gain);
        assert_eq!((gain as i64, v[0] as i16 as i64), (v[2], v[3]), "line {i}");
    }
    text.lines().count()
}

#[test]
fn cap_gain_excerpt() {
    assert_eq!(replay("tests/data/cap_gain.txt"), 120);
}
