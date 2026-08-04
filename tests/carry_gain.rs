//! The adaptive gain clamped for the next half-frame.

use mcelp::excitation::carry_gain;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(carry_gain(v[0]), v[1], "line {i}");
    }
    text.lines().count()
}

#[test]
fn carry_gain_excerpt() {
    assert_eq!(replay("tests/data/carry_gain.txt"), 120);
}
