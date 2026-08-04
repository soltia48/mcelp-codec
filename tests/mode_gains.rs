//! The gain each mode's extension uses.

use mcelp::pulses::mode_gains;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(mode_gains(v[0], v[1], v[2]).as_slice(), &v[3..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn mode_gains_excerpt() {
    assert_eq!(replay("tests/data/mode_gains.txt"), 120);
}
