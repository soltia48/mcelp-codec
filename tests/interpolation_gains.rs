//! The two extension gains for the coming half-frame.
//!
//! Columns: the first reflection coefficient, the previous gain, the lag, the
//! previous lag, then the two gains.

use mcelp::excitation::interpolation_gains;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(
            interpolation_gains(v[0], v[1], v[2], v[3]),
            (v[4], v[5]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn interpolation_gains_excerpt() {
    assert_eq!(replay("tests/data/interpolation_gains.txt"), 120);
}
