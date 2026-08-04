//! The adaptive contribution taken out of the target.
//!
//! Columns: the normalised target, the contribution, then the gain, the
//! reciprocal mantissa and its leftover shift, the scaled contribution and the
//! residual.

use mcelp::pulses::remove_adaptive;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut contribution: [i16; SPAN] = v[SPAN..2 * SPAN].try_into().unwrap();
        let mut residual = [0i16; SPAN];
        let (inverse, left, gain) = remove_adaptive(&v[..SPAN], &mut contribution, &mut residual);

        let at = 2 * SPAN;
        assert_eq!(
            (gain, inverse, left),
            (v[at], v[at + 1], v[at + 2]),
            "line {i} words"
        );
        assert_eq!(
            contribution.as_slice(),
            &v[at + 3..at + 3 + SPAN],
            "line {i} contribution"
        );
        assert_eq!(
            residual.as_slice(),
            &v[at + 3 + SPAN..],
            "line {i} residual"
        );
    }
    text.lines().count()
}

#[test]
fn remove_adaptive_excerpt() {
    assert_eq!(replay("tests/data/remove_adaptive.txt"), 120);
}
