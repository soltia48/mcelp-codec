//! The per-subframe coding mode.
//!
//! Columns: eighty target samples, eighty filtered ones, the adaptive gain, the
//! incoming state (previous energy and hold counter), then the chosen mode and
//! the state on the way out.

use mcelp::mode::Mode;

const SUBFRAME: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut state = Mode::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let w: Vec<i16> = v[..2 * SUBFRAME + 1].iter().map(|&x| x as i16).collect();
        state.restore(v[2 * SUBFRAME + 1], v[2 * SUBFRAME + 2] as i16);
        let chosen = state.choose(&w[..SUBFRAME], &w[SUBFRAME..2 * SUBFRAME], w[2 * SUBFRAME]);
        let base = 2 * SUBFRAME + 3;
        assert_eq!(
            (chosen as i64, state.energy(), state.hold() as i64),
            (v[base], v[base + 1], v[base + 2]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn mode_excerpt() {
    assert_eq!(replay("tests/data/mode.txt"), 120);
}
