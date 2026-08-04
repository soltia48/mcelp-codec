//! The filter memories carried to the next subframe.
//!
//! Columns: the two gains, ten input samples, ten synthesised, ten weighted,
//! ten adaptive, ten innovation, then the two ten-tap memories.

use mcelp::excitation::update_memories;

const N: usize = 10;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let s = |k: usize| &v[2 + k * N..2 + (k + 1) * N];
        let (error, residual) = update_memories(s(0), s(1), s(2), s(3), s(4), (v[0], v[1]));
        assert_eq!(error.as_slice(), s(5), "line {i} error");
        assert_eq!(residual.as_slice(), s(6), "line {i} residual");
    }
    text.lines().count()
}

#[test]
fn update_memories_excerpt() {
    assert_eq!(replay("tests/data/update_memories.txt"), 120);
}
