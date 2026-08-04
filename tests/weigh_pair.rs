//! The shape pair's result weighted and aligned.
//!
//! Columns: the weight's index in the weight table, the pair, the two shifts, then the three words
//! coming out.

use mcelp::pulses::weigh_pair;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i32> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let at = v[0] as usize;
        let mut pair = (v[1] as i16, v[2] as i16);
        let kept = weigh_pair(
            &mut pair,
            mcelp::tables::MODE_WEIGHTS[at],
            mcelp::tables::MODE_WEIGHTS[at - 5],
            v[3] as i16,
            v[4] as i16,
        );
        assert_eq!(
            (kept as i32, pair.0 as i32, pair.1 as i32),
            (v[5], v[6], v[7]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn weigh_pair_excerpt() {
    assert_eq!(replay("tests/data/weigh_pair.txt"), 120);
}
