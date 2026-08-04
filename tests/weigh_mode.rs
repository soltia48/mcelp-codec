//! A mode's score weighted and the pair aligned.
//!
//! Columns: the weight's index in the weight table, the pair, the two shifts,
//! then the next index and the pair coming out.

use mcelp::pulses::weigh_mode;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i32> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let at = v[0] as usize;
        let mut pair = (v[1] as i16, v[2] as i16);
        weigh_mode(
            &mut pair,
            mcelp::tables::MODE_WEIGHTS[at],
            v[3] as i16,
            v[4] as i16,
        );
        assert_eq!(at + 1, v[5] as usize, "line {i} next index");
        assert_eq!(
            (pair.0 as i32, pair.1 as i32),
            (v[6], v[7]),
            "line {i} pair"
        );
    }
    text.lines().count()
}

#[test]
fn weigh_mode_excerpt() {
    assert_eq!(replay("tests/data/weigh_mode.txt"), 120);
}
