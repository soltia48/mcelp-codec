//! The shortlisted positions' quantities gathered.
//!
//! Columns: the mode, sixty shortlist positions, the three eighty-word tables
//! indexed by position, then the three gathered arrays.

use mcelp::pulses::{gather, kept, tracks};

const POSITIONS: usize = 80;
const LIST: usize = 60;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let mut at = 1;
        let positions = &v[at..at + LIST];
        at += LIST;
        let value = &v[at..at + POSITIONS];
        at += POSITIONS;
        let adaptive = &v[at..at + POSITIONS];
        at += POSITIONS;
        let sign = &v[at..at + POSITIONS];
        at += POSITIONS;

        let n: usize = (0..tracks(mode)).map(|t| kept(mode, t)).sum();
        let (values, shares, signs) = gather(mode, positions, value, adaptive, sign);
        assert_eq!(values.as_slice(), &v[at..at + n], "line {i} values");
        assert_eq!(shares.as_slice(), &v[at + n..at + 2 * n], "line {i} shares");
        assert_eq!(
            signs.as_slice(),
            &v[at + 2 * n..at + 3 * n],
            "line {i} signs"
        );
    }
    text.lines().count()
}

#[test]
fn gather_excerpt() {
    assert_eq!(replay("tests/data/gather.txt"), 120);
}
