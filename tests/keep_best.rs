//! The better of this mode and the best so far.
//!
//! Columns: the mode, this mode's pair, the running best (pair and mode), the
//! chosen positions and flags, the kept positions and flags going in, then the
//! best and the kept arrays coming out.

use mcelp::pulses::keep_best;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let mut best = (v[3], v[4], v[5]);
        let mut positions: [i16; 3] = v[12..15].try_into().unwrap();
        let mut flags: [i16; 3] = v[15..18].try_into().unwrap();

        keep_best(
            mode,
            (v[1], v[2]),
            &mut best,
            &v[6..9],
            &v[9..12],
            &mut positions,
            &mut flags,
        );
        assert_eq!([best.0, best.1, best.2], v[18..21], "line {i} best");
        assert_eq!(positions, v[21..24], "line {i} positions");
        assert_eq!(flags, v[24..27], "line {i} flags");
    }
    text.lines().count()
}

#[test]
fn keep_best_excerpt() {
    assert_eq!(replay("tests/data/keep_best.txt"), 120);
}
