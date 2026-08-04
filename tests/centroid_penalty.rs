//! The mode penalised for misplacing its energy.
//!
//! Columns: the innovation's centre, the target's, the score and headroom going
//! in, then the two coming out.

use mcelp::pulses::centroid_penalty;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let (mut score, mut headroom) = (v[2] as i16, v[3] as i16);
        centroid_penalty(&mut score, &mut headroom, v[0], v[1] as i16);
        assert_eq!((score as i64, headroom as i64), (v[4], v[5]), "line {i}");
    }
    text.lines().count()
}

#[test]
fn centroid_penalty_excerpt() {
    assert_eq!(replay("tests/data/centroid_penalty.txt"), 120);
}
