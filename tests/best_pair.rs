//! One shape chosen from each codebook.
//!
//! Columns: the two blocks of eight shape vectors, sixteen correlations,
//! sixteen energies, the energy shift, the sixteen shape numbers, then the two
//! chosen numbers.

use mcelp::shape_pair::best_pair;

const SPAN: usize = 80;
const BLOCK: usize = SPAN * 8;
const CANDIDATES: usize = 16;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let vectors: Vec<i16> = v[..2 * BLOCK].iter().map(|&x| x as i16).collect();
        let (first, second) = vectors.split_at(BLOCK);
        let mut at = 2 * BLOCK;
        let correlation = &v[at..at + CANDIDATES];
        at += CANDIDATES;
        let energy = &v[at..at + CANDIDATES];
        at += CANDIDATES;
        let up = (v[at] + 1) as i16;
        at += 1;
        let numbers: Vec<i16> = v[at..at + CANDIDATES].iter().map(|&x| x as i16).collect();
        at += CANDIDATES;

        let got = best_pair(
            first,
            second,
            correlation,
            energy,
            up,
            (&numbers[..8], &numbers[8..]),
        );
        assert_eq!(
            (got.0.0 as i64, got.0.1 as i64),
            (v[at], v[at + 1]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn best_pair_excerpt() {
    assert_eq!(replay("tests/data/best_pair.txt"), 40);
}
