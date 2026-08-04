//! Each shape vector scored.
//!
//! Columns: the eighty target samples, the two blocks of eight shape vectors,
//! then sixteen correlations and their shift, sixteen energies and theirs.

use mcelp::shape_pair::measure_shapes;

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
        let words: Vec<i16> = v[..SPAN + 2 * BLOCK].iter().map(|&x| x as i16).collect();
        let (target, vectors) = words.split_at(SPAN);

        let (correlation, up, energy, over) = measure_shapes(target, vectors);
        let mut at = SPAN + 2 * BLOCK;
        assert_eq!(
            correlation.as_slice(),
            &v[at..at + CANDIDATES],
            "line {i} correlation"
        );
        at += CANDIDATES;
        assert_eq!(up as i64, v[at], "line {i} correlation shift");
        at += 1;
        assert_eq!(
            energy.as_slice(),
            &v[at..at + CANDIDATES],
            "line {i} energy"
        );
        at += CANDIDATES;
        assert_eq!(over as i64, v[at], "line {i} energy shift");
    }
    text.lines().count()
}

#[test]
fn measure_shapes_excerpt() {
    assert_eq!(replay("tests/data/measure_shapes.txt"), 40);
}
