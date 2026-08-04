//! The adaptive direction taken out of each shape.
//!
//! Columns: the gain and its shift, the eighty-sample adaptive contribution,
//! the eight shape vectors going in and the eight coming out.

use mcelp::shape_pair::orthogonalise;

const SPAN: usize = 80;
const SHORTLIST: usize = 8;
const BLOCK: usize = SPAN * SHORTLIST;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut vectors = v[2 + SPAN..2 + SPAN + BLOCK].to_vec();
        orthogonalise(&mut vectors, &v[2..2 + SPAN], v[0], v[1]);
        assert_eq!(vectors.as_slice(), &v[2 + SPAN + BLOCK..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn orthogonalise_excerpt() {
    assert_eq!(replay("tests/data/orthogonalise.txt"), 40);
}
