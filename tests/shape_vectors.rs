//! The shortlisted shapes laid out as vectors.
//!
//! Columns: the eighty response taps, the eight shape numbers, then the eight
//! eighty-sample vectors.

use mcelp::shape_pair::{quarter, shape_vectors};

const SPAN: usize = 80;
const SHORTLIST: usize = 8;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut response: [i16; SPAN] = v[..SPAN].try_into().unwrap();
        quarter(&mut response);
        let got = shape_vectors(
            &response,
            &v[SPAN..SPAN + SHORTLIST],
            &mcelp::tables::FCB_ALT_A,
        );
        assert_eq!(got.as_slice(), &v[SPAN + SHORTLIST..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn shape_vectors_excerpt() {
    assert_eq!(replay("tests/data/shape_vectors.txt"), 40);
}
