//! The eight best pulse shapes from one codebook.
//!
//! Columns: the eighty lag correlations, then the eight shape numbers.

use mcelp::shape_pair::shapes;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let shapes = shapes(&mcelp::tables::FCB_ALT_A, &v[..SPAN]);
        assert_eq!(shapes.as_slice(), &v[SPAN..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn shapes_excerpt() {
    assert_eq!(replay("tests/data/shapes.txt"), 120);
}
