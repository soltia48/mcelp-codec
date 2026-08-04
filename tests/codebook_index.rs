//! The fixed codebook index built from the chosen pulses.
//!
//! Columns: the mode, three chosen indices, three sign bits, then the index.

use mcelp::pulses::codebook_index;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(
            codebook_index(v[0] as usize, &v[1..4], &v[4..7]),
            v[7],
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn index_excerpt() {
    assert_eq!(replay("tests/data/codebook_index.txt"), 120);
}
