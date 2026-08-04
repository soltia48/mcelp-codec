//! The fourteen fields packed into nine words.
//!
//! Columns: fifteen words from the field array (the fourteen fields and the
//! suppression flag), then the nine packed words.

use mcelp::bitstream::{FIELDS, Params, pack};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let params = Params {
            field: v[..FIELDS].try_into().unwrap(),
            suppress: v[FIELDS] != 0,
        };
        assert_eq!(pack(&params).as_slice(), &v[FIELDS + 1..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn pack_excerpt() {
    assert_eq!(replay("tests/data/pack.txt"), 120);
}
