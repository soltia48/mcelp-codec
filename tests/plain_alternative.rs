//! The residual on its own against the shape pair.
//!
//! Columns: the kept weight, the pair going in, the eighty residual samples,
//! then the pair coming out.

use mcelp::pulses::plain_alternative;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut pair = (v[1], v[2]);
        plain_alternative(&mut pair, v[0], &v[3..3 + SPAN]);
        assert_eq!([pair.0, pair.1], v[3 + SPAN..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn plain_alternative_excerpt() {
    assert_eq!(replay("tests/data/plain_alternative.txt"), 120);
}
