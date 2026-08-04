//! Periodic extension of the adaptive codebook vector.
//!
//! Columns: the lag, the fractional part, the gain, the eighty samples going in
//! and the eighty coming out.

use mcelp::SUBFRAME;
use mcelp::pitch::extend;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut buffer: [i16; SUBFRAME] = v[3..3 + SUBFRAME].try_into().unwrap();
        extend(&mut buffer, v[0], v[1], v[2]);
        assert_eq!(buffer.as_slice(), &v[3 + SUBFRAME..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn extend_excerpt() {
    assert_eq!(replay("tests/data/extend.txt"), 120);
}
