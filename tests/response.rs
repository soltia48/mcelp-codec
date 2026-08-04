//! The adaptive codebook's filtered response.
//!
//! Columns: the lag, the fractional part, the gain, the eighty impulse response
//! taps and the eighty samples that come out.

use mcelp::shape_pair::response;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let got = response(&v[3..3 + SPAN], v[0], v[1], v[2]);
        assert_eq!(got.as_slice(), &v[3 + SPAN..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn response_excerpt() {
    assert_eq!(replay("tests/data/response.txt"), 120);
}
