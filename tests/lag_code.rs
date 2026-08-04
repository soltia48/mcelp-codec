//! The lag as it goes into the bitstream.
//!
//! Columns: the subframe offset (0 for the first, 80 for the second), the lag,
//! its fractional part, the previous subframe's lag, then the code.

use mcelp::pitch::{Lag, encode_absolute, encode_relative};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let lag = Lag {
            integer: v[1],
            frac: v[2],
        };
        let got = if v[0] == 0 {
            encode_absolute(&lag)
        } else {
            encode_relative(&lag, v[3])
        };
        assert_eq!(got, v[4], "line {i}");
    }
    text.lines().count()
}

#[test]
fn lag_code_excerpt() {
    assert_eq!(replay("tests/data/lag_code.txt"), 120);
}
