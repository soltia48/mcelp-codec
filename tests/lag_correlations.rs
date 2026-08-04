//! The closed-loop correlation of every lag.
//!
//! Each row is `lags`, the prescale shift, the eighty target samples, the
//! eighty filtered ones, the forty response taps, the `lags - 1` excitation
//! samples the slides bring in, and the `lags` correlations that come out.

use mcelp::convolve::SUBFRAME;
use mcelp::pitch_search::lag_correlations;

const RESPONSE: usize = 40;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i32> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let lags = v[0] as usize;
        let down = v[1];
        let w: Vec<i16> = v[2..].iter().map(|&x| x as i16).collect();
        let (target, rest) = w.split_at(SUBFRAME);
        let (filtered, rest) = rest.split_at(SUBFRAME);
        let (response, rest) = rest.split_at(RESPONSE);
        let (past, expected) = rest.split_at(lags - 1);

        let mut filtered: [i16; SUBFRAME] = filtered.try_into().unwrap();
        let got = lag_correlations(target, &mut filtered, past, response, down, lags);
        assert_eq!(got, expected, "line {i}");
    }
    text.lines().count()
}

#[test]
fn lag_excerpt() {
    assert_eq!(replay("tests/data/lag_correlations.txt"), 60);
}
