//! Shape filtering with a forty-tap shape, which is how the closed-loop search
//! starts its correlations.
//!
//! Columns: eighty excitation samples, the impulse response (of which the
//! first forty taps are used), then the eighty samples that come out.

use mcelp::convolve::filter;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let got = filter(&v[80..120], &v[..80]);
        assert_eq!(got.as_slice(), &v[160..240], "line {i}");
    }
    text.lines().count()
}

#[test]
fn filter40_excerpt() {
    assert_eq!(replay("tests/data/filter40.txt"), 120);
}
