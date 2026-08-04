//! Shape filtering with eighty taps: the excitation convolved with the response.

use mcelp::convolve::filter;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let got = filter(&v[..SPAN], &v[SPAN..2 * SPAN]);
        assert_eq!(got.as_slice(), &v[2 * SPAN..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn pulse_response_excerpt() {
    assert_eq!(replay("tests/data/pulse_response.txt"), 60);
}
