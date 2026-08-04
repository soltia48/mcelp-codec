//! The weighting filter's shape, decided once per half-frame.
//!
//! Columns: the two reflection coefficients, the incoming state (two smoothed
//! log areas and the flat flag), the two line spectra, then the two numerator
//! and two denominator factors and the state on the way out.

use mcelp::lpc::ORDER;
use mcelp::weight_lpc::Sharpness;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut state = Sharpness::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let reflection: [i16; 2] = v[0..2].try_into().unwrap();
        state.restore([v[2], v[3]], v[4]);
        let first: [i16; ORDER] = v[5..5 + ORDER].try_into().unwrap();
        let second: [i16; ORDER] = v[5 + ORDER..5 + 2 * ORDER].try_into().unwrap();

        let (numerator, denominator) = state.choose(&reflection, &[&first, &second]);
        let base = 5 + 2 * ORDER;
        assert_eq!(
            (numerator[0], numerator[1], denominator[0], denominator[1]),
            (v[base], v[base + 1], v[base + 2], v[base + 3]),
            "line {i} factors"
        );
        assert_eq!(
            state.saved(),
            ([v[base + 4], v[base + 5]], v[base + 6]),
            "line {i} state"
        );
    }
    text.lines().count()
}

#[test]
fn sharpness_excerpt() {
    assert_eq!(replay("tests/data/sharpness.txt"), 120);
}
