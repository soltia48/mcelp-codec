//! The LPC analysis composed: window to line spectrum.
//!
//! Columns: four hundred windowed samples — eighty from the previous frame's
//! tail and both of this frame's half-frames — then the two half-frames' ten
//! line spectral frequencies.  The window ends on the second half-frame, so it
//! is that one's line spectrum the chain reproduces.

use mcelp::lpc::{
    ANALYSIS, ORDER, autocorrelate, lag_window, levinson, line_spectrum, split_polynomials,
};

fn analyse(window: &[i16; ANALYSIS]) -> [i16; ORDER] {
    let mut r = autocorrelate(window);
    lag_window(&mut r);
    let (a, _) = levinson(&r);
    let (sum, difference) = split_polynomials(&a);
    line_spectrum(&sum, &difference)
}

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let window: [i16; ANALYSIS] = v[..ANALYSIS].try_into().unwrap();
        assert_eq!(analyse(&window), v[ANALYSIS + ORDER..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn lpc_chain_excerpt() {
    assert_eq!(replay("tests/data/lpc_chain.txt"), 40);
}
