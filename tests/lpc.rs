//! Autocorrelation and lag windowing.
//!
//! Each line is one frame: 400 signal samples, the eleven correlations the
//! reference encoder computed, then the same after the lag window.

use mcelp::lpc::{ANALYSIS, ORDER, autocorrelate, lag_window};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut signal = [0i16; ANALYSIS];
        for (s, &x) in signal.iter_mut().zip(&v[..ANALYSIS]) {
            *s = x as i16;
        }
        let mut r = autocorrelate(&signal);
        assert_eq!(
            r.to_vec(),
            v[ANALYSIS..ANALYSIS + ORDER + 1].to_vec(),
            "acf, line {i}"
        );
        lag_window(&mut r);
        assert_eq!(
            r.to_vec(),
            v[ANALYSIS + ORDER + 1..].to_vec(),
            "lag, line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn lpc_sequence() {
    let path = format!("{}/tests/data/lpc.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
