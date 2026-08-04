//! Evaluating a half-polynomial on the unit circle.
//!
//! Columns: five coefficients, the point, and the value the reference produced.

use mcelp::lpc::chebyshev;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut poly = [0i16; 5];
        for (p, &c) in poly.iter_mut().zip(&v[..5]) {
            *p = c as i16;
        }
        assert_eq!(chebyshev(&poly, v[5] as i16), v[6], "line {i}");
    }
    text.lines().count()
}

#[test]
fn chebyshev_sequence() {
    let path = format!("{}/tests/data/chebyshev.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} evaluations verified", replay(&path));
}
