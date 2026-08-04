//! The Levinson-Durbin recursion.
//!
//! Columns: the eleven lag-windowed correlations, the eleven 32-bit prediction
//! coefficients the reference left behind, then its ten reflection coefficients.

use mcelp::lpc::{ORDER, levinson};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut r = [0i64; ORDER + 1];
        r.copy_from_slice(&v[..ORDER + 1]);
        if r[0] == 0 {
            continue;
        }
        let (lpc, reflection) = levinson(&r);
        let want: Vec<i16> = v[2 * (ORDER + 1)..].iter().map(|&x| x as i16).collect();
        assert_eq!(reflection.to_vec(), want, "reflection, line {i}");
        // The reference leaves the coefficients at 32 bits and rounds them on
        // the way out; compare against that rounding.
        let want: Vec<i16> = v[ORDER + 1..2 * (ORDER + 1)]
            .iter()
            .map(|&x| ((x + (1 << 15)) >> 16) as i16)
            .collect();
        assert_eq!(lpc.to_vec(), want, "lpc, line {i}");
    }
    text.lines().count()
}

#[test]
fn levinson_sequence() {
    let path = format!("{}/tests/data/levinson.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
