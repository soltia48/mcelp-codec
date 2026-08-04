//! The interpolated correlation at a fractional lag.
//!
//! Columns: the fraction, the centre index within the window, nine
//! correlations around it, then the reference's result.

use mcelp::pitch_search::interpolate;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let window: Vec<i16> = v[2..11].iter().map(|&x| x as i16).collect();
        assert_eq!(
            interpolate(&window, v[1] as usize, v[0] as i16),
            v[11],
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn interpolate_sequence() {
    let path = format!("{}/tests/data/interpolate.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} values verified", replay(&path));
}
