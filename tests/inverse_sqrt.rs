//! The reciprocal of a square root.
//!
//! Columns: the value and the result the reference produced.

use mcelp::fixed::inverse_sqrt;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(inverse_sqrt(v[0]), v[1], "line {i}");
    }
    text.lines().count()
}

#[test]
fn inverse_sqrt_sequence() {
    let path = format!("{}/tests/data/inverse_sqrt.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} values verified", replay(&path));
}
