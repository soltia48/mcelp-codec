//! Square root by binary search.
//!
//! Columns: the value and the root the reference produced.

use mcelp::fixed::sqrt;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(sqrt(v[0]), v[1], "line {i}");
    }
    text.lines().count()
}

#[test]
fn sqrt_sequence() {
    let path = format!("{}/tests/data/sqrt.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} roots verified", replay(&path));
}
