//! The encoder's 32-bit divide.
//!
//! Columns: the divisor, the dividend, and the quotient the reference produced.

use mcelp::lpc::divide;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(divide(v[1], v[0]), v[2], "line {i}");
    }
    text.lines().count()
}

#[test]
fn divide_sequence() {
    let path = format!("{}/tests/data/divide.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} divisions verified", replay(&path));
}
