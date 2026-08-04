//! The midpoint of two line spectral pairs.
//!
//! Columns: the first vector, the second, then their midpoint.

use mcelp::lpc::ORDER;
use mcelp::weight_lpc::midpoint;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut a = [0i16; ORDER];
        let mut b = [0i16; ORDER];
        a.copy_from_slice(&v[..ORDER]);
        b.copy_from_slice(&v[ORDER..2 * ORDER]);
        assert_eq!(
            midpoint(&a, &b).to_vec(),
            v[2 * ORDER..].to_vec(),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn midpoint_sequence() {
    let path = format!("{}/tests/data/midpoint.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
