//! Bandwidth expansion of the prediction filter.
//!
//! Columns: eleven coefficients, the factor, then the eleven expanded ones.

use mcelp::lpc::ORDER;
use mcelp::weight_lpc::expand;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut a = [0i16; ORDER + 1];
        a.copy_from_slice(&v[..ORDER + 1]);
        assert_eq!(
            expand(&a, v[ORDER + 1]).to_vec(),
            v[ORDER + 2..].to_vec(),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn expand_sequence() {
    let path = format!("{}/tests/data/expand.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} expansions verified", replay(&path));
}
