//! The second-stage search over the lower half.
//!
//! Columns: five residual values, five weights, then the entry chosen.

use mcelp::lsf_weight::search_stage2;

fn replay(path: &str, offset: usize) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut residual = [0i16; 5];
        let mut weight = [0i16; 5];
        residual.copy_from_slice(&v[..5]);
        weight.copy_from_slice(&v[5..10]);
        assert_eq!(search_stage2(&residual, &weight, offset), v[10], "line {i}");
    }
    text.lines().count()
}

#[test]
fn stage2_sequence() {
    let root = env!("CARGO_MANIFEST_DIR");
    let lower = replay(&format!("{root}/tests/data/stage2.txt"), 0);
    let upper = replay(&format!("{root}/tests/data/stage2b.txt"), 5);
    println!("{lower} + {upper} searches verified");
}
