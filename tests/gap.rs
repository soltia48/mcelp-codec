//! The minimum-gap pass of the LSF quantiser.
//!
//! Columns: five values in, five values out.

use mcelp::lsf_weight::enforce_spacing;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let half = v.len() / 2;
        let mut values = v[..half].to_vec();
        enforce_spacing(&mut values);
        assert_eq!(values, v[half..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn gap_sequence() {
    let path = format!("{}/tests/data/gap.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} passes verified", replay(&path));
}
