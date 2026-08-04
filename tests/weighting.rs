//! The shaping signal's high-pass.
//!
//! Each line is one half-frame: 160 samples in, then the 160 the reference
//! encoder produced.  Runs are chained so the filter memory is checked too.

use mcelp::weighting::Highpass;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut filter = Highpass::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut block = v[..160].to_vec();
        filter.run(&mut block);
        assert_eq!(block, v[160..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn highpass_sequence() {
    let path = format!("{}/tests/data/weighting.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} half-frames verified", replay(&path));
}
