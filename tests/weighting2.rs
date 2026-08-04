//! The shaping path's second and third stages.
//!
//! Each line is one half-frame: the 160 samples the high-pass produced, the 160
//! the smoothing filter produced, then the 160 the differencing pass left.

use mcelp::weighting::{Smooth, Tilt};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut smooth = Smooth::default();
    let mut tilt = Tilt::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut block = v[..160].to_vec();
        smooth.run(&mut block);
        assert_eq!(block, v[160..320].to_vec(), "smooth, line {i}");
        tilt.run(&mut block);
        assert_eq!(block, v[320..].to_vec(), "tilt, line {i}");
    }
    text.lines().count()
}

#[test]
fn smooth_sequence() {
    let path = format!("{}/tests/data/weighting2.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} half-frames verified", replay(&path));
}
