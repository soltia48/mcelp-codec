//! The shaped spectrum on its way back to the time domain.
//!
//! Each line is one frame: 129 real bins, 129 imaginary bins, then the 256
//! samples the reference encoder produced.  Every frame stands on its own, so
//! nothing is carried between lines.

use mcelp::analysis::BINS;
use mcelp::shaping::{SPAN, inverse_transform};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut re = [0i16; BINS];
        let mut im = [0i16; BINS];
        re.copy_from_slice(&v[..BINS]);
        im.copy_from_slice(&v[BINS..2 * BINS]);
        let got = inverse_transform(&re, &im);
        assert_eq!(
            got.to_vec(),
            v[2 * BINS..2 * BINS + SPAN].to_vec(),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn inverse_sequence() {
    let path = format!("{}/tests/data/inverse.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
