//! The whole noise-suppression path in one go.
//!
//! Columns: headroom, the voicing measure, 129 magnitude bins, the 129 real and
//! 129 imaginary transform bins, then the 160 shaping samples the reference
//! encoder emitted.  Every stage's state is carried across lines, so this
//! checks the composition as much as the arithmetic.

use mcelp::analysis::{BINS, Spectrum};
use mcelp::frontend::Suppressor;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut suppressor = Suppressor::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let at = |n: usize| -> [i16; BINS] {
            let mut a = [0i16; BINS];
            a.copy_from_slice(&v[2 + n * BINS..2 + (n + 1) * BINS]);
            a
        };
        let magnitude = at(0);
        let mut spectrum = Spectrum {
            re: at(1),
            im: at(2),
        };
        let got = suppressor.run(&mut spectrum, &magnitude, v[0], v[1]);
        assert_eq!(got.to_vec(), v[2 + 3 * BINS..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn pipeline_sequence() {
    let path = format!("{}/tests/data/pipeline.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} half-frames verified", replay(&path));
}
