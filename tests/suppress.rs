//! The noise-suppression gain applied to the complex spectrum.
//!
//! Each line is one frame: the analysis headroom, the 129 real and 129
//! imaginary bins going in, the 129 magnitude bins, the 129 noise-floor bins,
//! then the real and imaginary bins the reference encoder produced.

use mcelp::analysis::{BINS, Spectrum};
use mcelp::shaping::suppress;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let at = |n: usize| -> [i16; BINS] {
            let mut a = [0i16; BINS];
            a.copy_from_slice(&v[1 + n * BINS..1 + (n + 1) * BINS]);
            a
        };
        let mut s = Spectrum {
            re: at(0),
            im: at(1),
        };
        suppress(&mut s, &at(2), &at(3), v[0]);
        assert_eq!(s.re.to_vec(), at(4).to_vec(), "re, line {i}");
        assert_eq!(s.im.to_vec(), at(5).to_vec(), "im, line {i}");
    }
    text.lines().count()
}

#[test]
fn suppression_sequence() {
    let path = format!("{}/tests/data/suppress.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
