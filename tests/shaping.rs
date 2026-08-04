//! Span stitching and the analysis tail: reference from male.orig.ulaw.
//!
//! Each line is one half-frame: the analysis headroom, the 256-sample span the
//! shaping transform produced, then the 160 samples the reference encoder
//! emitted.  The runs are chained, so the seam carried between spans — and its
//! requantisation when the headroom moves — is checked as well as the mixing.
//! The committed excerpt covers the frames where the headroom walks from 8 down
//! to 0, which is what exercises that requantisation.

use mcelp::shaping::{HOP, Overlap, SPAN};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut overlap = Overlap::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut span = [0i16; SPAN];
        for (s, &x) in span.iter_mut().zip(&v[1..1 + SPAN]) {
            *s = x as i16;
        }
        let got = overlap.emit(&span, v[0] as i16);
        let want: Vec<i16> = v[1 + SPAN..1 + SPAN + HOP]
            .iter()
            .map(|&x| x as i16)
            .collect();
        assert_eq!(got.to_vec(), want, "line {i}");
    }
    text.lines().count()
}

#[test]
fn overlap_sequence() {
    let path = format!("{}/tests/data/shaping.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} spans verified", replay(&path));
}
