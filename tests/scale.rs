//! Per-band smoothing of the noise scale.
//!
//! Columns: the voice-activity flag, the twenty scales this frame computed, the
//! twenty the smoother was holding, then both arrays after the update.

use mcelp::bands::BANDS;
use mcelp::shaping::Scale;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut smoother = Scale::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let at = |n: usize| &v[1 + n * BANDS..1 + (n + 1) * BANDS];
        let mut scale = [0i16; BANDS];
        scale.copy_from_slice(at(0));
        smoother.smooth(&mut scale, v[0] != 0);
        assert_eq!(scale.to_vec(), at(2).to_vec(), "line {i}");
        assert_eq!(scale.to_vec(), at(3).to_vec(), "carry, line {i}");
    }
    text.lines().count()
}

#[test]
fn scale_sequence() {
    let path = format!("{}/tests/data/scale.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
