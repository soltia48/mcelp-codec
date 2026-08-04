//! The per-band weights and scales.
//!
//! Columns: the voice-activity flag, the overall noise level, the level sum,
//! the twenty mix values, the twenty band levels, then the twenty weights and
//! twenty scales the reference encoder produced.  Runs are chained.

use mcelp::bands::BANDS;
use mcelp::weights::Shaping;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut state = Shaping::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut mix = [0i64; BANDS];
        mix.copy_from_slice(&v[3..3 + BANDS]);
        let mut level = [0i16; BANDS];
        for (l, &x) in level.iter_mut().zip(&v[3 + BANDS..]) {
            *l = x as i16;
        }
        let (weight, scale) = state.finish(&mut mix, &level, v[1] as i16, v[2] as i16, v[0] != 0);
        let want = |n: usize| -> Vec<i16> {
            v[3 + (2 + n) * BANDS..3 + (3 + n) * BANDS]
                .iter()
                .map(|&x| x as i16)
                .collect()
        };
        assert_eq!(weight.to_vec(), want(0), "weight, line {i}");
        assert_eq!(scale.to_vec(), want(1), "scale, line {i}");
    }
    text.lines().count()
}

#[test]
fn shape_sequence() {
    let path = format!("{}/tests/data/shape.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
