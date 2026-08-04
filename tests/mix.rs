//! The per-band mix of two exponentials.
//!
//! Columns: the level sum, the overall noise level, the voice-activity score
//! and flag, the twenty band levels, then the twenty 32-bit mix values the
//! reference encoder produced.

use mcelp::bands::BANDS;
use mcelp::weights::Mix;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut level = [0i16; BANDS];
        for (l, &x) in level.iter_mut().zip(&v[4..4 + BANDS]) {
            *l = x as i16;
        }
        let (got, _) = Mix::run(&level, v[0] as i16, v[1] as i16, v[2] as i16, v[3] != 0);
        assert_eq!(got.to_vec(), v[4 + BANDS..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn mix_sequence() {
    let path = format!("{}/tests/data/mix.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
