//! The per-bin noise floor the suppression gain works from.
//!
//! Columns: headroom, the two level measures, 129 magnitude bins, the twenty
//! tracked band energies, the twenty per-band scales and weights, then the 129
//! floor bins the reference encoder produced.

use mcelp::analysis::BINS;
use mcelp::bands::BANDS;
use mcelp::shaping::noise_floor;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut magnitude = [0i16; BINS];
        for (m, &x) in magnitude.iter_mut().zip(&v[3..3 + BINS]) {
            *m = x as i16;
        }
        let mut energy = [0i64; BANDS];
        energy.copy_from_slice(&v[3 + BINS..3 + BINS + BANDS]);
        let word = |n: usize| -> [i16; BANDS] {
            let mut a = [0i16; BANDS];
            for (s, &x) in a.iter_mut().zip(&v[3 + BINS + n * BANDS..]) {
                *s = x as i16;
            }
            a
        };
        let got = noise_floor(
            &magnitude,
            &energy,
            &word(1),
            &word(2),
            v[0] as i16,
            v[1] as i16,
            v[2] as i16,
        );
        let want: Vec<i16> = v[3 + BINS + 3 * BANDS..]
            .iter()
            .map(|&x| x as i16)
            .collect();
        assert_eq!(got.to_vec(), want, "line {i}");
    }
    text.lines().count()
}

#[test]
fn floor_sequence() {
    let path = format!("{}/tests/data/floor.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
