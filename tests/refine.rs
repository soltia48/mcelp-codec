//! Carry forward, clamp, fold into the mix.
//!
//! Columns: the voice-activity flag, the twenty depths, the twenty mix values,
//! the twenty tracked band energies, the twenty band levels, the twenty fixed
//! weights, then the depths and mix values the reference encoder ended up with.
//! Runs are chained so the carried-forward depths are checked too.

use mcelp::bands::BANDS;
use mcelp::weights::Depths;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut state = Depths::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let at = |n: usize| -> [i64; BANDS] {
            let mut a = [0i64; BANDS];
            a.copy_from_slice(&v[1 + n * BANDS..1 + (n + 1) * BANDS]);
            a
        };
        let words = |n: usize| -> [i16; BANDS] {
            let mut a = [0i16; BANDS];
            for (s, &x) in a.iter_mut().zip(&v[1 + n * BANDS..]) {
                *s = x as i16;
            }
            a
        };
        let (mut depth, mut mix) = (at(0), at(1));
        state.refine(
            &mut depth,
            &mut mix,
            &at(2),
            &words(3),
            &words(4),
            v[0] != 0,
        );
        assert_eq!(depth.to_vec(), at(5).to_vec(), "depth, line {i}");
        assert_eq!(mix.to_vec(), at(6).to_vec(), "mix, line {i}");
    }
    text.lines().count()
}

#[test]
fn refine_sequence() {
    let path = format!("{}/tests/data/refine.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
