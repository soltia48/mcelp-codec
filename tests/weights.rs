//! Per-band suppression depth and weight.
//!
//! Columns: the level sum, then the twenty 32-bit depths and twenty weights the
//! reference encoder laid down.  Each frame stands on its own.

use mcelp::bands::BANDS;
use mcelp::weights::depths;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let (depth, weight) = depths(v[0] as i16);
        assert_eq!(depth.to_vec(), v[1..1 + BANDS].to_vec(), "depth, line {i}");
        let want: Vec<i16> = v[1 + BANDS..].iter().map(|&x| x as i16).collect();
        assert_eq!(weight.to_vec(), want, "weight, line {i}");
    }
    text.lines().count()
}

#[test]
fn depth_sequence() {
    let path = format!("{}/tests/data/weights.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}

/// The overall noise level in dB.
///
/// Reuses the noise-floor fixture, whose second column is the level the
/// reference encoder handed down and whose band energies are the ones it was
/// computed from.
#[test]
fn overall_level_matches() {
    check_levels(&format!(
        "{}/tests/data/floor.txt",
        env!("CARGO_MANIFEST_DIR")
    ));
}

fn check_levels(path: &str) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut energy = [0i64; BANDS];
        energy.copy_from_slice(&v[3 + 129..3 + 129 + BANDS]);
        assert_eq!(
            mcelp::weights::overall_level(&energy) as i64,
            v[1],
            "line {i}"
        );
    }
    println!("{} frames verified", text.lines().count());
}
