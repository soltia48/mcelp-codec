//! Spectrum in, detector features out.
//!
//! Columns: headroom, the voicing measure, 129 magnitude bins, the twenty fast
//! and twenty slow band energies, then the nine features the reference encoder
//! handed to the detector.

use mcelp::analysis::BINS;
use mcelp::bands::BANDS;
use mcelp::frontend::{features, unnormalise};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut magnitude = [0i16; BINS];
        for (m, &x) in magnitude.iter_mut().zip(&v[2..2 + BINS]) {
            *m = x as i16;
        }
        let mut fast = [0i64; BANDS];
        let mut slow = [0i64; BANDS];
        fast.copy_from_slice(&v[2 + BINS..2 + BINS + BANDS]);
        slow.copy_from_slice(&v[2 + BINS + BANDS..2 + BINS + 2 * BANDS]);

        let spectrum = unnormalise(&magnitude, v[0] as i16);
        let f = features(&spectrum, &fast, &slow, v[1] as i16);

        let w = &v[2 + BINS + 2 * BANDS..];
        assert_eq!(f.summary, w[0], "summary, line {i}");
        for (r, base) in [(0usize, 1usize), (1, 4)] {
            assert_eq!(f.reference[r].tracked, w[base], "tracked {r}, line {i}");
            assert_eq!(
                f.reference[r].difference,
                w[base + 1],
                "difference {r}, line {i}"
            );
            assert_eq!(
                f.reference[r].variance,
                w[base + 2],
                "variance {r}, line {i}"
            );
            assert_eq!(
                f.reference[r].excess_sum as i64,
                w[7 + r],
                "excess {r}, line {i}"
            );
        }
    }
    text.lines().count()
}

#[test]
fn feature_sequence() {
    let path = format!("{}/tests/data/features.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
