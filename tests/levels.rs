//! Per-band levels of the mean magnitude against the tracked noise.
//!
//! Columns: headroom, the voice-activity score, 129 magnitude bins, the twenty
//! slow band energies, then the twenty levels and the weighted sum the
//! reference encoder produced.  Runs are chained: the blend state and the
//! per-band smoothing memory both carry across frames.

use mcelp::analysis::BINS;
use mcelp::bands::BANDS;
use mcelp::shaping::Levels;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut levels = Levels::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut magnitude = [0i16; BINS];
        for (m, &x) in magnitude.iter_mut().zip(&v[2..2 + BINS]) {
            *m = x as i16;
        }
        let mut slow = [0i64; BANDS];
        slow.copy_from_slice(&v[2 + BINS..2 + BINS + BANDS]);
        let (got, sum) = levels.run(&magnitude, &slow, v[1] as i16, v[0] as i16);
        let want: Vec<i16> = v[2 + BINS + BANDS..2 + BINS + 2 * BANDS]
            .iter()
            .map(|&x| x as i16)
            .collect();
        assert_eq!(got.to_vec(), want, "line {i}");
        assert_eq!(sum as i64, v[2 + BINS + 2 * BANDS], "sum, line {i}");
    }
    text.lines().count()
}

#[test]
fn level_sequence() {
    let path = format!("{}/tests/data/levels.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} frames verified", replay(&path));
}
