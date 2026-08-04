//! Scaling the target before the correlations.
//!
//! Columns: 80 samples in, 80 out.

use mcelp::convolve::SUBFRAME;
use mcelp::pitch_search::prescale;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut target = [0i16; SUBFRAME];
        target.copy_from_slice(&v[..SUBFRAME]);
        prescale(&mut target);
        assert_eq!(target.to_vec(), v[SUBFRAME..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn prescale_sequence() {
    let path = format!("{}/tests/data/prescale.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} subframes verified", replay(&path));
}
