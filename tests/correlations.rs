//! The target's correlation at every candidate pulse position.
//!
//! Columns: the mode, the target and the adaptive contribution (eighty words
//! each), the filtered responses, the gain and its flag, then the normalised
//! magnitudes, the signs, the adaptive correlations and the shift.

use mcelp::pulses::{correlations, grid};

const POSITIONS: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let tracks = if mode < 2 { 3 } else { 2 };
        let mut at = 1;
        let target = &v[at..at + POSITIONS];
        let contribution = &v[at + POSITIONS..at + 2 * POSITIONS];
        at += 2 * POSITIONS;
        let responses = &v[at..at + POSITIONS * tracks];
        at += POSITIONS * tracks;
        let (gain, active) = (v[at], v[at + 1]);
        at += 2;

        let mut adaptive = [0i16; POSITIONS];
        let (value, sign, up) = correlations(
            mode,
            target,
            contribution,
            responses,
            gain,
            active,
            &mut adaptive,
        );
        assert_eq!(value.as_slice(), &v[at..at + POSITIONS], "line {i} values");
        at += POSITIONS;
        assert_eq!(sign.as_slice(), &v[at..at + POSITIONS], "line {i} signs");
        at += POSITIONS;
        // Only the mode's own grid is written, and only when the adaptive
        // codebook contributed; the rest keeps the previous call's values.
        if active > 0 {
            let (stride, count) = grid(mode);
            for k in 0..count {
                assert_eq!(
                    adaptive[k * stride],
                    v[at + k * stride],
                    "line {i} adaptive {k}"
                );
            }
        }
        at += POSITIONS;
        assert_eq!(up, v[at], "line {i} shift");
    }
    text.lines().count()
}

#[test]
fn correlations_excerpt() {
    assert_eq!(replay("tests/data/correlations.txt"), 60);
}
