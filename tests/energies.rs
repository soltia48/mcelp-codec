//! The pulse response energies.
//!
//! Columns: the track count, the shortlist sizes, the filtered responses, the
//! shortlisted positions, then the reported shift, the scaled responses, the
//! pair energies and the self energies.

use mcelp::pulses::{TRACK_STRIDE, energies};

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let tracks = v[0] as usize;
        let keep: Vec<usize> = v[1..1 + tracks].iter().map(|&x| x as usize).collect();
        let mut at = 1 + tracks;
        let mut filtered = v[at..at + SPAN * tracks].to_vec();
        at += SPAN * tracks;
        let position = &v[at..at + TRACK_STRIDE * tracks];
        at += TRACK_STRIDE * tracks;

        let (shift, pairs, alone) = energies(&keep, &mut filtered, position);
        assert_eq!(shift, v[at], "line {i} shift");
        at += 1;
        assert_eq!(
            filtered.as_slice(),
            &v[at..at + SPAN * tracks],
            "line {i} responses"
        );
        at += SPAN * tracks;
        assert_eq!(pairs.as_slice(), &v[at..at + pairs.len()], "line {i} pairs");
        at += pairs.len();
        assert_eq!(
            alone.as_slice(),
            &v[at..at + alone.len()],
            "line {i} energies"
        );
    }
    text.lines().count()
}

#[test]
fn energies_excerpt() {
    assert_eq!(replay("tests/data/energies.txt"), 60);
}
