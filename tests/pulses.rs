//! The per-track shortlist of pulse positions.
//!
//! Columns: the mode, the eighty correlations, then the sixty shortlisted
//! positions and their sixty ranks.

use mcelp::pulses::{TRACK_STRIDE, TRACKS, kept, shortlist, tracks};

const VALUES: usize = 80;
const LIST: usize = TRACKS * TRACK_STRIDE;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let got = shortlist(mode, &v[1..1 + VALUES]);
        let base = 1 + VALUES;
        // Only the entries the routine writes are meaningful; the rest of each
        // track's room keeps whatever the previous call left there.
        for track in 0..tracks(mode) {
            let at = TRACK_STRIDE * track;
            let n = kept(mode, track);
            assert_eq!(
                &got.position[at..at + n],
                &v[base + at..base + at + n],
                "line {i} track {track} positions"
            );
            assert_eq!(
                &got.rank[at..at + n],
                &v[base + LIST + at..base + LIST + at + n],
                "line {i} track {track} ranks"
            );
        }
    }
    text.lines().count()
}

#[test]
fn pulses_excerpt() {
    assert_eq!(replay("tests/data/pulses.txt"), 120);
}
