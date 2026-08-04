//! The fixed codebook's back half composed: energies, correction, joint search.
//!
//! Columns: the mode, the shortlist sizes, the track responses, the shortlist
//! positions, the gain and its shift, the correlations, the adaptive
//! correlations, the signs, then the three chosen positions and three flags.

use mcelp::pulses::{TRACK_STRIDE, correct, energies, kept, search_pairs, search_triples, tracks};

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let n = tracks(mode);
        let keep: Vec<usize> = (0..n).map(|t| kept(mode, t)).collect();
        let total: usize = keep.iter().sum();

        let mut at = 1 + n;
        let mut responses = v[at..at + SPAN * n].to_vec();
        at += SPAN * n;
        let positions = &v[at..at + TRACK_STRIDE * n];
        at += TRACK_STRIDE * n;
        let (gain, gain_shift) = (v[at], v[at + 1]);
        at += 2;
        let correlation = &v[at..at + total];
        at += total;
        let adaptive = &v[at..at + total];
        at += total;
        let sign = &v[at..at + total];
        at += total;

        let (headroom, mut pairs, mut alone) = energies(&keep, &mut responses, positions);
        correct(
            &keep, &mut pairs, &mut alone, adaptive, sign, gain, gain_shift, headroom,
        );
        let got = if n == 3 {
            search_triples(
                [keep[0], keep[1], keep[2]],
                correlation,
                &alone,
                &pairs,
                sign,
            )
        } else {
            search_pairs([keep[0], keep[1]], correlation, &alone, &pairs, sign)
        };
        assert_eq!(got.index[..n], v[at..at + n], "line {i} index");
        assert_eq!(got.used[..n], v[at + 3..at + 3 + n], "line {i} used");
    }
    text.lines().count()
}

#[test]
fn codebook_back_excerpt() {
    assert_eq!(replay("tests/data/codebook_back.txt"), 60);
}
