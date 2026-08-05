//! The fixed-codebook search composed: the five pulse modes' scores and
//! the index chosen.
//!
//! Columns: the weight row, the lag and its fraction, the subframe's coding
//! mode, the two extension gains, then the target, the adaptive codebook's
//! filtered contribution, the impulse response and the carried per-position
//! correlations (eighty each), then the five modes' score pairs and the
//! codebook index.

use mcelp::pulses::{Search, search};

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut adaptive = [0i16; SPAN];
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let at = 6;
        adaptive.copy_from_slice(&v[at + 3 * SPAN..at + 4 * SPAN]);
        let input = Search {
            target: &v[at..at + SPAN],
            contribution: &v[at + SPAN..at + 2 * SPAN],
            impulse: &v[at + 2 * SPAN..at + 3 * SPAN],
            lag: v[1],
            fraction: v[2],
            extension: (v[4], v[5]),
            coding: v[3],
            weights: v[0] as usize,
        };
        let got = search(&input, &mut adaptive);
        let want: Vec<(i16, i16)> = (0..5)
            .map(|m| (v[at + 4 * SPAN + 2 * m], v[at + 4 * SPAN + 2 * m + 1]))
            .collect();
        assert_eq!(got.scores.as_slice(), want.as_slice(), "line {i} scores");
        assert_eq!(
            got.numbers.as_slice(),
            &v[at + 4 * SPAN + 13..at + 4 * SPAN + 29],
            "line {i} numbers"
        );
        assert_eq!(
            got.pair,
            (v[at + 4 * SPAN + 11], v[at + 4 * SPAN + 12]),
            "line {i} pair"
        );
        assert_eq!(got.index, v[at + 4 * SPAN + 10], "line {i} index");
    }
    text.lines().count()
}

#[test]
fn fixed_codebook_search_excerpt() {
    assert_eq!(replay("tests/data/fixed_codebook_search.txt"), 120);
}
