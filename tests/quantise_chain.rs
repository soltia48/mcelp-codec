//! The line spectrum quantiser composed: cosine domain in, transport fields out.
//!
//! The line spectrum arrives in the cosine domain, so it has to be taken back
//! to frequencies before the quantiser sees it.  The predictor memory is the
//! encoder's own copy, which is a different array from the decoder's.
//!
//! Columns: ten line spectral pairs, ten spare words, thirty words of predictor
//! memory, then the two transport fields.

use mcelp::lpc::ORDER;
use mcelp::lsf_weight::{quantise, weights};
use mcelp::lsp::lsp_to_lsf;

fn replay(path: &str, own_weights: bool) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let lsp: [i16; ORDER] = v[..ORDER].try_into().unwrap();
        let lsf = lsp_to_lsf(&lsp);
        let traced: [i16; ORDER] = v[ORDER..2 * ORDER].try_into().unwrap();
        let weight = if own_weights { weights(&lsf) } else { traced };

        let mut history = [[0i16; ORDER]; 3];
        for (j, row) in history.iter_mut().enumerate() {
            let at = 2 * ORDER + j * ORDER;
            row.copy_from_slice(&v[at..at + ORDER]);
        }
        let got = quantise(&lsf, &weight, &history).fields();
        assert_eq!([got.0, got.1], v[2 * ORDER + 3 * ORDER..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn quantise_chain_excerpt() {
    assert_eq!(replay("tests/data/quantise_chain.txt", true), 120);
}
