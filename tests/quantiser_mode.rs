//! Which mode the line spectrum is quantised in.
//!
//! Columns: eleven prediction coefficients, eighty target samples, ten line
//! spectral frequencies, the incoming state (two tracked spectra and the age),
//! then the chosen mode and the state on the way out.

use mcelp::lpc::ORDER;
use mcelp::lsp::Quantiser;

const SUBFRAME: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut state = Quantiser::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let lpc: [i16; ORDER + 1] = v[..ORDER + 1].try_into().unwrap();
        let target = &v[ORDER + 1..ORDER + 1 + SUBFRAME];
        let mut at = ORDER + 1 + SUBFRAME;
        let lsf: [i16; ORDER] = v[at..at + ORDER].try_into().unwrap();
        at += ORDER;
        state.restore(
            v[at..at + ORDER].try_into().unwrap(),
            v[at + ORDER..at + 2 * ORDER].try_into().unwrap(),
            v[at + 2 * ORDER],
        );
        at += 2 * ORDER + 1;

        let mode = state.choose(&lpc, target, &lsf);
        assert_eq!(mode, v[at], "line {i} mode");
        let (near, far, age) = state.saved();
        assert_eq!(near.as_slice(), &v[at + 1..at + 1 + ORDER], "line {i} near");
        assert_eq!(
            far.as_slice(),
            &v[at + 1 + ORDER..at + 1 + 2 * ORDER],
            "line {i} far"
        );
        assert_eq!(age, v[at + 1 + 2 * ORDER], "line {i} age");
    }
    text.lines().count()
}

#[test]
fn quantiser_excerpt() {
    assert_eq!(replay("tests/data/quantiser_mode.txt"), 120);
}
