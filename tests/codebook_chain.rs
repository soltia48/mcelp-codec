//! The fixed codebook's front half composed: correlations, shortlist, gather.
//!
//! Columns: the mode, eighty target samples and eighty of the adaptive
//! contribution, the track responses, the gain and its flag, then the three
//! gathered arrays.  Only two of the three are checked — see below.

use mcelp::pulses::{correlations, gather, kept, shortlist, tracks};

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    // The adaptive correlations live in memory the reference reuses, so they
    // have to be carried from one call to the next.
    let mut adaptive = [0i16; SPAN];
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mode = v[0] as usize;
        let n = SPAN * tracks(mode);
        let responses = &v[1 + 2 * SPAN..1 + 2 * SPAN + n];
        let (gain, active) = (v[1 + 2 * SPAN + n], v[2 + 2 * SPAN + n]);

        let (value, sign, _) = correlations(
            mode,
            &v[1..1 + SPAN],
            &v[1 + SPAN..1 + 2 * SPAN],
            responses,
            gain,
            active,
            &mut adaptive,
        );
        let list = shortlist(mode, &value);
        let (values, shares, signs) = gather(mode, &list.position, &value, &adaptive, &sign);

        let total: usize = (0..tracks(mode)).map(|t| kept(mode, t)).sum();
        let at = 3 + 2 * SPAN + n;
        assert_eq!(values.as_slice(), &v[at..at + total], "line {i} values");
        assert_eq!(
            signs.as_slice(),
            &v[at + 2 * total..at + 3 * total],
            "line {i} signs"
        );
        // The adaptive correlations are only written when the adaptive codebook
        // contributed; otherwise the array holds whatever the shape vectors
        // the shape search left there, and nothing downstream reads it.
        if active > 0 {
            assert_eq!(
                shares.as_slice(),
                &v[at + total..at + 2 * total],
                "line {i} shares"
            );
        }
    }
    text.lines().count()
}

#[test]
fn codebook_chain_excerpt() {
    assert_eq!(replay("tests/data/codebook_chain.txt"), 60);
}
