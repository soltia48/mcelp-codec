//! The joint search over shortlisted pairs.
//!
//! Columns: the two shortlist sizes, the correlations, the self energies, the
//! pair energies, the sign flags, then the three chosen indices and three used
//! flags.

use mcelp::pulses::search_pairs;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let keep = [v[0] as usize, v[1] as usize];
        let (both, pairs) = (keep[0] + keep[1], keep[0] * keep[1]);
        let mut at = 2;
        let correlation = &v[at..at + both];
        at += both;
        let energy = &v[at..at + both];
        at += both;
        let cross = &v[at..at + pairs];
        at += pairs;
        let sign = &v[at..at + both];
        at += both;

        let got = search_pairs(keep, correlation, energy, cross, sign);
        assert_eq!(got.index[..2], v[at..at + 2], "line {i} index");
        assert_eq!(got.used[..2], v[at + 3..at + 5], "line {i} used");
    }
    text.lines().count()
}

#[test]
fn pairs_excerpt() {
    assert_eq!(replay("tests/data/pairs.txt"), 120);
}
