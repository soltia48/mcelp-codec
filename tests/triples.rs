//! The joint search over shortlisted triples.
//!
//! Columns: the three shortlist sizes, the correlations, the self energies, the
//! pair energies, the sign flags, then the three chosen indices and used flags.

use mcelp::pulses::search_triples;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let keep = [v[0] as usize, v[1] as usize, v[2] as usize];
        let total = keep[0] + keep[1] + keep[2];
        let pairs = keep[0] * (keep[1] + keep[2]) + keep[1] * keep[2];
        let mut at = 3;
        let correlation = &v[at..at + total];
        at += total;
        let energy = &v[at..at + total];
        at += total;
        let cross = &v[at..at + pairs];
        at += pairs;
        let sign = &v[at..at + total];
        at += total;

        let got = search_triples(keep, correlation, energy, cross, sign);
        assert_eq!(got.index, v[at..at + 3], "line {i} index");
        assert_eq!(got.used, v[at + 3..at + 6], "line {i} used");
    }
    text.lines().count()
}

#[test]
fn triples_excerpt() {
    assert_eq!(replay("tests/data/triples.txt"), 120);
}
