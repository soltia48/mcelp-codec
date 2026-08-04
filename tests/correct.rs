//! The adaptive share removed and signs folded in.
//!
//! Columns: the track count, the shortlist sizes, the gain, its shift, the
//! response headroom, the pair and self energies going in, the adaptive
//! correlations, the signs, then the pair and self energies coming out.

use mcelp::pulses::correct;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let tracks = v[0] as usize;
        let keep: Vec<usize> = v[1..1 + tracks].iter().map(|&x| x as usize).collect();
        let total: usize = keep.iter().sum();
        let count: usize = (0..tracks)
            .map(|t| keep[t] * keep[t + 1..].iter().sum::<usize>())
            .sum();
        let mut at = 1 + tracks;
        let (gain, gain_shift, headroom) = (v[at], v[at + 1], v[at + 2]);
        at += 3;
        let mut pairs = v[at..at + count].to_vec();
        at += count;
        let mut alone = v[at..at + total].to_vec();
        at += total;
        let adaptive = &v[at..at + total];
        at += total;
        let sign = &v[at..at + total];
        at += total;

        correct(
            &keep, &mut pairs, &mut alone, adaptive, sign, gain, gain_shift, headroom,
        );
        assert_eq!(pairs.as_slice(), &v[at..at + count], "line {i} pairs");
        at += count;
        assert_eq!(alone.as_slice(), &v[at..at + total], "line {i} energies");
    }
    text.lines().count()
}

#[test]
fn correct_excerpt() {
    assert_eq!(replay("tests/data/correct.txt"), 120);
}
