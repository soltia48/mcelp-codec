//! The gain codebook search.
//!
//! Columns: the adaptive scale, five mantissas, five shifts, the gain shift,
//! then the chosen index and the two gains it stands for.

use mcelp::gain::{MEASURES, search};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mantissa: [i16; MEASURES] = v[1..1 + MEASURES].try_into().unwrap();
        let shifts: [i16; MEASURES] = v[1 + MEASURES..1 + 2 * MEASURES].try_into().unwrap();
        let got = search(v[0], &mantissa, &shifts, v[1 + 2 * MEASURES]);
        let base = 2 + 2 * MEASURES;
        assert_eq!(
            (got.0 as i16, got.1, got.2),
            (v[base], v[base + 1], v[base + 2]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn search_excerpt() {
    assert_eq!(replay("tests/data/gain_search.txt"), 120);
}
