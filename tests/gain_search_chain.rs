//! The gain search composed: five measures to the gain index.
//!
//! Columns: the adaptive scale, the gain exponent, five mantissas, five
//! exponents, then the chosen index and the two gains.

use mcelp::gain::{MEASURES, align, search};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut mantissa: [i16; MEASURES] = v[2..2 + MEASURES].try_into().unwrap();
        let exponent: [i16; MEASURES] = v[2 + MEASURES..2 + 2 * MEASURES].try_into().unwrap();

        let shifts = align(&mut mantissa, &exponent, v[1]);
        let (index, first, second) = search(v[0], &mantissa, &shifts, v[1]);

        let at = 2 + 2 * MEASURES;
        assert_eq!(
            (index as i16, first, second),
            (v[at], v[at + 1], v[at + 2]),
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn gain_search_chain_excerpt() {
    assert_eq!(replay("tests/data/gain_search_chain.txt"), 120);
}
