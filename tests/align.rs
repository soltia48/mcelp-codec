//! The five measures brought to one exponent.
//!
//! Columns: five mantissas, five exponents, the pitch gain's exponent, then the
//! shifted mantissas and the shifts left over.

use mcelp::gain::{MEASURES, align};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut mantissa: [i16; MEASURES] = v[..MEASURES].try_into().unwrap();
        let exponent: [i16; MEASURES] = v[MEASURES..2 * MEASURES].try_into().unwrap();
        let shifts = align(&mut mantissa, &exponent, v[2 * MEASURES]);
        let base = 2 * MEASURES + 1;
        assert_eq!(
            mantissa.as_slice(),
            &v[base..base + MEASURES],
            "line {i} mantissa"
        );
        assert_eq!(
            shifts.as_slice(),
            &v[base + MEASURES..base + 2 * MEASURES],
            "line {i} shifts"
        );
    }
    text.lines().count()
}

#[test]
fn align_excerpt() {
    assert_eq!(replay("tests/data/align.txt"), 120);
}
