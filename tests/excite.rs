//! The subframe's excitation from the two codebooks.
//!
//! Columns: the adaptive and fixed gains, the eighty innovation samples, the
//! eighty excitation samples going in and the eighty coming out.

use mcelp::excitation::excite;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut excitation = v[2 + SPAN..2 + 2 * SPAN].to_vec();
        excite(&mut excitation, &v[2..2 + SPAN], v[0], v[1]);
        assert_eq!(excitation.as_slice(), &v[2 + 2 * SPAN..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn excite_excerpt() {
    assert_eq!(replay("tests/data/excite.txt"), 120);
}
