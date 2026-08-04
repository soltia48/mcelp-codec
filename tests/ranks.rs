//! Chosen shortlist slots turned into track positions.
//!
//! Columns: the mode, three chosen slots, sixty recorded ranks, then the three
//! positions.

use mcelp::pulses::ranks;

const LIST: usize = 60;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut chosen: [i16; 3] = v[1..4].try_into().unwrap();
        ranks(v[0] as usize, &mut chosen, &v[4..4 + LIST]);
        assert_eq!(chosen.as_slice(), &v[4 + LIST..], "line {i}");
    }
    text.lines().count()
}

#[test]
fn ranks_excerpt() {
    assert_eq!(replay("tests/data/ranks.txt"), 120);
}
