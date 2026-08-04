//! The closed-loop search bracket.
//!
//! Columns: the open-loop lag less three, then the first and last lag.

use mcelp::pitch_search::search_range;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(search_range(v[0] + 3), (v[1], v[2]), "line {i}");
    }
    text.lines().count()
}

#[test]
fn search_range_excerpt() {
    assert_eq!(replay("tests/data/search_range.txt"), 120);
}
