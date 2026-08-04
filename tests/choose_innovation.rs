//! Pulse codebook or shape pair.
//!
//! Columns: the pulse codebook's pair, the shape pair's, the two shape numbers,
//! the index going in and the index coming out.

use mcelp::pulses::choose_innovation;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        let mut index = v[6];
        choose_innovation((v[0], v[1]), (v[2], v[3]), (v[4], v[5]), &mut index);
        assert_eq!(index, v[7], "line {i}");
    }
    text.lines().count()
}

#[test]
fn choose_innovation_excerpt() {
    assert_eq!(replay("tests/data/choose_innovation.txt"), 120);
}
