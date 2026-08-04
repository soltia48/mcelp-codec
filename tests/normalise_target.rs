//! The target scaled up for the search.

use mcelp::pulses::normalise_target;

const SPAN: usize = 80;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(
            normalise_target(&v[..SPAN]).as_slice(),
            &v[SPAN..],
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn normalise_target_excerpt() {
    assert_eq!(replay("tests/data/normalise_target.txt"), 120);
}
