//! The point that halves the subframe's rectified area.

use mcelp::SUBFRAME;
use mcelp::pulses::centroid;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(centroid(&v[..SUBFRAME]), v[SUBFRAME], "line {i}");
    }
    text.lines().count()
}

#[test]
fn centroid_excerpt() {
    assert_eq!(replay("tests/data/centroid.txt"), 120);
}
