//! A pulse shape filtered by the weighted impulse response.
//!
//! Columns: the shape's number in the shape table, 80 impulse-response
//! samples, then the 80 filtered samples the reference produced.

use mcelp::convolve::{SHAPE, SUBFRAME, filter};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let shape = &mcelp::tables::FCB_SHAPES[SHAPE * v[0] as usize..][..SHAPE];
        let response: Vec<i16> = v[1..1 + SUBFRAME].iter().map(|&x| x as i16).collect();
        let want: Vec<i16> = v[1 + SUBFRAME..].iter().map(|&x| x as i16).collect();
        assert_eq!(filter(shape, &response).to_vec(), want, "line {i}");
    }
    text.lines().count()
}

#[test]
fn convolve_sequence() {
    let path = format!("{}/tests/data/convolve.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} shapes verified", replay(&path));
}
