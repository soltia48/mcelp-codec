//! The leaky integrator the analysis runs over the inverse transform's output.
//!
//! Each line is one span: 256 samples in, then the 256 the reference encoder
//! produced.  The runs are chained so the recursion's carry-over is checked.

use mcelp::shaping::{SPAN, Smoother};

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut smoother = Smoother::default();
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i64>().unwrap() as i16)
            .collect();
        let mut span = [0i16; SPAN];
        span.copy_from_slice(&v[..SPAN]);
        smoother.run(&mut span);
        assert_eq!(span.to_vec(), v[SPAN..].to_vec(), "line {i}");
    }
    text.lines().count()
}

#[test]
fn integrator_sequence() {
    let path = format!("{}/tests/data/smoother.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} spans verified", replay(&path));
}
