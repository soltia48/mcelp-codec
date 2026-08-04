//! The weighted distortion the LSF quantiser minimises.
//!
//! Columns: the predictor's mode, ten target values, ten candidate values,
//! ten weights, then the distortion the reference computed.

use mcelp::lpc::ORDER;
use mcelp::lsf_weight::distortion;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let at = |n: usize| -> [i16; ORDER] {
            let mut a = [0i16; ORDER];
            for (s, &x) in a.iter_mut().zip(&v[1 + n * ORDER..]) {
                *s = x as i16;
            }
            a
        };
        let scale: &[i16; ORDER] = mcelp::tables::LSF_MA_GAIN[ORDER * v[0] as usize..][..ORDER]
            .try_into()
            .unwrap();
        assert_eq!(
            distortion(&at(0), &at(1), scale, &at(2)),
            v[1 + 3 * ORDER],
            "line {i}"
        );
    }
    text.lines().count()
}

#[test]
fn distortion_sequence() {
    let path = format!("{}/tests/data/distortion.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} evaluations verified", replay(&path));
}
