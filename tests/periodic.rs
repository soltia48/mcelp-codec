//! Is the band around the lag strongly periodic?
//!
//! Columns: the lag, the base offset, 64 open-loop correlations, then the
//! reference's verdict.

use mcelp::pitch_search::periodic;

fn replay(path: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let got = periodic(v[0] as i16, v[1] as i16, &v[2..66]);
        assert_eq!(got as i64, v[66], "line {i}");
    }
    text.lines().count()
}

#[test]
fn periodic_sequence() {
    let path = format!("{}/tests/data/periodic.txt", env!("CARGO_MANIFEST_DIR"));
    println!("{} decisions verified", replay(&path));
}
