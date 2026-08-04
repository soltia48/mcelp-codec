//! Long-term band-energy reference: male.orig.ulaw.
//!
//! Each line of a fixture is one silent frame: the analysis headroom, the 129
//! magnitude bins, then the twenty 32-bit energies the reference encoder ended
//! up with.  Frames the detector called active are left out — the estimate
//! ignores those completely, so replaying only the silent ones is equivalent to
//! replaying the whole file.  The committed fixtures run long enough to walk
//! through all three stages of the adaptation ramp.

use mcelp::noise::NoiseEstimate;

fn replay(path: &str, mut estimate: NoiseEstimate) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for (i, line) in text.lines().enumerate() {
        let v: Vec<i64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut spectrum = [0i16; 129];
        for (s, &x) in spectrum.iter_mut().zip(&v[1..130]) {
            *s = x as i16;
        }
        estimate.update(&spectrum, v[0] as i16, false);
        assert_eq!(
            estimate.energy.to_vec(),
            v[130..150].to_vec(),
            "{path} line {i}"
        );
    }
    text.lines().count()
}

fn run(dir: &str) {
    let fast = replay(&format!("{dir}/noise_fast.txt"), NoiseEstimate::fast());
    let slow = replay(&format!("{dir}/noise_slow.txt"), NoiseEstimate::slow());
    println!("{fast} + {slow} frames verified");
}

#[test]
fn long_term_energy() {
    run(&format!("{}/tests/data", env!("CARGO_MANIFEST_DIR")));
}
