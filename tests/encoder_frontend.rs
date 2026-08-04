//! The encoder's analysis front end from mu-law input to shaping signal.
//!
//! The reference output was dumped from the instrumented reference at the end of
//! the analysis front end; the input is the raw example file, so this exercises everything
//! from the expansion table through to the overlap.

use mcelp::frontend::Frontend;
use mcelp::preprocess::{FRAME, HALF};
use mcelp::shaping::HOP;

fn replay(ulaw: &str, expected: &str) -> usize {
    let pcm = std::fs::read(ulaw).unwrap_or_else(|e| panic!("{ulaw}: {e}"));
    let text = std::fs::read_to_string(expected).unwrap_or_else(|e| panic!("{expected}: {e}"));
    let mut want = text.lines();
    let mut frontend = Frontend::default();
    let mut n = 0;

    for chunk in pcm.chunks_exact(FRAME) {
        let mut frame = [0u8; FRAME];
        frame.copy_from_slice(chunk);
        let linear = frontend.condition(&frame);
        for half in 0..2 {
            let mut block = [0i16; HALF];
            block.copy_from_slice(&linear[half * HALF..(half + 1) * HALF]);
            let got = frontend.process(&block);
            let Some(line) = want.next() else { return n };
            let expect: Vec<i16> = line
                .split_whitespace()
                .map(|x| x.parse::<i64>().unwrap() as i16)
                .collect();
            assert_eq!(got.to_vec(), expect[..HOP].to_vec(), "half-frame {n}");
            n += 1;
        }
    }
    n
}

#[test]
fn frontend_sequence() {
    let root = env!("CARGO_MANIFEST_DIR");
    let ulaw = format!("{root}/examples/male.orig.ulaw");
    let out = format!("{root}/tests/data/frontend.txt");
    println!("{} half-frames verified", replay(&ulaw, &out));
}
