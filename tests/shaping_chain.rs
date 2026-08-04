//! The encoder's front end all the way to the shaping signal.
//!
//! From the raw mu-law input through the input high-pass, the noise
//! suppression and the three shaping filters, against what the reference
//! leaves in its shaping buffer at the end of each half-frame.

use mcelp::frontend::Frontend;
use mcelp::preprocess::{FRAME, HALF};
use mcelp::shaping::HOP;
use mcelp::weighting::{Highpass, Smooth, Tilt};

/// Everything the shaping path does, in order.
#[derive(Default)]
struct Shaping {
    highpass: Highpass,
    smooth: Smooth,
    tilt: Tilt,
}

impl Shaping {
    fn run(&mut self, block: &mut [i16]) {
        self.highpass.run(block);
        self.smooth.run(block);
        self.tilt.run(block);
    }
}

fn replay(ulaw: &str, expected: &str) -> usize {
    let pcm = std::fs::read(ulaw).unwrap_or_else(|e| panic!("{ulaw}: {e}"));
    let text = std::fs::read_to_string(expected).unwrap_or_else(|e| panic!("{expected}: {e}"));
    let mut want = text.lines();
    let (mut frontend, mut shaping) = (Frontend::default(), Shaping::default());
    let mut n = 0;

    for chunk in pcm.chunks_exact(FRAME) {
        let mut frame = [0u8; FRAME];
        frame.copy_from_slice(chunk);
        let linear = frontend.condition(&frame);
        for half in 0..2 {
            let mut block = [0i16; HALF];
            block.copy_from_slice(&linear[half * HALF..(half + 1) * HALF]);
            let mut got = frontend.process(&block).to_vec();
            shaping.run(&mut got);

            let Some(line) = want.next() else { return n };
            let expect: Vec<i16> = line
                .split_whitespace()
                .map(|x| x.parse::<i64>().unwrap() as i16)
                .collect();
            assert_eq!(got.len(), HOP);
            assert_eq!(got, expect, "half-frame {n}");
            n += 1;
        }
    }
    n
}

#[test]
fn shaping_chain_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/shaping_chain_{example}.txt"),
            ),
            150,
            "{example}"
        );
    }
}
