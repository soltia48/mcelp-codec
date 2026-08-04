//! The weighting filter of each half-frame's first subframe.

use mcelp::encoder::Encoder;
use mcelp::lpc::ORDER;
use mcelp::preprocess::FRAME;

const N: usize = ORDER + 1;

fn replay(ulaw: &str, expected: &str) -> usize {
    let pcm = std::fs::read(ulaw).unwrap_or_else(|e| panic!("{ulaw}: {e}"));
    let text = std::fs::read_to_string(expected).unwrap_or_else(|e| panic!("{expected}: {e}"));
    let mut want = text.lines();
    let mut encoder = Encoder::default();
    let mut n = 0;
    for chunk in pcm.chunks_exact(FRAME) {
        let mut frame = [0u8; FRAME];
        frame.copy_from_slice(chunk);
        encoder.frame(&frame);
        let Some(line) = want.next() else { return n };
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        for half in 0..2 {
            let (num, den) = &encoder.weighting()[half][1];
            let at = half * 2 * N;
            assert_eq!(
                num.as_slice(),
                &v[at..at + N],
                "frame {n} half {half} numerator"
            );
            assert_eq!(
                den.as_slice(),
                &v[at + N..at + 2 * N],
                "frame {n} half {half} denominator"
            );
        }
        n += 1;
    }
    n
}

#[test]
fn encoder_weighting2_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/encoder_weighting2_{example}.txt"),
            ),
            75,
            "{example}"
        );
    }
}
