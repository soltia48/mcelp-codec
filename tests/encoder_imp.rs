//! The impulse response each of the four subframes convolves with.

use mcelp::encoder::Encoder;
use mcelp::preprocess::FRAME;

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
            for sub in 0..2 {
                let k = 2 * half + sub;
                assert_eq!(
                    encoder.impulse()[half][sub].as_slice(),
                    &v[k * 80..(k + 1) * 80],
                    "frame {n} half {half} subframe {sub}"
                );
            }
        }
        n += 1;
    }
    n
}

#[test]
fn encoder_imp_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/encoder_imp_{example}.txt"),
            ),
            75,
            "{example}"
        );
    }
}
