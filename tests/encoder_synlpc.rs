//! The synthesis filter of the first half-frame's two subframes.

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
        for sub in 0..2 {
            assert_eq!(
                encoder.synthesis_lpc()[0][sub].as_slice(),
                &v[sub * N..(sub + 1) * N],
                "frame {n} subframe {sub}"
            );
        }
        n += 1;
    }
    n
}

#[test]
fn encoder_synlpc_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/encoder_synlpc_{example}.txt"),
            ),
            75,
            "{example}"
        );
    }
}
