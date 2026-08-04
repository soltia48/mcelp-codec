//! The encoder's line spectrum fields, from the raw mu-law input.
//!
//! Runs the frame loop over a whole example and checks the two fields the
//! quantiser contributes against what the reference wrote.

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
        let params = encoder.frame(&frame);
        let Some(line) = want.next() else { return n };
        let expect: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!([params.field[0], params.field[1]], expect[..2], "frame {n}");
        n += 1;
    }
    n
}

#[test]
fn encoder_lsf_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/encoder_lsf_{example}.txt"),
            ),
            75,
            "{example}"
        );
    }
}
