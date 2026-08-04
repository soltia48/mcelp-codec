//! The twelve fields the subframe searches contribute, per frame.
//!
//! Columns: lag, codebook index and gain index for each of the four subframes,
//! in coding order.

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
        let v: Vec<i16> = line
            .split_whitespace()
            .map(|x| x.parse::<i32>().unwrap() as i16)
            .collect();
        assert_eq!(&params.field[2..], v.as_slice(), "frame {n}");
        n += 1;
    }
    n
}

#[test]
fn encoder_fields_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/encoder_fields_{example}.txt"),
            ),
            75,
            "{example}"
        );
    }
}
