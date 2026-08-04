//! The encoder against the reference, frame by frame, as packed bitstream.
//!
//! Each expected line is one frame of nine hexadecimal words, which is what
//! the reference encoder writes.

use mcelp::bitstream::pack;
use mcelp::encoder::Encoder;
use mcelp::preprocess::FRAME;

fn line(words: &[i16; 9]) -> String {
    words.iter().map(|&w| format!("{:04x}", w as u16)).collect()
}

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
        let Some(expected) = want.next() else {
            return n;
        };
        assert_eq!(line(&pack(&params)), expected, "frame {n}");
        n += 1;
    }
    n
}

#[test]
fn encode_endtoend_examples() {
    let dir = env!("CARGO_MANIFEST_DIR");
    for example in ["male", "female"] {
        assert_eq!(
            replay(
                &format!("{dir}/examples/{example}.orig.ulaw"),
                &format!("{dir}/tests/data/encode_endtoend_{example}.txt"),
            ),
            75,
            "{example}"
        );
    }
}
