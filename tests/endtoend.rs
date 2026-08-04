//! End-to-end check against the reference decoder's own output.

use mcelp::{Decoder, bitstream};

fn decode_all(bits: &str) -> Vec<u8> {
    let mut decoder = Decoder::default();
    let mut out = Vec::new();
    for line in bits.lines() {
        let Some(frame) = bitstream::parse_hex_line(line) else {
            continue;
        };
        if let Some(pcm) = decoder.decode(&frame) {
            out.extend_from_slice(&pcm);
        }
    }
    out
}

#[test]
fn decodes_examples_bit_exactly() {
    for name in ["male", "female"] {
        let bits = std::fs::read_to_string(format!("examples/{name}.mcelp")).unwrap();
        let want = std::fs::read(format!("examples/{name}.decoded.ulaw")).unwrap();
        let got = decode_all(&bits);
        assert_eq!(got.len(), want.len(), "{name}: length");
        if got != want {
            let first = got.iter().zip(&want).position(|(a, b)| a != b).unwrap();
            panic!(
                "{name}: first mismatch at byte {first} (frame {}, sample {})",
                first / 320,
                first % 320
            );
        }
    }
}

/// The examples contain no erased frames, so build some: setting the
/// suppression bit switches the decoder onto its concealment path.
#[test]
fn concealment_is_self_consistent() {
    for name in ["male", "female"] {
        let bits = std::fs::read_to_string(format!("examples/{name}.mcelp")).unwrap();
        let mut erased = String::new();
        for (i, line) in bits.lines().take(600).enumerate() {
            let mut frame = bitstream::parse_hex_line(line).unwrap();
            // Two bursts, one long enough to reach the deep-attenuation branch.
            if (100..115).contains(&i) || (200..230).contains(&i) {
                frame[17] |= 0x20;
            }
            for b in frame {
                erased.push_str(&format!("{b:02x}"));
            }
            erased.push('\n');
        }
        let out = decode_all(&erased);
        assert_eq!(out.len(), 75 * 320);
    }
}
