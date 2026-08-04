//! The library surface: the byte-level round trip a caller actually uses.
//!
//! Everything here goes through the crate root rather than the stage modules,
//! so it also pins the entry points down against accidental renaming.

use mcelp::{Decoder, Encoder, FRAME, FRAME_BYTES, bitstream};

fn ulaw(name: &str) -> Vec<u8> {
    let dir = env!("CARGO_MANIFEST_DIR");
    std::fs::read(format!("{dir}/examples/{name}.orig.ulaw")).unwrap()
}

fn frames(name: &str) -> Vec<[u8; FRAME_BYTES]> {
    let dir = env!("CARGO_MANIFEST_DIR");
    let text = std::fs::read_to_string(format!("{dir}/examples/{name}.mcelp")).unwrap();
    text.lines().filter_map(bitstream::parse_hex_line).collect()
}

/// `Encoder::encode` must agree with the `.mcelp` the reference encoder wrote,
/// which is the same thing the hex container holds.
#[test]
fn encode_produces_the_reference_frames() {
    for name in ["male", "female"] {
        let pcm = ulaw(name);
        let want = frames(name);
        let mut encoder = Encoder::new();
        let mut n = 0;
        for (chunk, expected) in pcm.chunks_exact(FRAME).zip(&want) {
            let got = encoder.encode(chunk.try_into().unwrap());
            assert_eq!(&got, expected, "{name}: frame {n}");
            n += 1;
        }
        assert_eq!(n, 75, "{name}");
    }
}

/// And the bytes it produces must be exactly what the container round-trips.
#[test]
fn hex_container_round_trips() {
    for name in ["male", "female"] {
        for frame in frames(name) {
            let line = bitstream::to_hex_line(&frame);
            assert_eq!(line.len(), 2 * FRAME_BYTES);
            assert_eq!(bitstream::parse_hex_line(&line), Some(frame));
        }
    }
}

/// `to_bytes` is the inverse of `canonicalize` for frames this crate packs.
#[test]
fn packing_round_trips_through_bytes() {
    for name in ["male", "female"] {
        for frame in frames(name) {
            let words = bitstream::canonicalize(&frame);
            assert_eq!(bitstream::canonicalize(&bitstream::to_bytes(&words)), words);
        }
    }
}

/// Encoding then decoding a whole example must reproduce the bundled output,
/// using only the root API.
#[test]
fn encode_then_decode_matches_the_reference() {
    for name in ["male", "female"] {
        let dir = env!("CARGO_MANIFEST_DIR");
        let want = std::fs::read(format!("{dir}/examples/{name}.decoded.ulaw")).unwrap();

        let (mut encoder, mut decoder) = (Encoder::new(), Decoder::new());
        let mut out = Vec::with_capacity(want.len());
        for chunk in ulaw(name).chunks_exact(FRAME) {
            let frame = encoder.encode(chunk.try_into().unwrap());
            if let Some(pcm) = decoder.decode(&frame) {
                out.extend_from_slice(&pcm);
            }
        }
        assert_eq!(out, want, "{name}");
    }
}

/// The linear form differs from the mu-law one only by the output companding.
#[test]
fn decode_linear_agrees_with_decode() {
    let (mut a, mut b) = (Decoder::new(), Decoder::new());
    for frame in frames("male") {
        let companded = a.decode(&frame).unwrap();
        let linear = b.decode_linear(&frame).unwrap();
        let mut expected = [0u8; FRAME];
        mcelp::ulaw::frame_from_linear(&linear, &mut expected);
        assert_eq!(companded, expected);
    }
}
