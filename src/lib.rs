//! M-CELP speech codec.
//!
//! A 3.6 kbit/s CELP codec: 320 mu-law samples in, 18 bytes out, and back.
//! Each stage is written as the algorithm it implements, and the routine label
//! quoted in its doc comment names the reference behaviour the test suite
//! checks it against.
//!
//! # Using it
//!
//! [`Encoder`] and [`Decoder`] are the entry points.  Both work one frame at a
//! time: [`FRAME`] mu-law samples in, [`FRAME_BYTES`] transport bytes out, and
//! back.  Both carry state between frames, so one instance must see a stream in
//! order; a partial frame at the end of a stream is dropped, as the reference
//! does.
//!
//! ```
//! use mcelp::{Decoder, Encoder, FRAME};
//!
//! let speech: Vec<u8> = vec![0xff; 4 * FRAME]; // mu-law, 8 kHz mono
//!
//! let (mut encoder, mut decoder) = (Encoder::new(), Decoder::new());
//! let mut out = Vec::new();
//! for block in speech.chunks_exact(FRAME) {
//!     let frame = encoder.encode(block.try_into().unwrap());
//!     assert_eq!(frame.len(), mcelp::FRAME_BYTES);
//!     if let Some(pcm) = decoder.decode(&frame) {
//!         out.extend_from_slice(&pcm);
//!     }
//! }
//! assert_eq!(out.len(), speech.len());
//! ```
//!
//! [`Decoder::decode_linear`] gives the same samples as 16-bit linear PCM, and
//! [`bitstream`] holds the frame's parameter fields and the hex container the
//! bundled examples use.
//!
//! Everything else is the stages the two are built from.  They are public so
//! that each can be replayed against the reference on its own, which is what
//! the test suite does; they are not a stable interface, and a caller that only
//! wants to code speech never needs them.
//!
//! # The codec
//!
//! The bit stream is very close to ITU-T G.729 in its parameter set and in most
//! of its arithmetic — the LSF quantiser's `GAP1`/`GAP2`/`GAP3` spacings, the
//! `L_LIMIT`/`M_LIMIT` clamps, the 1/3-resolution pitch lag coding and the
//! conjugate two-stage gain codebook are all recognisably G.729 — but it runs
//! four times slower: one frame is **320 samples (40 ms) split into four
//! 80-sample subframes**, carried in 18 bytes.
//!
//! ## Decoding ([`Decoder`])
//!
//! | stage | module |
//! |-------|--------|
//! | frame unpacking            | [`bitstream`] |
//! | LSF dequantisation, LSP/LPC conversion | [`lsp`] |
//! | pitch lag, adaptive codebook | [`pitch`] |
//! | fixed codebook             | [`codebook`] |
//! | gain decoding              | [`gain`] |
//! | short-term synthesis       | [`synth`] |
//! | postfiltering              | [`postfilter`], [`ltp`] |
//!
//! ## Encoding ([`Encoder`])
//!
//! | stage | module |
//! |-------|--------|
//! | input conditioning         | [`preprocess`] |
//! | spectral analysis          | [`analysis`], [`bands`] |
//! | voice activity, noise floor | [`vad`], [`noise`], [`frontend`] |
//! | noise suppression          | [`weights`], [`shaping`] |
//! | LPC analysis, LSF quantisation | [`lpc`], [`lsf_weight`] |
//! | perceptual weighting       | [`weighting`], [`weight_lpc`] |
//! | closed-loop pitch search   | [`pitch_search`], [`convolve`] |
//! | fixed-codebook search      | [`pulses`], [`mode`] |
//!
//! # Numerics
//!
//! [`fixed`] reproduces the three arithmetic properties the code depends on:
//! fractional
//! multiplies, 40-bit wrapping accumulators, and sign
//! extension of memory operands.  Everything else is ordinary arithmetic.

/// Samples per frame: 40 ms at 8 kHz.
pub const FRAME: usize = 320;
/// Samples per half-frame, which is the unit both halves of the codec work in.
pub const HALF: usize = FRAME / 2;
/// Samples per subframe, the unit the two codebooks are searched over.
pub const SUBFRAME: usize = FRAME / 4;

pub mod analysis;
pub mod bands;
pub mod bitstream;
pub mod codebook;
pub mod convolve;
pub mod decoder;
pub mod encoder;
pub mod excitation;
pub mod fixed;
pub mod frontend;
pub mod gain;
pub mod lpc;
pub mod lsf_weight;
pub mod lsp;
pub mod ltp;
pub mod mode;
pub mod noise;
pub mod pitch;
pub mod pitch_search;
pub mod postfilter;
pub mod preprocess;
pub mod pulses;
pub mod shape_pair;
pub mod shaping;
pub mod synth;
pub mod tables;
pub mod ulaw;
pub mod vad;
pub mod weight_lpc;
pub mod weighting;
pub mod weights;

pub use bitstream::FRAME_BYTES;
pub use decoder::Decoder;
pub use encoder::Encoder;
