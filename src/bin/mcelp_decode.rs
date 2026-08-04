//! Decode an M-CELP bit stream on stdin to raw mu-law on stdout.
//!
//! Input is one hex-encoded 18-byte frame per line; output is 320 mu-law bytes
//! per frame.  This mirrors the reference `mcelp_codec` in its decode mode.

use std::io::{BufRead, Write};

use mcelp::{Decoder, bitstream};

fn main() -> std::io::Result<()> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    let mut decoder = Decoder::new();

    for line in stdin.lock().lines() {
        let Some(frame) = bitstream::parse_hex_line(&line?) else {
            continue;
        };
        if let Some(pcm) = decoder.decode(&frame) {
            out.write_all(&pcm)?;
        }
    }
    out.flush()
}
