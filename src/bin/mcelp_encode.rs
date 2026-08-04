//! Encode raw mu-law on stdin to an M-CELP bit stream on stdout.
//!
//! Input is 320 mu-law bytes per frame; output is one hex-encoded frame per
//! line, which is what `mcelp_decode` reads back.  A trailing partial frame is
//! dropped, as the reference `mcelp_codec` drops it.

use std::io::{Read, Write};

use mcelp::{Encoder, FRAME, bitstream};

fn main() -> std::io::Result<()> {
    let mut pcm = Vec::new();
    std::io::stdin().lock().read_to_end(&mut pcm)?;

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    let mut encoder = Encoder::new();

    for block in pcm.chunks_exact(FRAME) {
        let frame = encoder.encode(block.try_into().unwrap());
        out.write_all(bitstream::to_hex_line(&frame).as_bytes())?;
        out.write_all(b"\n")?;
    }
    out.flush()
}
