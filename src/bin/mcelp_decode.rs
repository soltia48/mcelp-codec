//! Simple .mcelp decoder CLI.
//! Reads 36-char hex frames (one per line) from stdin and writes raw μ-law PCM to stdout.

use mcelp_codec::{Codec, MCELP_FRAME_BYTES, bitstream::parse_frame_hex};
use std::io::{self, BufRead, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = stdout.lock();

    let mut codec = Codec::new();
    let mut frame_count = 0usize;
    let mut reset_count = 0usize;

    for (line_no, line) in stdin.lock().lines().enumerate() {
        let line = match line {
            Ok(s) => s,
            Err(e) => {
                eprintln!("read stdin: {e}");
                return ExitCode::FAILURE;
            }
        };
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let frame: [u8; MCELP_FRAME_BYTES] = match parse_frame_hex(trimmed) {
            Some(f) => f,
            None => {
                eprintln!("line {}: invalid hex frame: {}", line_no + 1, trimmed);
                return ExitCode::FAILURE;
            }
        };
        frame_count += 1;
        match codec.decode_frame(&frame) {
            Some(pcm) => {
                if let Err(e) = out.write_all(&pcm) {
                    eprintln!("write stdout: {e}");
                    return ExitCode::FAILURE;
                }
            }
            None => reset_count += 1,
        }
    }

    if let Err(e) = out.flush() {
        eprintln!("flush stdout: {e}");
        return ExitCode::FAILURE;
    }
    eprintln!("decoded {frame_count} frames ({reset_count} reset)");
    ExitCode::SUCCESS
}
