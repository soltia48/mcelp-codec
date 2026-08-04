//! Development aid: dump decoder state for one frame so it can be diffed
//! against the instrumented reference.  Ignored by default.

use mcelp::{Decoder, bitstream};

#[test]
#[ignore]
fn dump_frame_state() {
    let dir = env!("CARGO_MANIFEST_DIR");
    let ex = std::env::var("EXAMPLE").unwrap_or_else(|_| "male".into());
    let target: usize = std::env::var("FRAME")
        .map(|s| s.parse().unwrap())
        .unwrap_or(0);
    let half: usize = std::env::var("HALF")
        .map(|s| s.parse().unwrap())
        .unwrap_or(0);
    let bits = std::fs::read_to_string(format!("{dir}/examples/{ex}.mcelp")).unwrap();
    let mut decoder = Decoder::new();
    for (n, line) in bits.lines().enumerate() {
        let frame = bitstream::parse_hex_line(line).unwrap();
        let p = bitstream::unpack(&bitstream::canonicalize(&frame));
        if n == target {
            let (a, b) = decoder.dump_spectrum(&p);
            println!("MID {a:?}");
            println!("CUR {b:?}");
            println!("STATE suppress={} {}", p.suppress, decoder.dump_state());
            decoder.dump_half(&p, half, if half == 0 { &a } else { &b });
            println!("STATE2 {}", decoder.dump_state());
            return;
        }
        decoder.decode_params(&p);
    }
}
