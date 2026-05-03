# mcelp-codec

A Rust implementation of the **Mitsubishi CELP speech codec** (三菱CELP方式音声コーデック),
produced by reverse-engineering an existing implementation.

The current crate ships a **decoder only**. An encoder may be added in the future.

## Frame format

- **Input** — `.mcelp` text file, one frame per line, each line a 36-character hex
  string encoding an 18-byte MCELP frame.
- **Output** — raw 8 kHz mono μ-law PCM. Each 18-byte frame decodes to 320 bytes
  (40 ms) of μ-law samples.

## Build

```sh
cargo build --release
```

## Usage

The `mcelp_decode` binary reads hex frames from stdin and writes raw μ-law PCM to stdout.
Decode one of the bundled examples:

```sh
cargo run --release --bin mcelp_decode < examples/female.mcelp > female.ulaw
```

Or pipe directly into a player. With `sox`:

```sh
cargo run --release --bin mcelp_decode < examples/female.mcelp | play -t ul -r 8000 -c 1 -
```

…or with `ffplay`:

```sh
cargo run --release --bin mcelp_decode < examples/female.mcelp | ffplay -f mulaw -ar 8000 -ac 1 -
```

## Library

```rust
use mcelp_codec::{Codec, MCELP_FRAME_BYTES};

let mut codec = Codec::new();
let frame: [u8; MCELP_FRAME_BYTES] = /* 18 bytes */;
if let Some(pcm) = codec.decode_frame(&frame) {
    // pcm: [u8; 320] of μ-law
}
```

## License

MIT — see [LICENSE](LICENSE).
