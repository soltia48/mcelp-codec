# mcelp-codec

[![CI](https://github.com/soltia48/mcelp-codec/actions/workflows/ci.yml/badge.svg)](https://github.com/soltia48/mcelp-codec/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/mcelp.svg)](https://crates.io/crates/mcelp)
[![docs.rs](https://docs.rs/mcelp/badge.svg)](https://docs.rs/mcelp)

A Rust implementation of the **Mitsubishi CELP speech codec** (三菱CELP方式音声コーデック).

The crate ships both halves of the codec. Each stage is written as the
algorithm it implements, and both the encoder and the decoder reproduce the
reference bit for bit — the regression suite replays the reference's own
output, stage by stage.

## Frame format

- **Audio** — raw 8 kHz mono μ-law PCM, 320 bytes (40 ms) per frame.
- **Bit stream** — 18 bytes per frame, i.e. 3.6 kbit/s.
- **Container** — `.mcelp` text file, one frame per line, each line a
  36-character hex string encoding one 18-byte frame.

A trailing partial frame is dropped, as the reference implementation drops it.

## Library

```toml
[dependencies]
mcelp = "1.0"
```

`Encoder` and `Decoder` are the entry points. Both carry state between frames,
so one instance has to see a stream in order.

```rust
use mcelp::{Decoder, Encoder, FRAME};

let speech: Vec<u8> = read_ulaw();          // 8 kHz mono μ-law

let (mut encoder, mut decoder) = (Encoder::new(), Decoder::new());
let mut out = Vec::new();
for block in speech.chunks_exact(FRAME) {
    let frame = encoder.encode(block.try_into().unwrap());   // [u8; 18]
    if let Some(pcm) = decoder.decode(&frame) {              // [u8; 320]
        out.extend_from_slice(&pcm);
    }
}
```

`Decoder::decode` returns `None` only for an in-band reset frame, which resets
the decoder instead of producing audio; an erased frame is concealed and still
yields samples. `Decoder::decode_linear` gives the same samples as 16-bit
linear PCM instead of μ-law.

`mcelp::bitstream` exposes the frame's fourteen parameter fields
(`unpack`/`pack`) and the hex container (`parse_hex_line`/`to_hex_line`) for
callers that need to work below the frame.

The remaining modules are the stages the two entry points are built from. They
are public so that each can be replayed against the reference on its own, which
is what the test suite does; they are not a stable interface.

## Command line

```sh
cargo build --release
```

`mcelp_encode` reads raw μ-law on stdin and writes hex frames on stdout;
`mcelp_decode` does the reverse.

```sh
cargo run --release --bin mcelp_encode < speech.ulaw > speech.mcelp
cargo run --release --bin mcelp_decode < examples/female.mcelp > female.ulaw
```

Or pipe straight into a player — with `sox`:

```sh
cargo run --release --bin mcelp_decode < examples/female.mcelp | play -t ul -r 8000 -c 1 -
```

…or with `ffplay`:

```sh
cargo run --release --bin mcelp_decode < examples/female.mcelp | ffplay -f mulaw -ar 8000 -ac 1 -
```

## Tests

```sh
cargo test --release
```

Each stage is checked against data captured from the reference, and
`tests/api.rs` pins the byte-level round trip. `tests/endtoend.rs` and
`tests/encode_endtoend.rs` carry the whole of both bundled examples, so a
change that moves a single bit fails the suite.

## Continuous integration

Every push and pull request runs, in `.github/workflows/ci.yml`:

| job | what it checks |
|---|---|
| Lint | `cargo fmt --check`, `cargo clippy -D warnings`, `cargo doc -D warnings` |
| Test | the suite on Linux, macOS and Windows in release, plus Linux in debug for the overflow checks and `debug_assert`s release drops |
| MSRV | the whole suite on the `rust-version` declared in `Cargo.toml` |
| Package | `cargo package`, guarding the size so the excluded regression data cannot creep back in |

Each test job also runs the two command line tools over the bundled examples
and compares the bytes, so the codec is shown to be bit exact on every
platform, not merely to compile there.

Tagging a commit `v<version>` runs `.github/workflows/release.yml`, which
checks the tag against `Cargo.toml`, publishes to crates.io and attaches
prebuilt tools to a GitHub release.

## Documentation

`docs/` holds a technical description of the codec, stage by stage — frame
structure, the two pipelines, each coding stage, and the fixed-point model
everything runs on.

- [English](docs/en/README.md)
- [日本語](docs/ja/README.md)

## License

MIT — see [LICENSE](LICENSE).
