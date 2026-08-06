# How the implementation is verified

The codec is specified by its arithmetic: a decoder that rounds differently
produces different output from the same bits. So "it sounds right" is not a
useful test. The suite checks that this implementation reproduces the
reference **bit for bit**, and it does so at three levels.

```sh
cargo test --release
```

125 test files, 117 data files, 4.1 MB of captured reference values.

## Level 1 — stage tests

Most of the suite is one file per stage, replaying values captured from the
reference through a single function.

```
tests/levinson.rs       ← src/lpc.rs        levinson()
tests/chebyshev.rs      ← src/lsp.rs        the root search
tests/shortlist ...     ← src/pulses.rs     each block of the search
```

Each has a data file of the same name in `tests/data/`, one call per line:
inputs first, then the values the reference produced. The test's doc comment
names the columns.

This granularity is the point. Because encoder state flows forward, a wrong
value in one stage does not show up as one bad frame — it shows up as a
growing divergence several frames later. A stage test fails at the first call
where that stage's arithmetic differs, which localises the fault to a single
function.

## Level 2 — chain tests

A few tests compose two or three stages, checking that they hand values to
each other in the right form. These are named `*_chain`:

| test | chain |
|---|---|
| `quantise_chain` | cosine domain → frequencies → weights → quantiser → transport fields |
| `adaptive_gain_chain` | cross-correlation → adaptive gain |
| `gain_search_chain` | five measures → alignment → gain index |
| `lpc_chain`, `codebook_chain`, `shaping_chain` | the same idea elsewhere |

Chain tests catch a class of bug the stage tests cannot: two functions that are
each individually correct but disagree on a scale factor or an exponent
convention at their boundary.

## Level 3 — end to end

| test | what it pins |
|---|---|
| `encode_endtoend` | the encoder's output for both bundled examples, frame by frame |
| `endtoend` | the decoder's output for both, byte for byte, plus a concealment run |
| `api` | the byte-level round trip through the public API |

`examples/` carries a matched set — `*.orig.ulaw` (input), `*.mcelp` (the
reference encoder's bit stream), `*.decoded.ulaw` (the reference decoder's
output) — so both directions are anchored to the same material. Any change
that moves a single bit anywhere fails these.

## What the tests do not cover

Worth being explicit about:

- **Only two examples.** 75 frames each of one male and one female speaker.
  The paths those two recordings happen to exercise are well covered; a path
  neither of them takes is covered only by whichever stage test reaches it.
- **Data-indexed tables.** A codebook entry that neither example selects is
  never read. The tables are complete, but their far corners are untested.
- **Concealment is self-consistent, not reference-checked.** The examples
  contain no erased frames, so `endtoend` synthesises some by setting the
  suppression bit. That checks the concealment path runs and stays in step
  with itself; it does not check it against the reference's concealment.

## Working on the code

A few things follow from the above.

**Run the whole suite, in release.** Debug builds are perhaps twenty times
slower and the suite replays a lot of data. Both profiles are worth running
before committing — debug has overflow checks and `debug_assert`s that release
does not.

**Do not "fix" the arithmetic.** `acc` marks where a real machine wraps; it is
part of the specification, not a bug being papered over. Replacing it with
saturation, or reordering an accumulation, changes the output. See
[fixed-point.md](fixed-point.md).

**Prefer a stage test when adding coverage.** If you find a discrepancy,
capturing the inputs and outputs of the smallest function that shows it is
worth more than another end-to-end case.

**Watch for table boundaries.** Several tables share storage or are read at an
offset from a neighbour. Bounds are checked now, so an out-of-range read
panics rather than silently returning the neighbouring table's data — but if
you resize a table, check the code that indexes it rather than assuming its
nominal length.
