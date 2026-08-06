# Analysis by synthesis — the encoder's subframe loop

The phrase "analysis by synthesis" hides a specific and slightly awkward piece
of machinery. This document walks through one subframe of it, in the order the
encoder actually performs the steps.

## The problem

We want to choose an excitation whose *synthesised output* is closest to the
input. The obvious formulation — synthesise every candidate and compare — is
far too expensive. Two observations make it tractable:

1. The synthesis filter is **linear**. If a candidate is a sum of scaled
   pulses, its output is the same sum of scaled pulse responses. So the filter
   can be applied once, to an impulse, and every candidate scored by combining
   precomputed vectors.
2. The filter's **memory** is the same for every candidate — it comes from the
   subframes already coded. That part of the output can be computed once and
   subtracted from the input, leaving a *target* that depends only on the
   excitation being chosen.

Both are standard; the second is what the target computation below is doing.

## One subframe, step by step

### 1. Residual

The subframe's window of the shaping signal goes through `A(z)` — the inverse
of the synthesis filter — giving the LPC residual.

```
residual = A(z) · shaping_signal
```

### 2. The residual goes into the excitation buffer

Before anything else, that residual is written into the excitation buffer at
the current subframe's position. This looks premature — the excitation for
this subframe has not been chosen yet — but it is deliberate: a pitch lag
shorter than 80 samples reaches into the part of the buffer this subframe has
not filled in. Seeding it with the residual is what the reference does, and the
adaptive codebook's periodic extension then overwrites it.

### 3. Re-synthesis with the error memory

```
resynthesised = 1/A(z) · residual        starting from error_memory
```

`error_memory` holds the last ten samples of the difference between the input
and what the decoder will reconstruct. Starting the filter from that memory
rather than from rest is the trick: what comes out is the input signal minus
the part of it the already-coded subframes will produce.

### 4. Weighting

```
weighted = A(z/γ₁) · [error_memory ‖ resynthesised]
target   = 1/A(z/γ₂) · weighted          starting from residual_memory
```

The result is the **target**: the weighted signal with everything the previous
subframes already account for taken out. This is what all three searches are
matched against.

### 5. Impulse response

The searches need the impulse response of the whole weighted synthesis chain:

```
δ → A(z/γ₁) → 1/A(z) → 1/A(z/γ₂) → h[0…79]
```

computed from rest, once per subframe, 80 samples long. Every candidate
excitation is convolved with `h` rather than filtered, which is what makes the
codebook searches affordable.

### 6. The three searches, in order

```
closed-loop pitch lag        → the adaptive contribution
    ↓
mode decision                → how this subframe is to be coded
    ↓
fixed-codebook search        → the innovation
    ↓
gain search                  → both gains at once
```

The order is forced by dependency. The adaptive contribution has to be known
before the fixed codebook can be searched against what remains; both have to
be known before the gains can be chosen; and the mode decision needs the
adaptive contribution to judge how periodic the subframe is.

Each search is scored on the same `correlation² / energy` ratio, compared by
cross-multiplication so that no division is needed anywhere in the loop.

### 7. Building the excitation and updating memories

Once both indices and the gain index are chosen, the excitation is assembled
exactly as the decoder will assemble it, and put back through:

- the synthesis filter, updating `speech_memory`;
- the weighting denominator, updating `residual_memory`;
- the input-minus-reconstruction difference, updating `error_memory`.

## Why the encoder mirrors the decoder

Notice that step 7 dequantises what step 6 just quantised, and runs the result
through the decoder's filters. The encoder never uses the unquantised values
it computed — it uses the values the decoder will see.

This is not defensive programming; it is required for correctness. The next
subframe's target is derived from these memories, so if the encoder tracked
its own ideal signal instead of the decoder's reconstruction, the two would
drift apart and every subsequent subframe would be optimised against a state
the decoder is not in.

The same discipline appears at the frame level. The line spectrum quantiser
feeds its *decoded* output back into the predictor memory, and the
interpolation decision is taken on the decoded filter's reflection
coefficients — so that the decoder, making the same test on the same values,
reaches the same conclusion without being told.

## Consequences for debugging

Because state flows forward, an encoder bug rarely produces one bad frame. It
produces a frame that is slightly wrong and then a growing divergence, because
every later subframe is being optimised against a wrong memory.

That is why the test suite checks each stage against reference data
independently rather than only end to end: a stage test fails at the first
frame where that stage's arithmetic differs, while the end-to-end test would
only tell you that something, somewhere, drifted. See
[testing.md](testing.md).
