# Fixed-point arithmetic

Every value in the codec is an integer. There is no floating point anywhere,
and this is not an implementation detail that could be changed: the codec is
specified by its arithmetic. A decoder that rounds differently produces a
different bit stream from the same input, and the two implementations diverge.

`src/fixed.rs` holds the primitives; the rest of the code is written in terms
of them.

## Q notation

A Q*n* value is an integer to be read as a fraction scaled by 2ⁿ.

| format | range | used for |
|---|---|---|
| Q15 | −1.0 … 0.99997 | signals, gains, most coefficients |
| Q13 | −4.0 … 3.9999 | line spectral frequencies (π = 25736) |
| Q12 | −8.0 … 7.9998 | direct-form predictor coefficients |

Q15 cannot represent 1.0. Where a coefficient needs to be one, the code uses
16384 in Q14, or handles the case separately — this is why a few tables carry
values that look like half of what they should be.

## The three properties

The arithmetic reproduces three properties of a 16-bit fixed-point machine.
They are unusual enough that code written against ordinary integer semantics
will not match.

### 1. Fractional multiplies

Every product is shifted left by one. Q15 × Q15 would naturally land in Q30;
the extra shift puts it in Q31, so that the high half of the accumulator is
again a Q15 value that can be stored directly.

```rust
pub fn mul(x: i16, y: i16) -> i64 {
    (x as i64) * (y as i64) * 2
}
```

The doubling is why so many table values are halved: a coefficient that should
be 0.8 is stored as 0.4 in Q15, because the multiply will double it back.

### 2. Wide accumulators that wrap

The accumulator is 40 bits, and it **wraps** rather than saturating.
Saturation happens only where the code explicitly asks for it.

```rust
pub fn acc(v: i64) -> i64 { (v << 24) >> 24 }   // fold into 40 bits
pub fn sat(v: i64) -> i64 { v.clamp(i32::MIN as i64, i32::MAX as i64) }
```

This is the property most likely to trip up a reader. An intermediate sum that
overflows 32 bits is not clipped — it wraps, and the wrapped value is used. In
several places the algorithm *depends* on that: an energy accumulated over a
loud subframe can pass through an overflowed state and come back. Replacing
`acc` with saturation "to be safe" changes the output.

The eight extra bits above 32 are headroom: they let a long accumulation of
32-bit products run without any intermediate clipping, which is exactly what
the correlation and energy loops need.

### 3. Explicit rounding and truncation

Rounding is not the language's. Where the codec rounds, it adds a half at
bit 15 and discards what is below:

```rust
pub fn round(v: i64) -> i64 { acc(v + 0x8000) & !0xffff }
```

And where a value passes through 16-bit storage, the code says so — `trunc32`
keeps exactly the 32 bits a store-and-load round-trips, discarding accumulator
bits above them.

## Division

There is no divide instruction to model, so division is built from repeated
conditional subtract-and-shift:

- `restoring_divide` — the primitive, one quotient bit per step;
- `normalised_reciprocal` — normalise, then divide 16383 by the magnitude,
  returning quotient and exponent;
- `divide` — a signed normalised divide keeping the scale of a Q15 ratio;
- `sqrt`, `inverse_sqrt` — square root by binary search, and its reciprocal.

Most of the codec avoids division entirely. Every "which candidate is better"
comparison is a ratio, and ratios are compared by cross-multiplication:
`a/b > c/d` becomes `a·d > c·b`. The searches never divide.

## Normalisation

With no exponent field, dynamic range is managed by hand. The pattern
throughout is:

1. count the redundant sign bits of a value (`exp`);
2. shift it up by that much, so it uses the full word;
3. do the arithmetic at full precision;
4. carry the shift alongside the result and apply it at the end.

Several functions return a mantissa and an exponent as a pair for exactly this
reason, and the gain search's first job is to bring five such pairs onto a
common exponent before it can add them.

## Reading the code

Two habits make the arithmetic readable:

- **An `i64` is an accumulator, not a number.** Its low 16 bits are the "low
  half", bits 16–31 the "high half" that gets stored, and bits 32–39 headroom.
  `hi` and `low` name the two halves.
- **`acc` marks where a real machine would have wrapped.** It is not
  defensive; it is part of the specification. If you find yourself wanting to
  remove one, the value it is folding is one the algorithm expects to see
  folded.
