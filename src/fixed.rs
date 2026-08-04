//! Fixed-point primitives.
//!
//! The codec's arithmetic is that of a 16-bit fixed-point machine with
//! fractional multiplies, sign-extended memory operands and a wide accumulator
//! that wraps.  That gives three properties the code below reproduces:
//!
//! * multiplies are *fractional*: every product is shifted left by one so that
//!   Q15 x Q15 lands in Q31 rather than Q30 ([`mul`]);
//! * accumulators are 40 bit and **wrap** rather than saturate — saturation only
//!   happens where the code explicitly asks for it ([`sat`]);
//! * memory operands are sign extended into the accumulator.
//!
//! Everything is expressed as a plain `i64` holding the 40-bit accumulator
//! value; [`acc`] folds a result back into that range.

/// Fold a value into the 40-bit accumulator range (two's complement wrap).
#[inline]
pub fn acc(v: i64) -> i64 {
    (v << 24) >> 24
}

/// Saturate a 40-bit accumulator to 32 bits.
#[inline]
pub fn sat(v: i64) -> i64 {
    v.clamp(i32::MIN as i64, i32::MAX as i64)
}

/// High half of an accumulator.
#[inline]
pub fn hi(v: i64) -> i16 {
    (v >> 16) as i16
}

/// Low half of an accumulator.
#[inline]
pub fn low(v: i64) -> i16 {
    v as i16
}

/// Fractional multiply: Q15 x Q15 -> Q31, i.e. `x * y << 1`.
#[inline]
pub fn mul(x: i16, y: i16) -> i64 {
    (x as i64) * (y as i64) * 2
}

/// Fractional multiply of two unsigned 16-bit operands.
#[inline]
pub fn mul_uu(x: i16, y: i16) -> i64 {
    (x as u16 as i64) * (y as u16 as i64) * 2
}

/// Fractional multiply of an unsigned 16-bit operand by a signed one.
#[inline]
pub fn mul_us(x: i16, y: i16) -> i64 {
    (x as u16 as i64) * (y as i64) * 2
}

/// 32x16 fractional multiply, built from the operand's two halves.
/// `x` is a 32-bit accumulator value, `y` a Q15 coefficient.
#[inline]
pub fn mul32x16(x: i64, y: i16) -> i64 {
    let partial = acc(mul_us(low(x), y)) >> 16;
    acc(partial + mul(hi(x), y))
}

/// 32x16 fractional multiply with the low half halved before the multiply and
/// the partial product shifted back up, saturating at every step.
///
/// The reference spells the plain [`mul32x16`] this way wherever the two partial
/// products would otherwise collide in one accumulator.
#[inline]
pub fn scale32(x: i64, y: i16) -> i64 {
    let low_half = ((x as u32 as u16) >> 1) as i64;
    let partial = sat(shift(sat(shift(sat(low_half * (y as i64) * 2), -16)), 1));
    sat(acc(partial + (hi(x) as i64) * (y as i64) * 2))
}

/// Full 32x32 fractional product, keeping the upper half.
///
/// Accumulated a half at a time with the low halves read unsigned, which is how
/// the reference builds every 32-bit product that is not [`scale32`]'s shape.
/// Callers that need the result clipped wrap it in [`sat`].
pub fn mul32(x: i64, y: i64) -> i64 {
    let (xl, yl) = (low(x), low(y));
    let (xh, yh) = (hi(x), hi(y));
    let mut a = shift(acc(mul_uu(xl, yl)), -16);
    a = acc(a + mul_us(xl, yh));
    a = acc(a + mul_us(yl, xh));
    a = shift(a, -16);
    acc(a + mul(xh, yh))
}

/// Square of a 32-bit value, built out of its halves like any other product.
#[inline]
pub fn square32(v: i64) -> i64 {
    mul32(v, v)
}

/// Rounding: add a half at bit 15 and drop what is below it.
#[inline]
pub fn round(v: i64) -> i64 {
    acc(v + 0x8000) & !0xffff
}

/// A fractional multiply with the result rounded at bit 15.
#[inline]
pub fn mulr(x: i16, y: i16) -> i64 {
    round(sat(mul(x, y)))
}

/// Keep only the 32 bits a store/load pair round-trips through memory.
#[inline]
pub fn trunc32(v: i64) -> i64 {
    v as i32 as i64
}

/// Sign-extend the six-bit shift field a normalising shift is driven by.
#[inline]
pub fn norm_shift(amount: i32) -> i32 {
    let ts = amount & 0x3f;
    if ts & 0x20 != 0 { ts - 64 } else { ts }
}

/// Shift an accumulator by that field, saturating.
#[inline]
pub fn norm(v: i64, amount: i16) -> i64 {
    sat(shift(acc(v), norm_shift(amount as i32)))
}

/// Redundant-sign-bit count of a 40-bit accumulator: the left shift
/// that would normalise it against bit 31.  Zero for a zero accumulator.
pub fn exp(v: i64) -> i32 {
    let v = acc(v);
    if v == 0 {
        return 0;
    }
    let mag = if v < 0 { !v } else { v };
    let top = 63 - mag.leading_zeros() as i32; // index of the highest set bit
    30 - top
}

/// Arithmetic shift of an accumulator by a signed amount.
#[inline]
pub fn shift(v: i64, s: i32) -> i64 {
    if s >= 0 { acc(v << s) } else { acc(v >> (-s)) }
}

/// One step of a conditional subtract-and-shift restoring division.
#[inline]
fn restoring_step(num: i64, den: i16) -> i64 {
    let alu = num - ((den as i64) << 15);
    acc(if alu >= 0 { (alu << 1) + 1 } else { num << 1 })
}

/// Restoring division, one bit at a time, of `a / den`.
/// The quotient accumulates in the low bits of the returned accumulator.
pub fn restoring_divide(mut num: i64, den: i16, steps: u32) -> i64 {
    for _ in 0..steps {
        num = restoring_step(num, den);
    }
    num
}

/// Reciprocal-style divide used throughout: normalises `value`, then divides
/// 16383 by the normalised magnitude.  Returns `(quotient, shift)` where
/// `shift` is the exponent that has to be applied to the quotient afterwards.
pub fn normalised_reciprocal(value: i64) -> (i64, i32) {
    let e = exp(value);
    let shift_out = e + 1;
    let normalised = shift(value, e);
    let magnitude = sat(if normalised < 0 {
        -normalised
    } else {
        normalised
    }) >> 16;
    let quotient = restoring_divide((16383i64) << 16, magnitude as i16, 15);
    // The quotient sits in the low 16 bits; move it to the high half.
    let quotient = acc(((quotient as u64 as i64) & 0xffff) << 16);
    if hi(normalised) < 0 {
        (acc(-quotient), shift_out)
    } else {
        (quotient, shift_out)
    }
}

/// `numerator / value`, a signed normalised divide.
///
/// The quotient is left shifted by the divisor's exponent, so the result keeps
/// the scale of a Q15 ratio.
pub fn divide(value: i64, numerator: i16) -> i64 {
    let negative = value < 0;
    let magnitude = if negative { acc(-value) } else { value };
    let e = exp(magnitude);
    let normalised = shift(magnitude, e);
    let quotient = low(restoring_divide((16383i64) << 16, hi(normalised), 15));
    let product = acc((quotient as i64) * (numerator as i64) * 2);
    shift(if negative { acc(-product) } else { product }, e + 1)
}

/// Square root by binary search.
///
/// Fourteen refinements plus a final decrement, each testing the candidate by
/// squaring it; the candidate is held in sixteen bits, so everything wraps
/// there.  The result is returned in the high half.
pub fn sqrt(value: i64) -> i64 {
    let mut r: i16 = 16384;
    for step in 0..15 {
        let square = acc((r as i64) * (r as i64) * 2);
        let bit: i16 = if step == 14 { 0 } else { 0x2000 >> step };
        r = r.wrapping_add(bit);
        if acc(square - value) > 0 {
            r = r.wrapping_sub(if step == 14 { 1 } else { bit << 1 });
        }
    }
    (r as i64) << 16
}

/// Reciprocal of a square root.
///
/// Used to normalise a correlation by the energy of the segment it was taken
/// over.  A non-positive energy gives full scale, which makes the caller treat
/// the lag as perfectly correlated rather than dividing by nothing.
pub fn inverse_sqrt(value: i64) -> i64 {
    if acc(value) <= 0 {
        return 0x7fffffff;
    }
    let (quotient, shift_out) = normalised_reciprocal(sqrt(value));
    let amount = shift_out - 15;
    sat(shift(acc(quotient), amount))
}
