//! Fixed-point arithmetic primitives.

/// Saturate an i32 value to i16.
#[inline]
pub fn sat16(x: i32) -> i16 {
    x.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

/// Q15 × Q15 → Q15, arithmetic shift without rounding (matches the DSP `mpy + sftr`).
/// Used from Phase 1 onward.
#[inline]
pub fn mul_q15_trunc(a: i16, b: i16) -> i16 {
    sat16(((a as i32) * (b as i32)) >> 15)
}

/// Q15 × Q15 → Q15 with nearest-neighbor rounding. Used from Phase 1 onward.
#[inline]
pub fn mul_q15_round(a: i16, b: i16) -> i16 {
    let p = (a as i32) * (b as i32);
    sat16((p + (1 << 14)) >> 15)
}

/// i40 (= 40-bit signed integer in i64) → i32 saturation.
/// Equivalent to `sat32` in the gain pipeline; clamps the intermediate accumulator to the i32 range.
#[inline]
pub fn sat32(v: i64) -> i64 {
    v.clamp(i32::MIN as i64, i32::MAX as i64)
}

/// `acc32_16(v)` — extract a 17-bit signed value from the i40 accumulator.
///
/// Behavior: combines bits 16..32 (= h) with bit 32 (= LSB of g) to form a 17-bit
/// signed integer, sign-extending if bit 16 is set.
///
/// Mathematical meaning: when `v` is interpreted as the result of a Q-format
/// product, take its upper 17 bits (= 1 sign + 16 mantissa) as an i32.
#[inline]
pub fn acc32_16(v: i64) -> i32 {
    let g = ((v >> 32) & 0xff) as i32;
    let h = ((v >> 16) & 0xffff) as i32;
    let mut x = ((g & 1) << 16) | h;
    if (x & 0x10000) != 0 {
        x -= 0x20000;
    }
    x
}

/// `shift_acc40(v, n)` — arithmetic shift of a 40-bit signed value.
///
/// `n >= 0`: left shift (overflow wraps; saturate downstream if needed).
/// `n < 0`: sign-extending right shift after a 40-bit mask (= equivalent to `sxm=true`).
pub fn shift_acc40(src: i64, shift: i16) -> i64 {
    if shift >= 0 {
        src.wrapping_shl(shift as u32)
    } else {
        let rshift = -shift as u32;
        if rshift >= 40 {
            if src < 0 { -1 } else { 0 }
        } else {
            let mut out = ((src & 0x00ff_ffff_ffff) as u64 >> rshift) as i64;
            if src < 0 {
                out |= !((1i64 << (40 - rshift)) - 1);
            }
            out
        }
    }
}

/// `to_i40(v)` — truncate an i64 value to 40 bits and sign-extend from bit 39.
///
/// Reproduces i40 semantics after overflow: on `overflow`, the behavior is
/// "preserve the lower 40 bits and sign-extend".
#[inline]
pub fn to_i40(v: i64) -> i64 {
    let masked = v & 0x00ff_ffff_ffff;
    if (masked & 0x0080_0000_0000) != 0 {
        masked | !0x00ff_ffff_ffff
    } else {
        masked
    }
}

/// `exp_acc(v)` — return the exponent (= normalization shift amount) of an i40 value.
/// The returned `t` is the value such that `norm_acc_with_t(v, t)` maximizes the
/// absolute value after shifting (= concentrates the energy in the upper 16 bits).
/// Range: -8..7.
pub fn exp_acc(src: i64) -> i16 {
    if src == 0 {
        return 0;
    }
    let mut t = -8i16;
    let mut cnt = 38i32;
    if src < 0 {
        while cnt >= 0 && (src & (1i64 << cnt)) != 0 {
            cnt -= 1;
            t = t.wrapping_add(1);
        }
    } else {
        while cnt >= 0 && (src & (1i64 << cnt)) == 0 {
            cnt -= 1;
            t = t.wrapping_add(1);
        }
    }
    t
}

/// `norm_acc_with_t(v, t)` — normalize using the shift amount derived from `exp_acc` (with i40 truncation).
pub fn norm_acc_with_t(src: i64, t: i16) -> i64 {
    let mut ts = t & 0x003f;
    if (ts & 0x0020) != 0 {
        ts |= !0x003f;
    }
    to_i40(shift_acc40(src, ts))
}

/// `mpya_from_acc32_16(t, acc, frct_mul)` — compute `t * acc32_16(acc) * frct_mul`,
/// saturating to i32::MAX if the result exceeds it (= clamp only the positive
/// upper bound; the negative side is not clamped — asymmetric saturation).
///
/// Used in the postfilter input MAC.
#[inline]
pub fn mpya_from_acc32_16(t: i16, acc: i64, frct_mul: i64) -> i64 {
    let p = i64::from(t) * i64::from(acc32_16(acc)) * frct_mul;
    if p > 0x7fff_ffff { 0x7fff_ffff } else { p }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sat16_clamps() {
        assert_eq!(sat16(i32::MAX), i16::MAX);
        assert_eq!(sat16(i32::MIN), i16::MIN);
        assert_eq!(sat16(0), 0);
        assert_eq!(sat16(12345), 12345);
    }

    #[test]
    fn mul_q15_trunc_one_times_x() {
        // Q15 1.0 (= 32767, not 32768) × x ≈ x; allow a tiny error since there is no rounding.
        let half = i16::MAX; // ~ +0.99997
        assert!((mul_q15_trunc(half, 16384) - 16383).abs() <= 1);
    }

    #[test]
    fn acc32_16_known_values() {
        // v = 0x07FE0000 (positive small): upper 17 bits = 0x07FE → 2046
        assert_eq!(acc32_16(0x07FE0000), 2046);
        // v with bit 32 set: result negative-extended
        assert_eq!(acc32_16(0x100000000), -65536);
    }

    #[test]
    fn shift_acc40_left_and_right() {
        assert_eq!(shift_acc40(1, 8), 256);
        assert_eq!(shift_acc40(256, -4), 16);
        // Negative input arithmetic right shift
        assert_eq!(shift_acc40(-256, -4), -16);
    }

    #[test]
    fn to_i40_truncates_and_sign_extends() {
        // Within the 40-bit range: passes through unchanged.
        assert_eq!(to_i40(123), 123);
        // When bit 39 is set, sign-extend.
        assert_eq!(to_i40(0x0080_0000_0000), -0x0080_0000_0000);
    }
}
