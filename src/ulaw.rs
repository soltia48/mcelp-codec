//! μ-law companding.

use crate::tables::ulaw::ULAW_TO_LIN;

/// μ-law byte → linear PCM (Q15). LUT lookup.
#[inline]
pub fn ulaw_to_linear_i16(byte: u8) -> i16 {
    ULAW_TO_LIN[byte as usize]
}

/// linear PCM (Q15) → μ-law byte.
pub fn linear_i16_to_ulaw(sample: i16) -> u8 {
    // §C.1: one's-complement absolute value, >>2, bias 33, 13-bit saturation.
    let m0: i32 = if sample >= 0 {
        sample as i32
    } else {
        !(sample as i32)
    };
    let m: i32 = ((m0 >> 2) + 33).min(8191);

    // §C.2: segment decomposition.
    let seg_bits = m >> 6;
    let seg_shift: i32 = if seg_bits == 0 {
        0
    } else {
        seg_bits.ilog2() as i32 + 1
    };
    let mantissa = (m >> (seg_shift + 1)) & 0x0F;
    let segment = 7 - seg_shift;

    // §C.3: assembly.
    let mut ulaw = ((segment << 4) | (15 - mantissa)) as u8;
    if sample >= 0 {
        ulaw |= 0x80;
    }
    ulaw
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ulaw_to_linear_known_values() {
        assert_eq!(ulaw_to_linear_i16(0) as u16, 0x8284);
        assert_eq!(ulaw_to_linear_i16(1) as u16, 0x8684);
        assert_eq!(ulaw_to_linear_i16(255), 0);
    }

    #[test]
    fn linear_to_ulaw_known_values() {
        assert_eq!(linear_i16_to_ulaw(0), 0xff);
        assert_eq!(linear_i16_to_ulaw(-1), 0x7f);
        // saturation at extremes must not panic
        let _ = linear_i16_to_ulaw(i16::MAX);
        let _ = linear_i16_to_ulaw(i16::MIN);
    }

    /// All 256 bytes roundtrip without panicking (= linear_i16_to_ulaw is safe across all of i16).
    /// **A complete roundtrip is not guaranteed**: bytes 127 and 255 both decode
    /// to linear 0, and `linear_i16_to_ulaw(0) == 0xFF` collapses them to
    /// positive zero (G.711 negative-zero folding). This is a property of the
    /// table itself.
    #[test]
    fn roundtrip_no_panic() {
        for b in 0u8..=255 {
            let lin = ulaw_to_linear_i16(b);
            let _ = linear_i16_to_ulaw(lin);
        }
    }

    /// Verify the roundtrip matches except for linear 0 (= negative-zero folding only affects bytes 127/255).
    #[test]
    fn roundtrip_nonzero_linear() {
        for b in 0u8..=255 {
            let lin = ulaw_to_linear_i16(b);
            if lin == 0 {
                continue;
            }
            let back = linear_i16_to_ulaw(lin);
            assert_eq!(back, b, "byte {b}: lin={lin}");
        }
    }
}
