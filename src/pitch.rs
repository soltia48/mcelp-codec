//! Pitch lag decoding + adaptive codebook.

use crate::tables::pitch::INTERP_FILTER;

/// Pitch-lag state shared within a block (sub 0/1 or sub 2/3).
#[derive(Clone, Copy, Debug, Default)]
pub struct PitchLagState {
    pub prev_lag: i16,
}

/// Pitch-lag decoding result for a single subframe.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PitchLag {
    pub lag: i16,
    pub sub_lag: i16,
}

pub fn decode_lag_absolute(phase_input: i16) -> PitchLag {
    if phase_input < 197 {
        let q = ((phase_input as i32) + 2) / 3;
        let lag = (q + 19) as i16;
        let sub_lag = (phase_input as i32) - (lag as i32) * 3 + 58;
        PitchLag {
            lag,
            sub_lag: sub_lag as i16,
        }
    } else {
        PitchLag {
            lag: phase_input - 112,
            sub_lag: 0,
        }
    }
}

pub fn decode_lag_differential(phase_input: i16, prev_lag: i16) -> PitchLag {
    let mut anchor = (prev_lag as i32 - 5).max(20);
    if anchor + 9 - 143 > 0 {
        anchor = 134;
    }
    let q = ((phase_input as i32) + 2) / 3;
    let t1 = q - 1;
    let lag = (anchor + t1) as i16;
    let sub_lag = (phase_input as i32) - 2 - t1 * 3;
    PitchLag {
        lag,
        sub_lag: sub_lag as i16,
    }
}

pub fn decode_subframe_lag(state: &mut PitchLagState, phase: i16, phase_input: i16) -> PitchLag {
    let result = if phase == 0 {
        decode_lag_absolute(phase_input)
    } else {
        decode_lag_differential(phase_input, state.prev_lag)
    };
    state.prev_lag = result.lag;
    result
}

pub fn decode_lag_fract(sub_lag: i16) -> (usize, i16) {
    let mut a = -(sub_lag as i32);
    let mut lag_adjust = 0i16;
    if a < 0 {
        a += 3;
        lag_adjust = 1;
    }
    (a as usize, lag_adjust)
}

/// Generate 80 samples of fractional-pitch interpolation output.
/// `past_excitation` is written back via self-feedback (= within the same
/// call it is referenced as the future tap for subsequent samples).
pub fn pitch_adaptive_codebook(
    past_excitation: &mut [i16; 320],
    base_offset: usize,
    write_offset: usize,
    fract: usize,
    output: &mut [i16; 80],
) {
    debug_assert!(fract < 3);
    let h = INTERP_FILTER[fract];
    let h_rev: [i16; 10] = match fract {
        0 => {
            let h0 = INTERP_FILTER[0];
            [
                h0[1], h0[2], h0[3], h0[4], h0[5], h0[6], h0[7], h0[8], h0[9], 0,
            ]
        }
        1 => INTERP_FILTER[2],
        2 => INTERP_FILTER[1],
        _ => unreachable!(),
    };

    for n in 0..80 {
        let mut acc: i64 = 0;
        for k in 0..10 {
            acc += (past_excitation[base_offset + n - k] as i64) * (h[k] as i64);
            acc += (past_excitation[base_offset + n + 1 + k] as i64) * (h_rev[k] as i64);
        }
        // Q15 × Q15 = Q30; round to Q15.
        let q15 = ((acc + (1 << 14)) >> 15).clamp(i16::MIN as i64, i16::MAX as i64) as i16;
        output[n] = q15;
        past_excitation[write_offset + n] = q15;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_lag_absolute_known_values() {
        let r = decode_lag_absolute(98);
        assert_eq!(r.lag, 52);
        assert_eq!(r.sub_lag, 0);
        let r = decode_lag_absolute(103);
        assert_eq!(r.lag, 54);
        assert_eq!(r.sub_lag, -1);
    }

    #[test]
    fn decode_lag_differential_known_values() {
        let r = decode_lag_differential(6, 52);
        assert_eq!(r.lag, 48);
        assert_eq!(r.sub_lag, 1);
    }

    #[test]
    fn decode_lag_differential_clamps() {
        // Lower bound: prev_lag=10 → anchor=20, phase=0 → t1=-1 → lag=19
        assert_eq!(decode_lag_differential(0, 10).lag, 19);
        // Upper bound: prev_lag=200 → anchor=134, phase=31 → t1=10 → lag=144
        assert_eq!(decode_lag_differential(31, 200).lag, 144);
    }

    #[test]
    fn decode_lag_fract_table() {
        assert_eq!(decode_lag_fract(-1), (1, 0));
        assert_eq!(decode_lag_fract(0), (0, 0));
        assert_eq!(decode_lag_fract(1), (2, 1));
    }

    #[test]
    fn pitch_adaptive_codebook_zero_input_zero_output() {
        // All-zero excitation → all-zero output.
        let mut buf = [0i16; 320];
        let mut out = [0i16; 80];
        pitch_adaptive_codebook(&mut buf, 100, 200, 0, &mut out);
        assert!(out.iter().all(|&v| v == 0));
    }
}
