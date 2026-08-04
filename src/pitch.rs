//! Adaptive codebook: pitch lag decoding and fractional-delay prediction.
//!
//! Lags are carried at 1/3-sample resolution over the range 20..=143.  The
//! first subframe of each half-frame sends an absolute lag in 8 bits, the
//! second sends a 5-bit delta relative to it.

use crate::fixed::{acc, hi, low, mul, restoring_divide, sat, shift, trunc32};
use crate::tables::{LAG_INTERP, LAG_INTERP_BACKWARD as BACKWARD, LAG_INTERP_FORWARD as FORWARD};

/// Samples per subframe.
pub use crate::SUBFRAME;
/// Taps on each side of the fractional interpolation filter.
const INTERP_TAPS: usize = 10;
/// Sub-sample resolution of the pitch lag.
const RESOLUTION: i16 = 3;

/// Integer lag plus a fractional part in `-1..=1` thirds of a sample.
#[derive(Clone, Copy, Debug, Default)]
pub struct Lag {
    /// Integer part, 20..=143.
    pub integer: i16,
    /// Fractional part in thirds of a sample.
    pub frac: i16,
}

/// Decode the 8-bit absolute lag sent by the first subframe of a half-frame.
pub fn decode_absolute(code: i16) -> Lag {
    if code >= 197 {
        // The top of the range is coded at integer resolution only.
        return Lag {
            integer: code - 112,
            frac: 0,
        };
    }
    let quotient = low(restoring_divide((code as i64) + 2, RESOLUTION, 16));
    let integer = quotient + 19;
    let frac = code - low(mul(integer, RESOLUTION) >> 1) + 58;
    Lag { integer, frac }
}

/// Decode the 5-bit lag delta sent by the second subframe of a half-frame.
pub fn decode_relative(code: i16, prev: i16) -> Lag {
    let mut min = (prev - 5).max(20);
    if min + 9 > 143 {
        min = 134;
    }
    let quotient = low(restoring_divide((code as i64) + 2, RESOLUTION, 16)) - 1;
    Lag {
        integer: quotient + min,
        frac: code - 2 - low(mul(quotient, RESOLUTION) >> 1),
    }
}

/// Build the adaptive-codebook contribution by resampling the excitation
/// history at a fractional delay.
///
/// `exc[..at]` holds past excitation; `exc[at..at + SUBFRAME]` is written.
pub fn predict(exc: &mut [i16], at: usize, lag: &Lag) {
    // The filter is indexed by the *negated* fraction; a negative fraction is
    // folded into the integer delay.
    let mut phase = -lag.frac;
    let mut delay = lag.integer as usize;
    if phase < 0 {
        phase += RESOLUTION;
        delay += 1;
    }
    let phase = phase as usize;
    let forward = &crate::tables::PITCH_INTERP[phase..phase + 3 * INTERP_TAPS - 2];
    let mirrored = RESOLUTION as usize - phase;
    let backward = &crate::tables::PITCH_INTERP[mirrored..mirrored + 3 * INTERP_TAPS - 2];

    for n in 0..SUBFRAME {
        let base = at + n - delay;
        let mut a = 0i64;
        for k in 0..INTERP_TAPS {
            a = acc(a + mul(exc[base - k], forward[3 * k]));
        }
        for k in 0..INTERP_TAPS {
            a = acc(a + mul(exc[base + 1 + k], backward[3 * k]));
        }
        exc[at + n] = hi(sat(acc(a + 0x8000)));
    }
}

/// Periodic enhancement of the fixed-codebook vector ("pitch sharpening").
///
/// Only applied when the lag is shorter than a subframe: the
/// innovation is fed through the same fractional-delay predictor as the
/// adaptive codebook and mixed back in with the previous subframe's pitch gain.
pub fn sharpen(code: &mut [i16; SUBFRAME], lag: &Lag, gain: i16) {
    // The first output is produced ten samples before the delay itself, because
    // the interpolation filter straddles the delay point.
    let mut first = (lag.integer - INTERP_TAPS as i16) as i32;
    let mut phase = -lag.frac;
    if phase < 0 {
        phase += RESOLUTION;
        first += 1;
    }
    if first >= SUBFRAME as i32 {
        return;
    }
    let phase = phase as usize;
    let first = first as usize;

    // Guard region standing in for the (zero) innovation before this subframe,
    // followed by the innovation itself as it is progressively rewritten.
    const GUARD: usize = 2 * INTERP_TAPS + 1;
    let mut buf = [0i16; GUARD + SUBFRAME];
    buf[GUARD..GUARD + first].copy_from_slice(&code[..first]);

    let forward = &crate::tables::PITCH_INTERP[phase..phase + 3 * INTERP_TAPS - 2];
    let mirrored = RESOLUTION as usize - phase;
    let backward = &crate::tables::PITCH_INTERP[mirrored..mirrored + 3 * INTERP_TAPS - 2];

    for n in first..SUBFRAME {
        let centre = GUARD + (n - first) - INTERP_TAPS;
        let mut a = 0i64;
        for k in 0..INTERP_TAPS {
            a = acc(a + mul(buf[centre - k], forward[3 * k]));
        }
        for k in 0..INTERP_TAPS {
            a = acc(a + mul(buf[centre + 1 + k], backward[3 * k]));
        }
        let scaled = crate::fixed::mul32x16(acc(a), gain);
        let mixed = acc(crate::fixed::shift(scaled, 1) + ((code[n] as i64) << 16) + (1 << 15));
        code[n] = hi(mixed);
        buf[GUARD + n] = code[n];
    }
}

/// Code at which the absolute range stops carrying a fractional part.
const COARSE_FROM: i16 = 197;
/// Offset the coarse codes start at.
const COARSE_BIAS: i16 = 112;

/// Code the first subframe's lag absolutely (the inverse of
/// [`decode_absolute`]).
pub fn encode_absolute(lag: &Lag) -> i16 {
    // The switch is on the code, not the lag: an integer that would just reach
    // the coarse range still codes finely if its fraction pulls it back below.
    let fine = RESOLUTION * lag.integer - 58 + lag.frac;
    if fine >= COARSE_FROM {
        lag.integer + COARSE_BIAS
    } else {
        fine
    }
}

/// Code the second subframe's lag against the first (the inverse of
/// [`decode_relative`]).
pub fn encode_relative(lag: &Lag, previous: i16) -> i16 {
    let mut min = (previous - 5).max(20);
    if min + 9 > 143 {
        min = 134;
    }
    RESOLUTION * (lag.integer - min) + 2 + lag.frac
}

/// Fractional-delay interpolation tables, read with a stride of three: one for
/// the taps ahead of the sample and one for those behind it.
const TAPS: usize = 8;
/// Zeros the interpolation starts against, so that the first few outputs see a
/// clean history.
const PAD: usize = 9;

/// Extend the adaptive codebook vector periodically.
///
/// When the lag is shorter than a subframe the vector runs off the end of the
/// history, so the part that is missing is filled in by repeating what is there
/// one period earlier, resampled at the fractional part of the lag by an
/// eight-tap interpolator.  The result is scaled by `gain` and added into what
/// the buffer already holds.
///
/// The extension feeds on itself: once more than one period has been produced,
/// the interpolator's window reaches into samples this same loop wrote.
pub fn extend(buffer: &mut [i16; SUBFRAME], lag: i16, fraction: i16, gain: i16) {
    // Only the `sub #4,a` sits in the branch's delay slot; the other two
    // adjustments belong to the fall-through alone.
    let backwards = acc(-(fraction as i64));
    let (period, phase) = if backwards >= 0 {
        (acc((lag as i64) - 4), backwards)
    } else {
        (acc((lag as i64) - 3), acc(backwards + 3))
    };
    if acc(period - SUBFRAME as i64) >= 0 {
        return;
    }
    let (period, phase) = (period as usize, phase as usize);

    let mut weight = [0i16; TAPS];
    for k in 0..TAPS / 2 {
        weight[k] = LAG_INTERP[FORWARD + phase + 3 * k];
        weight[TAPS / 2 + k] = LAG_INTERP[BACKWARD - phase + 3 * k];
    }

    // Nine zeros, one period of history, then room for what is produced.
    let mut history = [0i16; PAD + SUBFRAME];
    history[PAD..PAD + period].copy_from_slice(&buffer[..period]);

    for k in 0..SUBFRAME - period {
        let centre = PAD - 4 + k;
        let mut total = 0i64;
        for j in 0..TAPS / 2 {
            total = acc(total + crate::fixed::mul(weight[j], history[centre - j]));
        }
        for j in 0..TAPS / 2 {
            total = acc(total + crate::fixed::mul(weight[TAPS / 2 + j], history[centre + 1 + j]));
        }
        let scaled = shift(crate::fixed::mul32x16(trunc32(total), gain), 1);
        let sum = acc(scaled + ((buffer[period + k] as i64) << 16) + (1 << 15));
        buffer[period + k] = hi(sum);
        history[PAD + period + k] = hi(sum);
    }
}
