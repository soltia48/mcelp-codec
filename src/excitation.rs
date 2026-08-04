//! Building a subframe's excitation once both codebooks have been searched.
//!
//! The closed loop leaves an adaptive gain, an innovation and a fixed gain
//! behind.  What is left is bookkeeping: clamping the adaptive gain to what a
//! periodic frame is allowed and to what the next half-frame may carry over,
//! deciding the two gains the coming periodic extension will use, summing the
//! two contributions into the excitation, and putting the result back through
//! the analysis filters so the next subframe starts from the right memories.

use crate::fixed::{acc, hi, low, mul, sat, shift};

/// Largest adaptive gain a periodic half-frame is allowed.
const PERIODIC_CAP: i64 = 15565;

/// Cap the adaptive gain when the half-frame looks periodic.
///
/// A strongly periodic half-frame would otherwise let the adaptive codebook
/// take almost all of the target, leaving the fixed codebook nothing to correct
/// with; holding the gain below about a half keeps some of the residual for it.
pub fn cap_gain(periodic: i16, gain: &mut i16) {
    if periodic == 0 {
        return;
    }
    *gain = low(acc((*gain as i64).min(PERIODIC_CAP)));
}

/// Range the adaptive gain is held to when it is carried to the next
/// half-frame.
const CARRY_LOW: i64 = 3277;
const CARRY_HIGH: i64 = 13017;

/// Hold the adaptive gain inside a working range before carrying it forward.
///
/// The gain the next half-frame's periodic extension uses is this one, clamped:
/// too small and the extension contributes nothing, too large and it rings.
pub fn carry_gain(gain: i16) -> i16 {
    low(acc((gain as i64).clamp(CARRY_LOW, CARRY_HIGH)))
}

/// A half, the value both extension gains are set to when they are set at all.
const HALF: i64 = 16384;
/// The previous gain has to be above this for the second gain to be allowed.
const WAS_STRONG: i64 = 13107;
/// How far the lag may have moved and still count as the same pitch.
const SAME_PITCH: i64 = 24576;

/// The two extension gains for the coming half-frame.
///
/// The first is always a half.  The second is raised to a half when the
/// spectrum is not peaky, the previous gain was already strong, and the lag has
/// stayed within three quarters and one and a half of what it was — that is,
/// when the pitch has carried over rather than restarted.  A peaky spectrum
/// clears it; anything else leaves it where it was.
pub fn interpolation_gains(
    reflection: i16,
    previous_gain: i16,
    lag: i16,
    previous_lag: i16,
) -> (i16, i16) {
    let peaky = acc(shift((reflection as i64) << 16, 3) - (HALF << 16)) > 0;
    let strong = acc(shift((previous_gain as i64) << 16, 1) - (WAS_STRONG << 16)) > 0;
    let stretched = acc(((lag as i64) << 16) - acc(SAME_PITCH * (previous_lag as i64) * 2)) > 0;
    let squeezed =
        acc(shift(acc((previous_lag as i64) * SAME_PITCH * 2), 1) - ((lag as i64) << 16)) > 0;

    // Only the peaky test clears the gain; failing one of the others leaves
    // whatever the previous half-frame carried over.
    let second = if peaky {
        0
    } else if strong && stretched && squeezed {
        HALF
    } else {
        previous_gain as i64
    };
    (low(HALF), low(second))
}

/// Build the subframe's excitation from the two codebooks.
///
/// The adaptive codebook's contribution is what the buffer already holds — the
/// history read back at the chosen lag — and the fixed codebook's is the
/// innovation.  Each is taken at its own gain, the fixed one counting double,
/// and the sum is rounded and saturated back into the buffer.
pub fn excite(excitation: &mut [i16], innovation: &[i16], adaptive: i16, fixed: i16) {
    for k in 0..crate::convolve::SUBFRAME {
        let shaped = acc(mul(fixed, innovation[k]));
        let mut total = acc(mul(adaptive, excitation[k]));
        total = acc(total + shift(shaped, 1));
        total = acc(shift(total, 1) + 0x8000);
        excitation[k] = hi(sat(total));
    }
}

/// Taps of the filter memories carried between subframes.
const MEMORY: usize = 10;

/// Update the filter memories after a subframe.
///
/// Two are kept: the difference between the input and what the synthesis filter
/// produced, and the weighted target with both codebooks' contributions taken
/// out.  Both are what the next subframe's filters start from.
pub fn update_memories(
    input: &[i16],
    synthesised: &[i16],
    weighted: &[i16],
    adaptive: &[i16],
    innovation: &[i16],
    gains: (i16, i16),
) -> ([i16; MEMORY], [i16; MEMORY]) {
    let (mut error, mut residual) = ([0i16; MEMORY], [0i16; MEMORY]);
    for k in 0..MEMORY {
        error[k] = hi(sat(acc(
            ((input[k] as i64) << 16) - ((synthesised[k] as i64) << 16)
        )));

        let mut total = (weighted[k] as i64) << 16;
        total = acc(total - shift(acc(mul(gains.0, adaptive[k])), 1));
        total = acc(total - shift(acc(mul(gains.1, innovation[k])), 2));
        residual[k] = hi(sat(total));
    }
    (error, residual)
}
