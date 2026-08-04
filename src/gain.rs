//! The gain codebook: the encoder's search over it and the decoder's read of
//! the index that comes back.
//!
//! Both subframe gains come from a single 7-bit index into a conjugate
//! two-stage codebook: the 3 high bits pick a row of the first stage, the 4 low
//! bits a row of the second, and the two rows are summed.  The first component
//! of the sum is the pitch gain; the second is a *correction factor* applied to
//! a predicted code gain, which is itself derived from the innovation's energy
//! and a 4-tap moving average of past gain energies.
//!
//! [`search`] is the encoder's side of the same table: the five
//! measures the closed loop leaves behind are brought onto one exponent by
//! [`align`] and every entry is scored against them.

use crate::fixed::{acc, exp, hi, low, mul, mul32x16, sat, shift, trunc32};
use crate::pitch::SUBFRAME;
use crate::tables::{GAIN_CB1, GAIN_CB2, LOG2, POW2};

/// Taps of the code-gain energy predictor.
const ENERGY_TAPS: usize = 4;
/// 20*log10(2) in Q12 — converts a log2 value to dB.
const DB_PER_OCTAVE: i16 = 24660;
/// The reciprocal of [`DB_PER_OCTAVE`] in Q15.
const OCTAVES_PER_DB: i16 = 5439;
/// Energy floor written into the predictor memory after a reset.
pub const ENERGY_RESET: i16 = -17254;

/// Gains in force for the current subframe.
#[derive(Clone, Copy, Debug, Default)]
pub struct Gains {
    /// Adaptive-codebook (pitch) gain, Q14.
    pub pitch: i16,
    /// Fixed-codebook gain.
    pub code: i16,
}

/// Persistent gain-decoder state.
#[derive(Clone)]
pub struct GainState {
    /// Log energies of the last four innovations, newest first.
    pub energy: [i16; ENERGY_TAPS],
    /// Gains decoded for the previous subframe.
    pub gains: Gains,
    /// Consecutive suppressed frames, used to ramp the gains down.
    pub erasures: i16,
    /// Steps left in the ramp that restores the gains after a long erasure.
    pub recovery: i16,
    /// Pitch gain of the previous subframe, clamped, used for pitch sharpening.
    pub sharpen_gain: i16,
    /// One-subframe-delayed copy of [`Self::sharpen_gain`].
    pub sharpen_gain_prev: i16,
}

impl Default for GainState {
    fn default() -> Self {
        GainState {
            energy: [ENERGY_RESET; ENERGY_TAPS],
            gains: Gains::default(),
            erasures: 0,
            recovery: 0,
            sharpen_gain: 3277,
            sharpen_gain_prev: 3277,
        }
    }
}

/// Base-2 logarithm of a positive accumulator.
///
/// Returns the fractional part as an accumulator and the negated exponent, so
/// that `log2(x) = -exponent + fraction`.
pub fn log2(value: i64) -> (i64, i16) {
    if value <= 0 {
        return (0, 0);
    }
    let scaled = shift(value, 5);
    let e = exp(scaled);
    let normalised = shift(scaled, e);
    let mut a = shift(acc(normalised - (16384i64 << 16)), 5);
    let index = hi(shift(a, -14)) as usize;
    a = acc(a - ((index as i64) << 30));
    let frac = hi(shift(a, -1));

    let b = (LOG2[index] as i64) << 16;
    let slope = shift(acc(b - ((LOG2[index + 1] as i64) << 16)), 2);
    let interpolated = acc(b - (hi(slope) as i64) * (frac as i64) * 2);
    ((hi(interpolated) as i64) << 16, -(e as i16))
}

/// `2^x` for a positive accumulator holding the exponent.
pub fn pow2(value: i64) -> i64 {
    let mut a = shift(value, 5);
    let index = hi(shift(a, -15)) as usize;
    a = acc(a - ((index as i64) << 31));
    let frac = hi(shift(a, -2));

    let b = (POW2[index] as i64) << 16;
    let slope = shift(acc(b - ((POW2[index + 1] as i64) << 16)), 3);
    let doubled = shift(b, 1);
    let interpolated = acc(doubled - (hi(slope) as i64) * (frac as i64) * 2);
    (hi(shift(interpolated, -1)) as i64) << 16
}

/// Push the log energy of a decoded gain into the predictor memory.
fn update_energy(state: &mut GainState, value: i64) {
    for k in (1..ENERGY_TAPS).rev() {
        state.energy[k] = state.energy[k - 1];
    }
    let (frac, neg_exp) = log2(value);
    let combined = sat(shift(acc(frac + ((neg_exp as i64) << 31)), -3));
    state.energy[0] = hi(mul32x16(combined, DB_PER_OCTAVE));
}

/// Predict the code gain from the innovation's energy and the gain history.
///
/// Returns the predicted gain and the exponent it has to be scaled by.
pub fn predict_code_gain(state: &GainState, code: &[i16; SUBFRAME]) -> (i16, i16) {
    let mut energy = 10i64 << 17;
    for &c in code.iter() {
        energy = acc(energy + (c as i64) * (c as i64) * 2);
    }
    let mean = (hi(shift(energy, -2)) as i64) << 16;

    let (frac, neg_exp) = log2(mean);
    let combined = sat(shift(acc(frac + ((neg_exp as i64) << 31)), -3));
    let db = shift(mul32x16(combined, DB_PER_OCTAVE), 3);
    let mut target = acc((18432i64 << 20) - acc(db - (19488i64 << 19)));

    let pred = &crate::tables::GAIN_PRED[..ENERGY_TAPS];
    let mut history = 0i64;
    for (&p, &e) in pred.iter().zip(state.energy.iter()) {
        history = acc(history + mul(p, e));
    }
    target = sat(shift(acc(target + shift(history, 6)), -5));

    let octaves = shift(mul32x16(target, OCTAVES_PER_DB), 5);
    let exponent = hi(shift(octaves, -13));
    let fraction = acc(octaves - ((exponent as i64) << 29));
    (hi(pow2((hi(shift(fraction, 2)) as i64) << 16)), exponent)
}

/// Decode both gains for a subframe.
pub fn decode(state: &mut GainState, index: i16, code: &[i16; SUBFRAME]) -> Gains {
    let stage1 = &GAIN_CB1[2 * ((index >> 4) as usize)..][..2];
    let stage2 = &GAIN_CB2[2 * ((index & 15) as usize)..][..2];

    let pitch = hi(((stage1[0] as i64) + (stage2[0] as i64)) << 16);
    let correction = hi(((stage1[1] as i64) + (stage2[1] as i64)) << 16);

    let (gcode0, exponent) = predict_code_gain(state, code);
    let scaled = shift(
        (correction as i64) * (gcode0 as i64) * 2,
        (exponent - 7) as i32,
    );
    let code_gain = hi(sat(shift(scaled, -1)));

    update_energy(state, (correction as i64) << 16);

    let mut gains = Gains {
        pitch,
        code: code_gain,
    };

    // Coming out of a long run of erasures the gains are restored over four
    // subframes rather than all at once.
    let mut ramp = if state.erasures >= 8 {
        4
    } else {
        state.recovery
    };
    state.erasures = 0;
    state.recovery = ramp;
    if ramp > 0 {
        ramp -= 1;
        state.recovery = ramp;
        let scale = (3277i64) << (3 - ramp);
        gains.pitch = hi(acc(scale * (gains.pitch as i64) * 2));
        gains.code = hi(acc(scale * (gains.code as i64) * 2));
    }
    state.gains = gains;
    gains
}

/// Gain handling for a suppressed frame: the gains decay and the energy memory
/// is pulled down towards the reset floor.
pub fn decode_suppressed(state: &mut GainState) -> Gains {
    let (factor, ceiling, floor) = if state.erasures - 4 >= 0 {
        (26542i16, 29491i16, 4096i16)
    } else {
        (31470i16, 32113i16, 0i16)
    };
    state.erasures += 1;
    let pitch = hi(mul(state.gains.pitch, factor).min((ceiling as i64) << 15));
    let code = hi(mul(state.gains.code, factor));
    state.gains = Gains { pitch, code };
    decay_energy(state, floor);
    state.gains
}

/// Pull the gain-energy memory towards the reset floor.
fn decay_energy(state: &mut GainState, offset: i16) {
    let mut sum = 0i64;
    for &e in state.energy.iter() {
        sum = acc(sum + ((e as i64) << 16));
    }
    let mut a = shift(shift(sum, -2), 1);
    a = acc(a - ((offset as i64) << 16));
    a = a.max((ENERGY_RESET as i64) << 17);
    for k in (1..ENERGY_TAPS).rev() {
        state.energy[k] = state.energy[k - 1];
    }
    state.energy[0] = hi(shift(a, -1));
}

/// Terms the five measures are worth relative to one another, before they are
/// brought to a common exponent.
///
/// The first two are energies and the last three correlations, so they carry
/// different numbers of doublings; the pitch gain's exponent enters twice for
/// the energy term and once for each of the cross terms.
fn bias(index: usize, gain_exponent: i64) -> i64 {
    match index {
        0 => 12,
        1 => 13,
        2 => -2 * gain_exponent,
        3 => 7 - gain_exponent,
        _ => 6 - gain_exponent,
    }
}

/// Measures taken to a common exponent, ready for the codebook search.
pub const MEASURES: usize = 5;

/// Bring the five measures onto one exponent.
///
/// Each is held as a mantissa and an exponent; the search that follows adds
/// them together, so they first have to be lined up.  A measure more than
/// sixteen bits below the largest has its mantissa shifted down and the rest of
/// the gap left in the shift, and one more than thirty-two bits below is
/// dropped altogether.
pub fn align(
    mantissa: &mut [i16; MEASURES],
    exponent: &[i16; MEASURES],
    gain_exponent: i16,
) -> [i16; MEASURES] {
    let mut adjusted = [0i64; MEASURES];
    for i in 0..MEASURES {
        adjusted[i] = acc((exponent[i] as i64) + bias(i, gain_exponent as i64));
    }
    let smallest = adjusted.iter().copied().fold(i64::MAX, i64::min);

    let mut shift_out = [0i16; MEASURES];
    for i in 0..MEASURES {
        let gap = acc((low(adjusted[i]) as i64) - (low(smallest) as i64)).min(32);
        if gap == 32 {
            mantissa[i] = 0;
            shift_out[i] = 0;
        } else if acc(gap - 16) <= 0 {
            shift_out[i] = -(gap as i16);
        } else {
            mantissa[i] = hi(shift((mantissa[i] as i64) << 16, (16 - gap) as i32));
            shift_out[i] = -16;
        }
    }
    shift_out
}

/// The gain codebook, split into a coarse table of eight entries and a fine one
/// of sixteen.  Each entry is a pair: one term for the adaptive contribution and
/// one for the fixed one.
const COARSE_ENTRIES: usize = 8;
const FINE_ENTRIES: usize = 16;

/// Search the gain codebook.
///
/// Every one of the 128 entries is scored against the five aligned measures.
/// The two gains an entry stands for are `p`, applied to the fixed codebook, and
/// `v`, the adaptive gain scaled by `scale`; the score is the coding error
/// written out in the five terms `p^2`, `p`, `v^2`, `v` and `p v`, each weighted
/// by its measure.  The entry with the smallest error wins.
pub fn search(
    scale: i16,
    mantissa: &[i16; MEASURES],
    shifts: &[i16; MEASURES],
    gain_shift: i16,
) -> (usize, i16, i16) {
    let coarse = &crate::tables::GAIN_CB1;
    let fine = &crate::tables::GAIN_CB2;

    let mut best = 0x7fff_ffffi64;
    let mut chosen = 0usize;
    for outer in 0..COARSE_ENTRIES {
        for inner in 0..FINE_ENTRIES {
            let p = hi(acc(
                ((coarse[2 * outer] as i64) + (fine[2 * inner] as i64)) << 16
            ));
            let q = hi(acc(((coarse[2 * outer + 1] as i64)
                + (fine[2 * inner + 1] as i64))
                << 16));

            let v = trunc32(acc(crate::fixed::mul(scale, q)));
            let square = trunc32(crate::fixed::square32(v));
            let cross = trunc32(crate::fixed::mul32x16(v, p));
            let energy = acc(crate::fixed::mul(p, p));

            let term = |value: i64, i: usize| shift(value, shifts[i] as i32);
            let mut cost = term(crate::fixed::mul32x16(energy, mantissa[0]), 0);
            cost = acc(cost + term(acc(crate::fixed::mul(p, mantissa[1])), 1));
            cost = acc(cost + term(crate::fixed::mul32x16(square, mantissa[2]), 2));
            cost = acc(cost + term(crate::fixed::mul32x16(v, mantissa[3]), 3));
            cost = acc(cost + term(crate::fixed::mul32x16(cross, mantissa[4]), 4));

            if cost < best {
                chosen = FINE_ENTRIES * outer + inner;
            }
            best = trunc32(best.min(cost));
        }
    }

    let (outer, inner) = (chosen / FINE_ENTRIES, chosen % FINE_ENTRIES);
    let p = hi(acc(
        ((coarse[2 * outer] as i64) + (fine[2 * inner] as i64)) << 16
    ));
    let q = hi(acc(((coarse[2 * outer + 1] as i64)
        + (fine[2 * inner + 1] as i64))
        << 16));
    let scaled = shift(acc(crate::fixed::mul(scale, q)), (gain_shift - 7) as i32);
    (chosen, p, hi(shift(scaled, -1)))
}
