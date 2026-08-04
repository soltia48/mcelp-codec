//! Line spectral frequencies: dequantisation, conversions and interpolation.
//!
//! Representations used here:
//!
//! * **LSF** — line spectral *frequencies* in radians, Q13 (`pi` = 25736).
//! * **LSP** — line spectral *pairs*, i.e. `cos(lsf)`, Q15.
//! * **LPC** — direct-form predictor `a[0..=10]`, Q12 with `a[0] = 4096`.
//!
//! The quantiser is a two-stage split vector quantiser feeding a third-order
//! moving-average predictor, with one of two predictor sets selected per frame.

use crate::fixed::{acc, hi, mul, mul32x16, sat, shift, trunc32};

/// Predictor order.
pub const ORDER: usize = 10;
/// Number of past frames the LSF moving-average predictor looks at.
pub const MA_TAPS: usize = 3;

/// Smallest allowed LSF, Q13.
const LSF_MIN: i16 = 41;
/// Largest allowed LSF, Q13.
const LSF_MAX: i16 = 25682;
/// Minimum spacing enforced between adjacent LSFs, Q13.
const LSF_GAP: i16 = 321;
/// The two spacings enforced on the quantiser output before prediction.
const REARRANGE_GAPS: [i16; 2] = [10, 5];

/// Indices carried by the two LSF fields of a frame.
#[derive(Clone, Copy, Debug)]
pub struct LsfIndices {
    /// Which of the two MA predictor sets to use.
    pub mode: usize,
    /// First-stage index, 7 bit.
    pub stage1: usize,
    /// Second-stage index for coefficients 0..4, 6 bit.
    pub stage2_low: usize,
    /// Second-stage index for coefficients 5..9, 6 bit.
    pub stage2_high: usize,
}

impl LsfIndices {
    /// Split the two received LSF fields into codebook indices.
    pub fn unpack(field0: i16, field1: i16) -> Self {
        LsfIndices {
            mode: ((field0 >> 7) & 1) as usize,
            stage1: (field0 & 0x7f) as usize,
            stage2_high: ((field1 >> 6) & 0x3f) as usize,
            stage2_low: (field1 & 0x3f) as usize,
        }
    }
}

/// Persistent state of the LSF decoder.
#[derive(Clone)]
pub struct LsfState {
    /// Quantiser output of the three previous frames, newest first.
    pub history: [[i16; ORDER]; MA_TAPS],
    /// Last decoded LSF vector, replayed when a frame is suppressed.
    pub prev_lsf: [i16; ORDER],
    /// Predictor set used by the last good frame.
    pub prev_mode: usize,
}

impl Default for LsfState {
    /// Reset state: the predictor memory and the "previous frame" are seeded
    /// with the flat mean LSF vector.
    fn default() -> Self {
        let mean = &crate::tables::LSF_MEAN[..ORDER];
        let mut init = [0i16; ORDER];
        init.copy_from_slice(mean);
        LsfState {
            history: [init; MA_TAPS],
            prev_lsf: init,
            prev_mode: 0,
        }
    }
}

/// Look up one row of a codebook laid out as `rows x ORDER`.
fn codebook_row(base: &'static [i16], index: usize) -> &'static [i16] {
    &base[index * ORDER..][..ORDER]
}

/// Decode one frame's LSF vector.
///
/// Sum the two quantiser stages, force a minimum spacing on the
/// result, run it through the MA predictor, push it into the predictor memory
/// and finally stabilise the reconstructed LSFs.
pub fn decode(state: &mut LsfState, idx: LsfIndices) -> [i16; ORDER] {
    let stage1 = codebook_row(&crate::tables::LSF_CB1, idx.stage1);
    let low = codebook_row(&crate::tables::LSF_CB2, idx.stage2_high);
    let high = codebook_row(&crate::tables::LSF_CB2, idx.stage2_low);

    let mut residual = [0i16; ORDER];
    for i in 0..ORDER / 2 {
        residual[i] = hi(((stage1[i] as i64) + (low[i] as i64)) << 16);
        residual[i + 5] = hi(((stage1[i + 5] as i64) + (high[i + 5] as i64)) << 16);
    }

    for gap in REARRANGE_GAPS {
        rearrange(&mut residual, gap);
    }

    let lsf = ma_predict(&residual, &state.history, idx.mode);
    push_history(&mut state.history, &residual);

    let mut lsf = lsf;
    stabilize(&mut lsf);
    state.prev_mode = idx.mode;
    lsf
}

/// Recover the quantiser residual that a given LSF vector implies.
///
/// The moving-average prediction is subtracted off and the remainder divided by
/// the predictor's residual gain — spelled as a multiply by its reciprocal.
/// The encoder uses this to find what it has to quantise; the decoder uses it
/// to keep the predictor memory consistent through an erasure.
pub fn residual(
    lsf: &[i16; ORDER],
    history: &[[i16; ORDER]; MA_TAPS],
    mode: usize,
) -> [i16; ORDER] {
    let pred = &crate::tables::LSF_MA_PRED
        [mode * MA_TAPS * ORDER..mode * MA_TAPS * ORDER + MA_TAPS * ORDER];
    let gain_inv = &crate::tables::LSF_MA_GAIN_INV[mode * ORDER..mode * ORDER + ORDER];
    let mut out = [0i16; ORDER];
    for i in 0..ORDER {
        let mut a = (lsf[i] as i64) << 16;
        for j in 0..MA_TAPS {
            a = acc(a - mul(history[j][i], pred[j * ORDER + i]));
        }
        out[i] = hi(shift(mul32x16(acc(a), gain_inv[i]), 3));
    }
    out
}

/// Rebuild the LSF vector for a suppressed frame.
///
/// The erasure branch: the previous LSF vector is
/// repeated, and the quantiser residual that *would* have produced it is
/// recovered by inverting the predictor so the MA memory stays consistent.
pub fn decode_suppressed(state: &mut LsfState) -> [i16; ORDER] {
    let lsf = state.prev_lsf;
    let residual = residual(&lsf, &state.history, state.prev_mode);
    push_history(&mut state.history, &residual);
    lsf
}

/// Force a minimum spacing of `gap` between neighbouring LSFs by pushing the
/// offending pair symmetrically apart.
pub fn rearrange(lsf: &mut [i16; ORDER], gap: i16) {
    for k in 0..ORDER - 1 {
        let cur = (lsf[k] as i64) << 16;
        let excess = shift(cur - ((lsf[k + 1] as i64) << 16) + ((gap as i64) << 16), -1);
        if excess > 0 {
            lsf[k] = hi(cur - excess);
            lsf[k + 1] = hi(excess + ((lsf[k + 1] as i64) << 16));
        }
    }
}

/// Reconstruct LSFs from the quantiser residual and the predictor memory.
fn ma_predict(
    residual: &[i16; ORDER],
    history: &[[i16; ORDER]; MA_TAPS],
    mode: usize,
) -> [i16; ORDER] {
    let pred = &crate::tables::LSF_MA_PRED
        [mode * MA_TAPS * ORDER..mode * MA_TAPS * ORDER + MA_TAPS * ORDER];
    let gain = &crate::tables::LSF_MA_GAIN[mode * ORDER..mode * ORDER + ORDER];

    let mut lsf = [0i16; ORDER];
    for i in 0..ORDER {
        let mut a = mul(residual[i], gain[i]);
        for j in 0..MA_TAPS {
            a = acc(a + mul(history[j][i], pred[j * ORDER + i]));
        }
        lsf[i] = hi(a);
    }
    lsf
}

/// Shift the predictor memory and insert this frame's residual.
fn push_history(history: &mut [[i16; ORDER]; MA_TAPS], residual: &[i16; ORDER]) {
    for j in (1..MA_TAPS).rev() {
        history[j] = history[j - 1];
    }
    history[0] = *residual;
}

/// Make an LSF vector monotonic and keep it inside the usable band.
pub fn stabilize(lsf: &mut [i16; ORDER]) {
    // One bubble pass: swap any inverted neighbours.
    for k in 0..ORDER - 1 {
        if lsf[k + 1] < lsf[k] {
            lsf.swap(k, k + 1);
        }
    }
    if lsf[0] < LSF_MIN {
        lsf[0] = LSF_MIN;
    }
    for k in 0..ORDER - 1 {
        if (lsf[k + 1] as i64) - (lsf[k] as i64) < LSF_GAP as i64 {
            lsf[k + 1] = hi(((lsf[k] as i64) << 16) + ((LSF_GAP as i64) << 16));
        }
    }
    if lsf[ORDER - 1] > LSF_MAX {
        lsf[ORDER - 1] = LSF_MAX;
    }
}

/// `cos()` of an LSF vector: LSF (Q13 radians) -> LSP (Q15), by table lookup
/// with linear interpolation.
pub fn lsf_to_lsp(lsf: &[i16; ORDER]) -> [i16; ORDER] {
    let mut lsp = [0i16; ORDER];
    for i in 0..ORDER {
        // 20861 in Q15, doubled by the fractional multiply, is 4/pi: it
        // scales a Q13 radian value to a 6.8 fixed-point table index.
        let scaled = mul(lsf[i], 20861);
        let index = hi(shift(scaled, -8).min(63 << 16)) as usize;
        let frac = (hi(scaled) & 0xff) as i64;
        let interp = acc(((crate::tables::COS[index] as i64) << 16)
            + (frac << 3) * (crate::tables::COS_SLOPE[index] as i64) * 2);
        lsp[i] = hi(interp);
    }
    lsp
}

/// `acos()` of an LSP vector: LSP (Q15) -> LSF (Q13 radians).
///
/// Walks the cosine table downwards from the top, so the descending LSP values
/// are matched by a single sweep.
pub fn lsp_to_lsf(lsp: &[i16; ORDER]) -> [i16; ORDER] {
    const TABLE_LAST: usize = 63;
    let mut lsf = [0i16; ORDER];
    // The sweep position carries over between coefficients, so two neighbours
    // may land on the same table entry.
    let mut index = TABLE_LAST;
    for i in (0..ORDER).rev() {
        let mut diff;
        loop {
            diff = (lsp[i] as i64) - (crate::tables::COS[index] as i64);
            if diff <= 0 || index == 0 {
                break;
            }
            index -= 1;
        }
        let slope = crate::tables::ACOS[index] as i64;
        // index.frac in 6.9 format, then scaled by pi (Q13) with rounding.
        let pos = acc(((index as i64) << 25) + ((diff * slope * 2) << 4));
        let angle = acc((pos >> 16) * 25736 * 2 + 0x8000);
        lsf[i] = hi(angle);
    }
    lsf
}

/// Symmetric/antisymmetric polynomial of five LSPs.
///
/// Produces `f[0..=5]` in Q24; `lsp` is read with a stride of two.
fn lsp_polynomial(lsp: &[i16], first: usize) -> [i64; 6] {
    let mut f = [0i64; 6];
    f[0] = 1 << 24;
    f[1] = trunc32(-((lsp[first] as i64) << 10));
    for i in 2..=5 {
        let l = -lsp[first + 2 * (i - 1)];
        f[i] = trunc32(acc(2 * f[i - 2] + (mul32x16(f[i - 1], l) << 1)));
        for j in (2..i).rev() {
            f[j] = trunc32(acc(f[j] + (mul32x16(f[j - 1], l) << 1) + f[j - 2]));
        }
        f[1] = trunc32(acc(f[1] + ((l as i64) << 10)));
    }
    f
}

/// LSP -> direct-form LPC.
///
/// `A(z) = (F1(z)(1 + z^-1) + F2(z)(1 - z^-1)) / 2` with `F1`/`F2` built from
/// the even and odd LSPs.
pub fn lsp_to_lpc(lsp: &[i16; ORDER]) -> [i16; ORDER + 1] {
    let mut f1 = lsp_polynomial(lsp, 0);
    let mut f2 = lsp_polynomial(lsp, 1);
    for i in (1..=5).rev() {
        f1[i] = trunc32(acc(f1[i] + f1[i - 1]));
        f2[i] = trunc32(acc(f2[i] - f2[i - 1]));
    }

    let mut a = [0i16; ORDER + 1];
    a[0] = 4096;
    for i in 1..=5 {
        a[i] = hi(acc(shift(shift(acc(f1[i] + f2[i]), -1), 4) + (1 << 15)));
        a[ORDER + 1 - i] = hi(acc(shift(shift(acc(f1[i] - f2[i]), -1), 4) + (1 << 15)));
    }
    a
}

/// Build the two LSP sets a half-frame needs and convert both to LPC.
///
/// The first subframe of the pair uses the midpoint between the previous and
/// current LSP sets, the second uses the current set unchanged.
pub fn interpolate_pair(
    prev: &[i16; ORDER],
    cur: &[i16; ORDER],
) -> (
    [i16; ORDER],
    [i16; ORDER],
    [i16; ORDER + 1],
    [i16; ORDER + 1],
) {
    let mid = midpoint(prev, cur);
    (mid, *cur, lsp_to_lpc(&mid), lsp_to_lpc(cur))
}

/// Midpoint of two line spectral pairs.
///
/// The subframe filters are built from the two half-frame vectors and the point
/// halfway between them; this is that midpoint.
pub fn midpoint(first: &[i16; ORDER], second: &[i16; ORDER]) -> [i16; ORDER] {
    let mut out = [0i16; ORDER];
    for i in 0..ORDER {
        // The sum, halved as it is stored.
        out[i] = hi(acc(((first[i] as i64) + (second[i] as i64)) << 15));
    }
    out
}

/// LPC -> reflection coefficients, by running the Levinson recursion backwards.
///
/// Each step peels the last coefficient off as a reflection coefficient and
/// divides the rest by `1 - k^2`, updating the array from both ends at once.
pub fn reflection_coefficients(lpc: &[i16; ORDER + 1]) -> [i16; ORDER] {
    let mut k = [0i16; ORDER];
    k.copy_from_slice(&lpc[1..]);

    for step in 0..ORDER - 1 {
        let order = ORDER - step;
        let top = order - 1;
        let kv = k[top];

        // 1 - k^2 in Q30, then its reciprocal.
        let residual = acc((16384i64 << 16) - shift(acc((kv as i64) * (kv as i64) * 2), 5));
        let (quotient, exponent) = crate::fixed::normalised_reciprocal(residual);
        let scale = (exponent << 27) >> 27; // a signed 5-bit shift field
        let inverse = hi(shift(quotient, -1));

        let pairs = ((order - 2) >> 1) + 1;
        for p in 0..pairs {
            let (i, j) = (p, top - 1 - p);
            let low = acc(((k[i] as i64) << 16) - shift(acc(mul(k[j], kv)), 3));
            let updated_low = hi(shift(mul32x16(sat(low), inverse), scale));
            let high = acc(((k[j] as i64) << 16) - shift(acc(mul(kv, k[i])), 3));
            let updated_high = hi(shift(mul32x16(sat(high), inverse), scale));
            k[i] = updated_low;
            k[j] = updated_high;
        }
    }
    k
}

/// Decide whether this frame's LSP sets may be interpolated across the two
/// half-frames.
///
/// Interpolation is vetoed while the spectrum looks flat: either a large first
/// reflection coefficient or a large prediction-error product is enough.
/// Returns the new flag together with whether the current LSP set should be
/// used for *both* halves, which happens on the transition into the peaky
/// state.
pub fn interpolation_control(reflection: &[i16; ORDER], previous_flag: bool) -> (bool, bool) {
    let term = |k: i16| acc((16384i64 << 16) - shift(acc((k as i64) * (k as i64) * 2), 5));
    let mut product = sat(term(reflection[0]));
    for &k in reflection[1..].iter() {
        product = mul32x16(product, hi(term(k)));
    }
    let error = hi(shift(product, 9));

    let peaky = acc(shift((reflection[0] as i64) << 16, 3) - (16384i64 << 16)) > 0
        || acc(shift((error as i64) << 16, 1) - (16384i64 << 16)) > 0;
    if peaky {
        (false, false)
    } else {
        (true, !previous_flag)
    }
}

/// Level the target is compared against when counting sign changes, and the
/// two halves the count is taken over.
const CROSSING_GAIN: i16 = 26214;
const FIRST_HALF: usize = 40;

/// Sign changes below which the frame is not treated as noisy at all, and above
/// which it is treated as fully noisy.
const QUIET: i64 = 10;
const NOISY: i64 = 20;

/// Frames of start-up during which the slower pair of weights is used.
const SETTLED: i16 = 50;

/// Interpolation weights and squared-distance thresholds, before and after the
/// start-up period.  Each pair is a weight and its complement against unity.
struct Weights {
    slow: [i16; 2],
    fast: [i16; 2],
    near: i64,
    far: i64,
}

const STARTING: Weights = Weights {
    slow: [5898, 26870],
    fast: [13435, 19333],
    near: 0x000d_4563,
    far: 0x0018_9375,
};
const RUNNING: Weights = Weights {
    slow: [27853, 4915],
    fast: [30147, 2621],
    near: 0x0000_29f1,
    far: 0x0043_9581,
};

/// How the line spectrum is to be quantised.
///
/// The encoder keeps two slowly moving copies of the line spectrum, one tracking
/// it closely and one lagging further behind, and decides the mode from how far
/// the current set has moved away from them.  A set that has barely moved can be
/// coded against the running average; one that has jumped needs the full
/// codebook.  Two spectral measures taken from the reflection coefficients then
/// override the decision when the spectrum is flat enough that the residual is
/// noise-like.
#[derive(Clone, Copy, Default)]
pub struct Quantiser {
    /// Frames since the start, up to [`SETTLED`].
    age: i16,
    /// The closely tracking copy.
    near: [i16; ORDER],
    /// The lagging copy.
    far: [i16; ORDER],
}

impl Quantiser {
    /// Choose the mode for one half-frame.
    pub fn choose(&mut self, lpc: &[i16; ORDER + 1], target: &[i16], lsf: &[i16; ORDER]) -> i16 {
        let mut reflection = reflection_coefficients(lpc);

        // Prediction error, as a fraction of the signal's own energy.
        let term = |k: i16| acc((16384i64 << 16) - shift(acc((k as i64) * (k as i64) * 2), 5));
        let mut product = sat(term(reflection[0]));
        for &k in reflection[1..].iter() {
            product = crate::fixed::mul32x16(product, hi(term(k)));
        }
        let error = hi(shift(product, 9));

        // A noisy frame is talked up towards flatness so that the tests at the
        // end catch it even when the reflection coefficient alone would not.
        let (first, second) = crossings(target);
        if first > QUIET && second > QUIET {
            reflection[0] = if first > NOISY && second > NOISY {
                hi(shift(
                    acc(shift((reflection[0] as i64) << 16, 3) + (16384i64 << 16)),
                    -3,
                ))
            } else {
                let excess = acc((first + second - NOISY) << 16);
                let graded = shift(acc(crate::fixed::mul(hi(excess), CROSSING_GAIN)), 7);
                hi(acc(graded + ((reflection[0] as i64) << 16)))
            };
        }

        let weights = if self.age < SETTLED {
            self.age += 1;
            &STARTING
        } else {
            &RUNNING
        };

        let previous = self.near;
        for i in 0..ORDER {
            self.near[i] = hi(acc(crate::fixed::mul(previous[i], weights.slow[0])
                + crate::fixed::mul(lsf[i], weights.slow[1])));
        }

        let moved = acc(distance(&self.near, &previous) - weights.near) >= 0
            && acc(distance(lsf, &self.far) - weights.far) >= 0;

        let mut mode = 1;
        if moved {
            if acc(shift((error as i64) << 16, 1) - (16384i64 << 16)) <= 0 {
                mode = 2;
            }
        } else {
            for (far, &current) in self.far.iter_mut().zip(lsf.iter()) {
                *far = hi(acc(crate::fixed::mul(*far, weights.fast[0])
                    + crate::fixed::mul(current, weights.fast[1])));
            }
        }

        if acc(shift((reflection[0] as i64) << 16, 3) - (16384i64 << 16)) > 0 {
            mode = 0;
        }
        mode
    }
}

/// Sign changes of the target about its own mean level, counted separately over
/// the two halves of the subframe.
fn crossings(target: &[i16]) -> (i64, i64) {
    let mut total = 0i64;
    for &x in target.iter() {
        total = acc(total + crate::fixed::mul(x, CROSSING_GAIN));
    }
    let level = hi(shift(total, -6));

    let mut counts = [0i64; 2];
    let mut previous = hi(acc(((target[0] as i64) - (level as i64)) << 16));
    for (i, &x) in target.iter().enumerate().skip(1) {
        let difference = hi(acc(((x as i64) - (level as i64)) << 16));
        if acc(crate::fixed::mul(previous, difference)) < 0 {
            counts[usize::from(i >= FIRST_HALF)] += 1;
        }
        previous = difference;
    }
    (counts[0], counts[1])
}

/// Squared distance between two line spectra.
fn distance(x: &[i16], y: &[i16]) -> i64 {
    let mut total = 0i64;
    for i in 0..ORDER {
        let difference = hi(acc(((x[i] as i64) - (y[i] as i64)) << 16));
        total = acc(total + crate::fixed::mul(difference, difference));
    }
    total
}

/// Accessors used by the trace replay to line the state up with the reference.
impl Quantiser {
    pub fn restore(&mut self, near: [i16; ORDER], far: [i16; ORDER], age: i16) {
        self.near = near;
        self.far = far;
        self.age = age;
    }
    pub fn saved(&self) -> ([i16; ORDER], [i16; ORDER], i16) {
        (self.near, self.far, self.age)
    }
}
