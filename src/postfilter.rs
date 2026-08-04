//! Adaptive postfilter.
//!
//! One subframe at a time:
//!
//! 1. two bandwidth-expanded copies of the subframe's LPC set are derived by
//!    pulling its LSFs towards a fixed bias vector by two different amounts —
//!    these are the numerator `A(z/g2)` and denominator `A(z/g1)`;
//! 2. the synthesised speech is run through `A(z/g2)` to get a residual;
//! 3. a long-term (pitch) postfilter is applied to that residual;
//! 4. the result goes back through `1/A(z/g1)`, is gain-normalised, tilt
//!    compensated, and finally matched in energy to the unfiltered speech.

use crate::fixed::{
    acc, divide, exp, hi, low, mul, mul32x16, normalised_reciprocal, sat, shift, trunc32,
};
use crate::lsp::{self, ORDER};
use crate::pitch::SUBFRAME;
use crate::synth;

/// Length of the impulse response the postfilter measures itself with.
pub const PROBE_LEN: usize = 20;

/// Bandwidth-expand an LSF set towards the fixed bias vector and convert back
/// to LPC.
pub fn weighted_lpc(lsp_in: &[i16; ORDER], weights: &[i16; ORDER]) -> [i16; ORDER + 1] {
    let lsf = lsp::lsp_to_lsf(lsp_in);
    let bias = &crate::tables::PF_LSF_BIAS[..ORDER];
    let weight = weights;

    let mut pulled = [0i16; ORDER];
    for i in 0..ORDER {
        let mut a = (lsf[i] as i64) << 16;
        a = acc(a - mul(lsf[i], weight[i]));
        a = acc(a + mul(bias[i], weight[i]));
        pulled[i] = hi(a);
    }
    lsp::lsp_to_lpc(&lsp::lsf_to_lsp(&pulled))
}

/// Inverse filter `A(z)`: `out[n] = s[n] + sum a[i] s[n-i]`.
///
/// `speech` holds `ORDER` samples of history followed by the block to filter.
pub fn inverse_filter(a: &[i16; ORDER + 1], speech: &[i16], out: &mut [i16]) {
    for n in 0..out.len() {
        let mut sum = 0i64;
        for i in 1..=ORDER {
            sum = acc(sum + mul(a[i], speech[ORDER + n - i]));
        }
        out[n] = hi(sat(acc(shift(sum, 3)
            + ((speech[ORDER + n] as i64) << 16)
            + 0x8000)));
    }
}

/// Impulse response of `A(z/g2)/A(z/g1)`, truncated to [`PROBE_LEN`] samples.
///
/// The reference gets it by feeding the numerator coefficients, zero padded, through
/// the denominator's synthesis filter.
pub fn probe_response(a1: &[i16; ORDER + 1], a2: &[i16; ORDER + 1]) -> [i16; PROBE_LEN] {
    let mut input = [0i16; PROBE_LEN];
    input[..ORDER + 1].copy_from_slice(a2);
    let mut probe = [0i16; PROBE_LEN];
    let mut mem = [0i16; ORDER];
    synth::synthesis(a1, &input, &mut probe, &mut mem);
    probe
}

/// Spectral tilt of the probe: minus its first normalised autocorrelation
/// coefficient, or zero when the signal is degenerate.
pub fn spectral_tilt(probe: &[i16; PROBE_LEN]) -> i16 {
    let mut energy = 0i64;
    for &v in probe.iter() {
        energy = acc(energy + (v as i64) * (v as i64) * 2);
    }
    if energy == 0 {
        return 0;
    }
    let e = exp(energy);
    let energy = shift(energy, e);

    let mut corr = 0i64;
    for k in 0..PROBE_LEN - 1 {
        corr = acc(corr + mul(probe[k], probe[k + 1]));
    }
    let corr = shift(corr, e);
    let magnitude = if corr < 0 { acc(-corr) } else { corr };
    let energy_hi = (hi(energy) as i64) << 16;
    if acc(magnitude - energy_hi) > 0 {
        return 0;
    }
    hi(acc(-divide(energy_hi, hi(corr))))
}

/// Scale a block so the short-term postfilter does not change its level.
pub fn normalise_gain(probe: &[i16; PROBE_LEN], block: &mut [i16; SUBFRAME]) {
    let mut sum = 0i64;
    for &v in probe.iter() {
        let v = (v as i64) << 16;
        sum = acc(sum + if v < 0 { -v } else { v });
    }
    if acc(shift(sum, 2) - (16384i64 << 16)) <= 0 {
        return;
    }
    let (quotient, e) = normalised_reciprocal(sum);
    let gain = shift(shift(quotient, e), -3);
    for v in block.iter_mut() {
        *v = hi(mul32x16(gain, *v));
    }
}

/// Tilt compensation `g * (1 + k z^-1)` applied to the postfiltered block.
///
/// `prev` is the sample immediately before the block.  The filter gain is
/// chosen so the compensation cannot amplify.
pub fn compensate_tilt(tilt: i16, prev: i16, block: &[i16; SUBFRAME], out: &mut [i16; SUBFRAME]) {
    let strength = if tilt > 0 { 6554i16 } else { 16384i16 };
    let scaled = acc((tilt as i64) * (strength as i64) * 2);
    let k = hi(shift(scaled, 1));

    let magnitude = if scaled < 0 { acc(-scaled) } else { scaled };
    let gain = divide(acc((32768i64 << 16) - magnitude), 16384);
    let combined = acc((hi(gain) as i64) * (k as i64) * 2);
    let direct = hi(gain);

    for n in 0..SUBFRAME {
        let previous = if n == 0 { prev } else { block[n - 1] };
        let a = mul32x16(combined, previous);
        let b = acc((block[n] as i64) * (direct as i64) * 2);
        out[n] = hi(sat(acc(a + shift(b, 1) + (1 << 15))));
    }
}

/// Automatic gain control: scale the postfiltered block back to the energy of
/// the unfiltered speech, with a first-order smoother across subframes.
///
/// `gain` carries the smoother's state from one subframe to the next.
pub fn agc(speech: &[i16; SUBFRAME], block: &mut [i16; SUBFRAME], gain: &mut i16) {
    /// Smoothing factor, Q15.
    const AGC_FAC: i16 = 32358;
    /// `1 - AGC_FAC`, Q15.
    const AGC_ONE_MINUS: i16 = 410;

    fn abs_sum(v: &[i16; SUBFRAME]) -> i64 {
        let mut s = 0i64;
        for &x in v.iter() {
            let x = (x as i64) << 16;
            s = acc(s + if x < 0 { -x } else { x });
        }
        s
    }

    let mut target = 0i16;
    let reference = abs_sum(speech);
    if reference != 0 {
        let e_ref = exp(reference);
        let scaled_ref = shift(shift(reference, e_ref), -1);

        let filtered = abs_sum(block);
        if filtered == 0 {
            // Nothing to scale, and the smoother restarts from zero.
            *gain = 0;
            return;
        }
        let e_filt = exp(filtered);
        let normalised = shift(filtered, e_filt);
        let (quotient, e_recip) = normalised_reciprocal(normalised);

        let ratio = shift(mul32x16(scaled_ref, hi(quotient)), e_recip);
        let ratio = sat(ratio);
        let scaled = mul32x16(ratio, AGC_ONE_MINUS);
        target = hi(sat(shift(shift(scaled, e_filt - e_ref + 1), -1)));
    }

    let mut running = (*gain as i64) << 16;
    for v in block.iter_mut() {
        let smoothed = acc((hi(running) as i64) * (AGC_FAC as i64) * 2);
        running = sat(acc(smoothed + ((target as i64) << 16)));
        let scaled = acc((hi(running) as i64) * (*v as i64) * 2);
        *v = hi(sat(acc(shift(scaled, 1) + 0x8000)));
    }
    *gain = hi(running);
}

/// State of the output reconstruction filter.
#[derive(Clone, Default)]
pub struct OutputFilter {
    /// Feed-forward delay lines, three samples per section.
    fir: [i16; 6],
    /// Feedback delay lines, three 32-bit values per section.
    iir: [i64; 6],
}

/// Coefficients of the feed-forward part, shared by both sections.
/// Coefficients of the feedback part, shared by both sections.
/// Final output scaling, Q15.
const OUT_SCALE: i16 = 19661;

impl OutputFilter {
    /// Rebuild the filter state from its flat memory image: three
    /// 16-bit and three 32-bit delay values per section.
    pub fn from_words(fir: &[i16; 6], iir: &[i16; 12]) -> Self {
        let mut f = OutputFilter::default();
        f.fir.copy_from_slice(fir);
        for k in 0..6 {
            f.iir[k] =
                (((iir[2 * k] as u16 as i64) << 16) | (iir[2 * k + 1] as u16 as i64)) as i32 as i64;
        }
        f
    }

    /// The same image back out, so a replayed filter's state can be compared
    /// against the reference's.
    #[doc(hidden)]
    pub fn to_words(&self) -> ([i16; 6], [i16; 12]) {
        let mut iir = [0i16; 12];
        for k in 0..6 {
            iir[2 * k] = hi(self.iir[k]);
            iir[2 * k + 1] = low(self.iir[k]);
        }
        (self.fir, iir)
    }

    /// Run one half-frame through the two cascaded sections.
    pub fn run(&mut self, block: &mut [i16]) {
        let ff = &crate::tables::OUT_FIR[..4];
        let fb = &crate::tables::OUT_IIR[..3];
        for sample in block.iter_mut() {
            let mut value = *sample;
            for section in 0..2 {
                let fir = &mut self.fir[section * 3..section * 3 + 3];
                let iir = &mut self.iir[section * 3..section * 3 + 3];

                let mut b = acc((value as i64) * (ff[0] as i64) * 2);
                for k in 0..3 {
                    b = acc(b + mul(fir[k], ff[k + 1]));
                }
                fir[2] = fir[1];
                fir[1] = fir[0];
                fir[0] = value;

                for k in 0..3 {
                    b = acc(b - shift(mul32x16(iir[k], fb[k]), 1));
                }
                iir[2] = iir[1];
                iir[1] = iir[0];
                iir[0] = trunc32(shift(b, 1));

                value = hi(sat(acc(shift(b, 2) + (1 << 15))));
            }
            let scaled = acc((value as i64) * (OUT_SCALE as i64) * 2);
            *sample = hi(sat(acc(shift(scaled, 1) + (1 << 15))));
        }
    }
}
