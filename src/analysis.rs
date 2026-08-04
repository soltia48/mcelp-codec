//! Spectral analysis front end.
//!
//! Every half-frame the encoder assembles a 256-sample analysis window out of
//! 96 samples of history and the 160 new ones, normalises it, tapers both ends,
//! pre-emphasises it and hands the result to a 128-point FFT.  The magnitude
//! spectrum that comes back drives the LPC/LSF estimation.

use crate::fixed::{acc, exp, hi, round, sat, shift};
use crate::preprocess::HALF;
use crate::tables::SINE as SINE_TABLE;

/// Samples in one analysis window.
pub const WINDOW: usize = 256;
/// Samples of the previous half-frame the window reaches back over.
pub const OVERLAP: usize = WINDOW - HALF;
/// Length of the taper applied to each end of the window.
const TAPER: usize = 10;
/// Largest left shift the normaliser will apply.
const MAX_HEADROOM: i16 = 8;
/// Pre-emphasis coefficient, Q15 (0.8).
const PREEMPHASIS: i16 = 26214;

/// State the analysis front end carries between half-frames.
#[derive(Clone)]
pub struct Analysis {
    /// Tail of the previous half-frame.
    pub history: [i16; OVERLAP],
    /// Last windowed sample, kept for the pre-emphasis of the next window.
    pub last_windowed: i16,
    /// Headroom the previous window was normalised by.
    pub prev_headroom: i16,
}

impl Default for Analysis {
    fn default() -> Self {
        Analysis {
            history: [0; OVERLAP],
            last_windowed: 0,
            prev_headroom: 0,
        }
    }
}

impl Analysis {
    /// Build, normalise and taper the analysis window.
    ///
    /// Returns the tapered window and the headroom it was normalised by.
    pub fn window(&mut self, block: &[i16; HALF]) -> ([i16; WINDOW], i16) {
        let mut buf = [0i16; WINDOW];
        buf[..OVERLAP].copy_from_slice(&self.history);
        buf[OVERLAP..].copy_from_slice(block);
        self.history.copy_from_slice(&block[HALF - OVERLAP..]);

        // Normalise the whole window against its largest sample, but never by
        // more than eight bits.
        let mut peak = 0i64;
        for &v in buf.iter() {
            // `ABS` saturates, so a sample at the negative rail peaks at full
            // scale rather than one count past it.
            let m = sat(((v as i64) << 16).abs());
            if m > peak {
                peak = m;
            }
        }
        let headroom = if peak == 0 {
            MAX_HEADROOM
        } else {
            (exp(peak) as i16).min(MAX_HEADROOM)
        };
        for v in buf.iter_mut() {
            *v = hi(shift((*v as i64) << 16, headroom as i32));
        }

        // Taper both ends and drop eight bits so the FFT cannot overflow.
        let taper = &crate::tables::TAPER[..TAPER];
        let mut out = [0i16; WINDOW];
        for n in 0..WINDOW {
            out[n] = if n < TAPER {
                hi(shift(acc((buf[n] as i64) * (taper[n] as i64) * 2), -8))
            } else if n >= WINDOW - TAPER {
                hi(shift(
                    acc((buf[n] as i64) * (taper[WINDOW - 1 - n] as i64) * 2),
                    -8,
                ))
            } else {
                hi(shift((buf[n] as i64) << 16, -8))
            };
        }
        (out, headroom)
    }

    /// Pre-emphasise the tapered window in place.
    ///
    /// The sample carried over from the previous window is re-scaled first, so
    /// that a change of normalisation between windows does not create a step.
    pub fn preemphasise(&mut self, window: &mut [i16; WINDOW], headroom: i16) {
        let realign = -((headroom - self.prev_headroom) as i32).abs();
        let mut previous = hi(shift((self.last_windowed as i64) << 16, realign));
        self.last_windowed = window[WINDOW - 1];

        for sample in window.iter_mut() {
            let x = *sample;
            let a = acc(((x as i64) << 16) - (previous as i64) * (PREEMPHASIS as i64) * 2);
            previous = x;
            *sample = hi(a);
        }
    }
}

/// Points in the complex transform: the 256 real samples are packed as 128
/// complex ones, even samples in the real part and odd ones in the imaginary.
pub const FFT_POINTS: usize = 128;
/// Sine table; the cosine of the same angle is 32 entries further on.
/// Offset from sine to cosine in that table.
const QUARTER: usize = 32;

/// Unique bins a 256-point real transform produces.
pub const BINS: usize = FFT_POINTS + 1;

/// One complex spectrum, split the way the reference stores it.  The buffers hold one
/// entry more than the transform needs, because unpacking the real spectrum
/// writes a mirrored bin at each end.
pub struct Spectrum {
    /// Real parts.
    pub re: [i16; BINS],
    /// Imaginary parts.
    pub im: [i16; BINS],
}

/// Split a real block into the complex sequence the transform works on.
pub fn deinterleave(block: &[i16; WINDOW]) -> Spectrum {
    let mut s = Spectrum {
        re: [0; BINS],
        im: [0; BINS],
    };
    for n in 0..FFT_POINTS {
        s.re[n] = block[2 * n];
        s.im[n] = block[2 * n + 1];
    }
    s
}

/// Decimation-in-frequency FFT with bit-reversed output reordering.
///
/// Sums are stored without rounding and are allowed to wrap, which is why the
/// input is scaled down by 256 first; only the twiddled outputs are rounded.
pub fn fft(s: &mut Spectrum) {
    let sin = |k: usize| crate::tables::SINE[k] as i64;
    let cos = |k: usize| crate::tables::SINE[QUARTER + k] as i64;

    let mut span = FFT_POINTS / 2;
    let mut blocks = 1usize;
    while span >= 1 {
        for group in 0..span {
            let step = blocks;
            for k in 0..blocks {
                let i = group + 2 * span * k;
                let j = i + span;
                let (re_i, re_j) = (s.re[i] as i64, s.re[j] as i64);
                let (im_i, im_j) = (s.im[i] as i64, s.im[j] as i64);
                let dre = hi((re_i - re_j) << 16) as i64;
                let dim = im_i - im_j;

                s.re[i] = hi((re_i + re_j) << 16);
                s.im[i] = hi((im_i + im_j) << 16);

                let w = group * step;
                s.re[j] = hi(round(acc(cos(w) * dre * 2 + sin(w) * dim * 2)));
                s.im[j] = hi(round(acc(cos(w) * dim * 2 - sin(w) * dre * 2)));
            }
        }
        span /= 2;
        blocks *= 2;
        if blocks > FFT_POINTS {
            break;
        }
    }
    bit_reverse(s);
}

/// Reorder the transform output, by the reverse-carry address arithmetic the
/// reference uses.
fn bit_reverse(s: &mut Spectrum) {
    let mut j = 0u16;
    for i in 0..FFT_POINTS - 1 {
        if (i as u16) < j {
            s.re.swap(i, j as usize);
            s.im.swap(i, j as usize);
        }
        j = bitrev16(bitrev16(j).wrapping_add(bitrev16(FFT_POINTS as u16 / 2)));
    }
}

/// Reverse all sixteen bits of a value.
fn bitrev16(v: u16) -> u16 {
    v.reverse_bits()
}

/// Step of the rotator: `cos(pi/128)` and `-sin(pi/128)`.
const ROTATOR_STEP: (i16, i16) = (32758, -804);
/// Bins between two re-seeds of the rotator.
const ANCHOR_PERIOD: usize = 32;

/// Which way the rotator turns.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Pack the 128-point complex transform into the 129 real bins.
    Forward,
    /// The conjugate operation, used on the way back to the time domain.
    Inverse,
}

/// Unpack the 128-point complex transform of the even/odd split into the 129
/// unique bins of the 256-point real transform.
///
/// The rotator is advanced by a complex multiply and renormalised at every
/// step, and re-anchored from the sine table every 32 bins.  Both directions
/// share this code; they differ only in which way the rotator turns and in
/// whether bin 128 has to be filled in first.
pub fn unpack_real(s: &mut Spectrum, direction: Direction) {
    let forward = direction == Direction::Forward;
    let mut sin = 0i16;
    let mut cos = if forward { 32767 } else { -32768 };
    // The rotator's first re-seed is one step into the table, not at its head.
    let mut anchor = crate::tables::SINE_ANCHOR_STEP;
    let mut countdown = ANCHOR_PERIOD - 1;

    if forward {
        // The forward direction wraps bin 128 round onto bin 0; on the way
        // back the caller has already supplied all 129 bins.
        s.re[FFT_POINTS] = s.re[0];
        s.im[FFT_POINTS] = s.im[0];
    }

    for k in 0..=FFT_POINTS / 2 {
        let (lo, hi_idx) = (k, FFT_POINTS - k);
        let sum_re = (s.re[lo] as i64) + (s.re[hi_idx] as i64);
        let diff_re = hi(((s.re[lo] as i64) - (s.re[hi_idx] as i64)) << 16) as i64;
        let sum_im = hi(((s.im[lo] as i64) + (s.im[hi_idx] as i64)) << 16);
        let diff_im = (s.im[lo] as i64) - (s.im[hi_idx] as i64);

        let p = hi(round(acc(
            (sum_im as i64) * (cos as i64) * 2 + diff_re * (sin as i64) * 2
        ))) as i64;
        let q = hi(round(acc(
            (sum_im as i64) * (sin as i64) * 2 - diff_re * (cos as i64) * 2
        ))) as i64;

        s.im[hi_idx] = hi(shift(acc((q - diff_im) << 16), -1));
        s.im[lo] = hi(shift(acc((diff_im + q) << 16), -1));
        s.re[hi_idx] = hi(shift(acc((sum_re - p) << 16), -1));
        s.re[lo] = hi(shift(acc((sum_re + p) << 16), -1));

        if countdown == 0 {
            // Re-anchor straight from the table rather than letting the
            // rotator drift.  The forward transform turns clockwise, so only
            // the sine is negated; the inverse turns the other way and negates
            // both.
            sin = negate(SINE_TABLE[anchor]);
            cos = SINE_TABLE[anchor + QUARTER];
            if !forward {
                cos = negate(cos);
            }
            anchor += crate::tables::SINE_ANCHOR_STEP;
            countdown = ANCHOR_PERIOD - 1;
        } else {
            countdown -= 1;
            let (new_sin, new_cos) = advance(sin, cos, forward);
            let (n_sin, n_cos) = renormalise(new_sin, new_cos);
            sin = n_sin;
            cos = n_cos;
        }
    }
}

/// Negate a word in a saturating accumulator.
pub(crate) fn negate(v: i16) -> i16 {
    hi(sat(-((v as i64) << 16)))
}

/// Rotate the twiddle by one bin, in whichever direction is being used.
///
/// The rounding the accumulate-and-round steps apply saturates, which matters
/// at the very last
/// bin of a quarter turn: the sine there lands just outside the 32-bit range
/// and has to pin at -1 rather than wrap round to +1.
fn advance(sin: i16, cos: i16, forward: bool) -> (i16, i16) {
    let (cs, sn) = ROTATOR_STEP;
    let turn = if forward { 1i64 } else { -1 };
    let s = hi(round(sat(
        (cs as i64) * (sin as i64) * 2 + turn * (sn as i64) * (cos as i64) * 2
    )));
    let c = hi(round(sat(
        (cs as i64) * (cos as i64) * 2 - turn * (sn as i64) * (sin as i64) * 2
    )));
    (s, c)
}

/// Pull the rotator back onto the unit circle with one Newton step.
///
/// The analysis front end runs with overflow mode enabled, so the
/// sum below saturates instead of wrapping; that saturation is what keeps the
/// correction at unity when the rotator is already normalised.
fn renormalise(sin: i16, cos: i16) -> (i16, i16) {
    // With overflow mode on, each of the three steps saturates rather than
    // wrapping; a rotator that has drifted outwards therefore pins at 1.0
    // instead of flipping sign.
    let mut n = sat((cos as i64) * (cos as i64) * 2);
    n = sat(acc(n + (sin as i64) * (sin as i64) * 2));
    let norm = hi(sat(acc(n + (1 << 15))));
    let quotient = crate::analysis::ratio(16384, norm);
    let e = exp((norm as i64) << 16);
    let scaled = shift((quotient as i64) << 16, e);
    let sum = sat(acc((16384i64 << 16) + ((hi(scaled) as i64) << 16)));
    let mid = sat(acc(shift(sum, -1) + (1 << 15)));
    let apply = |v: i16| {
        let prod = acc((v as i64) * (hi(mid) as i64) * 2) & !0xffff;
        hi(shift(prod, 1))
    };
    (apply(sin), apply(cos))
}

/// `numerator / denominator`, clipped, as the reference computes it.
pub fn ratio(numerator: i16, denominator: i16) -> i16 {
    if numerator < 0 || denominator < 0 || numerator > denominator {
        return 0;
    }
    if numerator == denominator {
        return 32767;
    }
    crate::fixed::low(crate::fixed::restoring_divide(
        (numerator as i64) << 16,
        denominator,
        15,
    ))
}

/// Shortest lag the open-loop pitch measure considers.
const CORR_MIN_LAG: usize = 16;
/// Number of lags it scans.
const CORR_LAGS: usize = 113;

/// Normalised peak autocorrelation of the analysis window.
///
/// Returns a voicing measure: the largest correlation over lags 16..=128
/// divided by the window's energy, or zero if the window is silent.  The
/// correlation length shrinks as the lag grows so both operands stay inside
/// the window.
pub fn voicing(window: &[i16; WINDOW]) -> i16 {
    let mut energy = 0i64;
    for &v in window.iter() {
        energy = sat(acc(energy + (v as i64) * (v as i64) * 2));
    }
    let energy = hi(sat(shift(energy, 8)));
    if energy == 0 {
        return 0;
    }

    let mut peak = 0i64;
    for j in 0..CORR_LAGS {
        let lag = CORR_MIN_LAG + j;
        let taps = WINDOW - lag;
        let mut c = 0i64;
        for i in 0..taps {
            c = sat(acc(c + (window[i] as i64) * (window[lag + i] as i64) * 2));
        }
        let c = sat(shift(c, 8));
        if c > peak {
            peak = c;
        }
    }
    ratio(hi(peak), energy)
}

/// Magnitude spectrum: `|Re| + |Im|` per bin.
pub fn magnitude(s: &Spectrum) -> [i16; BINS] {
    let mut mag = [0i16; BINS];
    for (m, (&real, &imaginary)) in mag.iter_mut().zip(s.re.iter().zip(s.im.iter())) {
        let re = ((real as i64) << 16).abs();
        let im = ((imaginary as i64) << 16).abs();
        *m = hi(sat(acc(im + ((hi(re) as i64) << 16))));
    }
    mag
}
