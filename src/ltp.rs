//! Long-term (pitch) postfilter.
//!
//! Applied to the LPC residual before it is re-synthesised.  The reference
//! spreads it over five routines; the shape is:
//!
//! 1. rescale the residual history so the search has full precision;
//! 2. refine the transmitted pitch lag by +-1 sample, then over seven
//!    fractional phases in both directions, always maximising
//!    `correlation^2 / energy`;
//! 3. build the delayed signal a second time with a longer 16-tap fractional
//!    filter and keep whichever version correlates better;
//! 4. mix the delayed signal into the residual with a gain derived from the
//!    winning correlation.

use crate::decoder::RES_HISTORY;
use crate::fixed::{acc, exp, hi, low, mul, sat, shift};
use crate::pitch::SUBFRAME;

/// Samples of residual the search looks at: history plus current subframe.
pub const WINDOW: usize = RES_HISTORY + SUBFRAME;
/// Fractional phases tried on each side of the integer lag.
const PHASES: usize = 7;
/// Samples produced per interpolated phase (one more than a subframe so both
/// the leading and trailing window can be scored).
const INTERP_LEN: usize = SUBFRAME + 1;
/// Taps of the short interpolation filter used during the search.
const SEARCH_TAPS: usize = 4;
/// Taps of the long interpolation filter used once the phase is settled.
const REFINE_TAPS: usize = 16;

/// A normalised `correlation^2 / energy` score, kept as mantissa/exponent pairs
/// so candidates at different scales can be compared by cross multiplication.
#[derive(Clone, Copy, Default)]
struct Score {
    /// Numerator, `correlation^2`.
    numerator: i64,
    /// Denominator, the delayed signal's energy.
    denominator: i64,
    /// Correlation mantissa of the winner.
    correlation: i16,
    /// Fractional phase, 0 when the integer lag was best.
    phase: i16,
    /// Which of the two 80-sample windows the phase selects.
    direction: i16,
}

/// Outcome of the search.
pub struct Search {
    /// Integer lag the postfilter should use.
    pub lag: usize,
    /// Fractional phase, `0` meaning "no fractional filtering".
    pub phase: i16,
    /// Window selector that goes with [`Self::phase`].
    pub direction: i16,
    /// Correlation mantissa and exponent of the winner.
    pub correlation: i16,
    /// Energy mantissa of the winner.
    pub energy: i16,
    /// Exponents belonging to the two mantissas.
    pub corr_exp: i16,
    /// Exponent of the energy.
    pub energy_exp: i16,
    /// The interpolated signal that won; all zero when the winning phase was
    /// the integer one, which never reads it.
    pub interpolated: [i16; SUBFRAME],
}

/// Energy of `n` samples ending at `end` (exclusive) of `sig`.
fn energy(sig: &[i16]) -> i64 {
    let mut e = 0i64;
    for &v in sig {
        e = acc(e + (v as i64) * (v as i64) * 2);
    }
    e
}

/// Cross correlation of two equal-length blocks.
fn correlate(a: &[i16], b: &[i16]) -> i64 {
    let mut c = 0i64;
    for (&x, &y) in a.iter().zip(b) {
        c = acc(c + mul(x, y));
    }
    c
}

/// Search for the best pitch delay to run the postfilter at.
///
/// `res` is the rescaled residual window; the subframe being filtered starts at
/// `RES_HISTORY`.  Returns `None` when the residual is too weak or no candidate
/// correlates positively, which switches the postfilter off.
pub fn search(res: &[i16; WINDOW], lag0: i16) -> Option<Search> {
    let here = RES_HISTORY;
    let current = &res[here..here + SUBFRAME];

    let e0 = energy(current);
    if acc(e0 - (1i64 << 16)) < 0 {
        return None;
    }
    let e0_exp = exp(e0) as i16;
    let e0_mant = hi(shift(e0, e0_exp as i32));

    // Integer refinement over lag-1, lag, lag+1.
    let mut best_corr = -32768i64 << 16;
    let mut best_k = 0usize;
    for k in 0..3 {
        let delay = (lag0 as usize) - 1 + k;
        let c = correlate(current, &res[here - delay..here - delay + SUBFRAME]);
        if c >= best_corr {
            best_corr = c;
            best_k = k;
        }
        best_corr = sat(best_corr);
    }
    if best_corr <= 0 {
        return None;
    }
    let lag = (lag0 as usize) - 1 + best_k;

    let delayed = &res[here - lag..here - lag + SUBFRAME];
    let e1 = energy(delayed);
    if acc(e1 - (1i64 << 16)) < 0 {
        return None;
    }
    let e1 = sat(e1);

    // Seven interpolated versions of the delayed segment, and the energy of the
    // two 80-sample windows each of them offers.
    let mut interp = [0i16; PHASES * INTERP_LEN];
    let mut window_energy = [[0i16; PHASES]; 2];
    let mut peak = e1;
    for p in 0..PHASES {
        let h = &crate::tables::SEARCH_FILTER[SEARCH_TAPS * p..SEARCH_TAPS * p + SEARCH_TAPS];
        let base = here - lag + 1;
        for n in 0..INTERP_LEN {
            let mut a = 0i64;
            for (k, &c) in h.iter().enumerate() {
                a = acc(a + mul(c, res[base + n - k]));
            }
            interp[p * INTERP_LEN + n] = hi(acc(a + 0x8000));
        }
        let block = &interp[p * INTERP_LEN..(p + 1) * INTERP_LEN];
        let lead = sat(energy(&block[..SUBFRAME]));
        let trail = sat(energy(&block[1..]));
        window_energy[0][p] = hi(lead);
        window_energy[1][p] = hi(trail);
        // Only the stored 16-bit mantissa takes part in the running maximum.
        let edge = ((if (block[0] as i64).abs() > (block[SUBFRAME] as i64).abs() {
            window_energy[0][p]
        } else {
            window_energy[1][p]
        }) as i64)
            << 16;
        if edge > peak {
            peak = edge;
        }
    }
    if acc(peak - (1i64 << 16)) < 0 {
        return None;
    }

    // Everything is now compared at a common scale.
    let peak_exp = exp(peak) as i16;
    let corr_shift = peak_exp.min(e0_exp);
    let mut best = Score {
        numerator: {
            let m = hi(shift(best_corr, corr_shift as i32)) as i64;
            acc(m * m * 2)
        },
        denominator: shift(e1, peak_exp as i32),
        correlation: hi(shift(best_corr, corr_shift as i32)),
        phase: PHASES as i16,
        direction: 1,
    };

    for p in 0..PHASES {
        let block = &interp[p * INTERP_LEN..(p + 1) * INTERP_LEN];
        for (direction, energies) in window_energy.iter().enumerate() {
            let candidate = if direction == 0 {
                correlate(current, &block[..SUBFRAME])
            } else {
                correlate(current, &block[1..])
            };
            let c = shift(candidate, corr_shift as i32).max(0);
            let mantissa = hi(c);
            let numerator = acc((mantissa as i64) * (mantissa as i64) * 2);
            let den_mant = energies[p];
            if better_ratio(
                numerator,
                best.denominator,
                best.numerator,
                den_mant,
                peak_exp,
            ) {
                best = Score {
                    numerator,
                    denominator: shift((den_mant as i64) << 16, peak_exp as i32),
                    correlation: mantissa,
                    phase: (PHASES - 1 - p) as i16,
                    direction: direction as i16,
                };
            }
        }
    }

    if best.correlation == 0 || acc(best.denominator - (1i64 << 16)) <= 0 {
        return None;
    }

    // Final gate: the winner must beat the plain integer lag by enough.
    let scaled = shift(
        shift(acc(crate::fixed::mul32x16(best.denominator, e0_mant)), -1),
        (2 * corr_shift - peak_exp - e0_exp) as i32,
    );
    if acc(scaled - best.numerator) > 0 {
        return None;
    }

    let phase = (PHASES as i16) - best.phase;
    let mut interpolated = [0i16; SUBFRAME];
    if phase != 0 {
        let p = (phase - 1) as usize;
        let off = best.direction as usize;
        interpolated.copy_from_slice(&interp[p * INTERP_LEN + off..][..SUBFRAME]);
    }

    Some(Search {
        lag: lag + 1 - best.direction as usize,
        phase,
        direction: best.direction,
        correlation: best.correlation,
        energy: hi(best.denominator),
        corr_exp: corr_shift,
        energy_exp: peak_exp,
        interpolated,
    })
}

/// Is `new_num / new_den` a better ratio than `best_num / best_den`?
///
/// Compared by cross multiplication so no division is needed.  `best_den` is
/// held pre-scaled by `peak_exp` while the challenger's denominator is still a
/// raw mantissa, which is what the shift undoes.
fn better_ratio(new_num: i64, best_den: i64, best_num: i64, new_den: i16, peak_exp: i16) -> bool {
    let lhs = shift(mul32(best_den, new_num), -(peak_exp as i32));
    let rhs = crate::fixed::mul32x16(best_num, new_den);
    acc(lhs - rhs) > 0
}

/// 32x32 fractional product, spelled the way the reference builds it out of one
/// one unsigned multiply and three multiply-accumulates.
fn mul32(x: i64, y: i64) -> i64 {
    let xh = hi(x) as i64;
    let xl = low(x) as u16 as i64;
    let yh = hi(y) as i64;
    let yl = low(y) as u16 as i64;
    let mut a = acc(xl * yl);
    a = shift(a, -16);
    a = acc(a + yl * xh * 2);
    a = acc(a + xl * yh * 2);
    a = shift(a, -16);
    acc(a + xh * yh * 2)
}

/// Rebuild the delayed signal with the long filter and report its correlation
/// and energy.
pub fn refine(res: &[i16; WINDOW], lag: usize, phase: i16, out: &mut [i16; SUBFRAME]) -> Stats {
    let here = RES_HISTORY;
    let at = REFINE_TAPS * ((phase - 1) as usize);
    let h = &crate::tables::REFINE_FILTER[at..at + REFINE_TAPS];
    let base = here + 8 - lag;
    for n in 0..SUBFRAME {
        let mut a = 0i64;
        for (k, &c) in h.iter().enumerate() {
            a = acc(a + mul(c, res[base + n - k]));
        }
        out[n] = hi(sat(acc(a + 0x8000)));
    }

    let c = correlate(out, &res[here..here + SUBFRAME]);
    let mut stats = Stats::default();
    if c >= 0 {
        let e = exp(c) as i16;
        stats.corr_exp = e;
        stats.correlation = hi(shift(c, e as i32));
    }
    let en = energy(out);
    let e = exp(en) as i16;
    stats.energy_exp = e;
    stats.energy = hi(shift(en, e as i32));
    stats
}

/// A correlation/energy pair with the exponents needed to compare it.
#[derive(Clone, Copy, Default)]
pub struct Stats {
    /// Correlation mantissa.
    pub correlation: i16,
    /// Energy mantissa.
    pub energy: i16,
    /// Exponent of [`Self::correlation`].
    pub corr_exp: i16,
    /// Exponent of [`Self::energy`].
    pub energy_exp: i16,
}

/// Decide whether the refined delayed signal beats the one the search picked
///.  Returns true when the refinement wins.
pub fn refinement_wins(new: &Stats, old: &Stats) -> bool {
    if new.energy == 0 {
        return false;
    }
    let num_new = acc((new.correlation as i64) * (new.correlation as i64) * 2);
    let num_old = acc((old.correlation as i64) * (old.correlation as i64) * 2);
    let challenger = crate::fixed::mul32x16(num_new, old.energy);
    let incumbent = crate::fixed::mul32x16(num_old, new.energy);

    // Line the two products up before comparing them.
    let diff = 2 * (new.corr_exp as i64) - 2 * (old.corr_exp as i64) + (old.energy_exp as i64)
        - (new.energy_exp as i64);
    let (challenger, incumbent) = if diff > 0 {
        (shift(challenger, -(diff as i32)), incumbent)
    } else {
        (challenger, shift(incumbent, diff as i32))
    };
    acc(challenger - incumbent) > 0
}

/// Mix the delayed signal into the residual.
///
/// `gain` is the weight given to the original residual; the remainder goes to
/// the pitch-delayed copy.
pub fn mix(gain: i16, residual: &[i16], delayed: &[i16], out: &mut [i16; SUBFRAME]) {
    let complement = hi(acc((32768i64 << 16) - ((gain as i64) << 16)));
    for n in 0..SUBFRAME {
        let a = acc(mul(gain, residual[n]) + mul(complement, delayed[n]) + 0x8000);
        out[n] = hi(a);
    }
}

/// Weight given to the *original* residual when mixing.
///
/// Two thirds unless the delayed signal's energy dominates its correlation with
/// the residual, in which case the weight follows `E / (E + correlation)`.
pub fn mix_gain(stats: &Stats) -> i16 {
    /// Default weight, 2/3 in Q15.
    const DEFAULT: i16 = 21845;

    let shift_amount = (stats.corr_exp - stats.energy_exp) as i32;
    let scaled_energy = shift((stats.energy as i64) << 16, shift_amount);
    if acc(((stats.correlation as i64) << 16) - scaled_energy) > 0 {
        return DEFAULT;
    }
    let denominator = acc(scaled_energy + (stats.correlation as i64) * 16384 * 2);
    hi(shift(
        crate::fixed::divide(denominator, stats.energy),
        shift_amount,
    ))
}

/// Run the long-term postfilter over one subframe.
///
/// `residual` is the unscaled residual window; the current subframe occupies
/// its last [`SUBFRAME`] samples.  Returns the pitch lag the filter settled on,
/// which the excitation builder of the next half-frame reads as a voicing hint.
pub fn filter(residual: &[i16; WINDOW], lag0: i16, out: &mut [i16; SUBFRAME]) -> i16 {
    let here = RES_HISTORY;

    // Rescale the whole window so the correlations use the full word.
    let mut peak = 0i64;
    for &v in residual.iter() {
        let m = ((v as i64) << 16).abs();
        if m > peak {
            peak = m;
        }
    }
    peak = sat(peak);
    let headroom = if peak == 0 { 12 } else { exp(peak) - 3 };

    let mut scaled = [0i16; WINDOW];
    for (s, &v) in scaled.iter_mut().zip(residual.iter()) {
        *s = hi(shift((v as i64) << 16, headroom));
    }

    let found = match search(&scaled, lag0) {
        None => {
            out.copy_from_slice(&residual[here..]);
            return 0;
        }
        Some(s) => s,
    };

    // Optionally rebuild the delayed signal with the longer filter and keep
    // whichever of the two correlates better with the residual.
    let mut stats = Stats {
        correlation: found.correlation,
        energy: found.energy,
        corr_exp: found.corr_exp,
        energy_exp: found.energy_exp,
    };
    let mut delayed = [0i16; SUBFRAME];
    if found.phase == 0 {
        // No fractional filtering: the delayed signal is read straight out of
        // the unscaled residual, so no rescaling is needed either.
        delayed.copy_from_slice(&residual[here - found.lag..here - found.lag + SUBFRAME]);
    } else {
        delayed = found.interpolated;
        let mut refined = [0i16; SUBFRAME];
        let refined_stats = refine(&scaled, found.lag, found.phase, &mut refined);
        if refinement_wins(&refined_stats, &stats) {
            stats = refined_stats;
            delayed = refined;
        }
        // Undo the search scaling before mixing.
        for v in delayed.iter_mut() {
            *v = hi(sat(shift((*v as i64) << 16, -headroom)));
        }
    }

    let gain = mix_gain(&stats);
    mix(gain, &residual[here..], &delayed, out);
    found.lag as i16
}
