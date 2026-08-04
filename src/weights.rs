//! Per-band suppression weights.
//!
//! This is the part of the noise-suppression path that decides *how much* of
//! each critical band to take away: the suppression depth comes from how
//! speech-like the frame looked, is laid out over the bands with a fixed tilt,
//! smoothed against the previous frame and finally turned into the per-band
//! weight the shaping applies.  The routine's last loop, which smooths the
//! scale itself, lives in [`crate::shaping`] beside the signal it acts on.

use crate::bands::BANDS;
use crate::fixed::{acc, hi, low, norm, round, sat, scale32, shift};

/// Three sets of per-band suppression depths, in dB, ordered from the deepest
/// to the shallowest.  Which one is used depends on the frame's level sum.
/// Fixed per-band weight, Q15, falling off towards the top of the band; it is
/// carried around divided by thirty-two.
/// Depth used when the frame looks like nothing but noise.
const DEEPEST: i16 = -10240;
/// Level sums at which the next depth set takes over.  The comparisons are
/// cumulative, so these are differences rather than absolute thresholds.
const DEPTH_STEPS: [i16; 3] = [2048, 4096, 3072];

/// Opening step of the suppression weights: the per-band depth and weight.
///
/// `sum` is the weighted level sum the noise floor returned.  A frame that scored
/// low gets the deepest suppression across every band; as the sum climbs the
/// depths shrink towards zero, and the deepest set is skipped entirely.
pub fn depths(sum: i16) -> ([i64; BANDS], [i16; BANDS]) {
    let mut remaining = sum as i64 - DEPTH_STEPS[0] as i64;
    let mut depth = [0i64; BANDS];
    if remaining <= 0 {
        depth = [shift(acc((DEEPEST as i64) << 16), -2); BANDS];
    } else {
        remaining -= DEPTH_STEPS[1] as i64;
        let set = if remaining <= 0 {
            0
        } else {
            remaining -= DEPTH_STEPS[2] as i64;
            if remaining <= 0 { 1 } else { 2 }
        };
        for (band, d) in depth.iter_mut().enumerate() {
            let v = crate::tables::DEPTH_TABLES[BANDS * set + band] as i64;
            *d = shift(acc(v << 16), -2);
        }
    }

    let mut weight = [0i16; BANDS];
    for (band, v) in weight.iter_mut().enumerate() {
        *v = low(shift(
            acc((crate::tables::BAND_WEIGHT[band] as i64) << 16),
            -21,
        ));
    }
    (depth, weight)
}

/// Weight the summed noise energy is scaled by before the logarithm.
const SUM_WEIGHT: i16 = 26214;
/// Offset subtracted from the log, and the dB conversion.
const LOG_OFFSET: i16 = 6;
use crate::bands::DB_PER_OCTAVE;

/// Overall level of the tracked noise, in dB.
///
/// The twenty band energies are summed, scaled down a little and turned into
/// decibels the same way the band levels are, except that no floor is applied.
/// This single number is what decides how hard the suppression is allowed to
/// push.
pub fn overall_level(energy: &[i64; BANDS]) -> i16 {
    let mut sum = 0i64;
    for &e in energy.iter() {
        sum = sat(acc(sum + e));
    }
    let scaled = shift(sat((SUM_WEIGHT as i64) * (hi(sum) as i64) * 2), -4);
    let log = crate::bands::log2(scaled);
    let level = acc(shift(log, -1) - ((LOG_OFFSET as i64) << 16));
    hi(sat(shift(scale32(level, DB_PER_OCTAVE), 12)))
}

/// The same product with the final accumulation rounded.
fn mul32r(x: i64, y: i16) -> i64 {
    round(scale32(x, y))
}

/// Two to the power of a value split into integer and fraction.
///
/// The fraction indexes a 32-entry table which is then interpolated linearly,
/// and the integer part becomes a shift.  The result is rounded up by a half.
fn pow2(integer: i16, fraction: i16) -> i64 {
    let scaled = sat((fraction as i64) * 32 * 2);
    let index = shift(scaled, -16) as usize;
    let frac = low(shift(scaled, -1)) & 32767;
    let base = (crate::tables::POW2_CURVE[index] as i64) << 16;
    let slope = acc(base - ((crate::tables::POW2_CURVE[index + 1] as i64) << 16));
    let interpolated = sat(acc(base - (frac as i64) * (hi(slope) as i64) * 2));
    let shifted = norm(interpolated, integer - 29);
    shift(acc(shifted + 1), -1)
}

/// Split a magnitude into the integer and fraction `pow2` wants.
fn split(value: i64) -> (i16, i16) {
    let magnitude = sat(if value < 0 { -value } else { value });
    let integer = low(acc(shift(shift(magnitude, -16), -12) + 15));
    let fraction = low(shift(magnitude, -13)) & 32767;
    (integer, fraction)
}

/// Scale factor applied before every exponential in this block.
const EXP_SCALE: i16 = 5443;
/// Weight the second exponential is trimmed by.
const TRIM: i16 = 19661;
/// Weight the band levels are summed with, in two halves.
const HALF_SUM_WEIGHT: i16 = 3277;
/// Per-band mixing coefficients between the two exponentials.
/// Values written into the mix array when the frame is silent.
const IDLE_MIX: i64 = 1024 << 16;
/// Total the two mixing weights always add up to.
const MIX_UNITY: i16 = 16384;

/// Level above which the amount of suppression is derived from the score
/// rather than from the level sum.
const LOUD_NOISE: i16 = 15360;

/// Middle step: the per-band mix of two exponentials.
pub struct Mix;

impl Mix {
    /// Returns the twenty 32-bit mix values and the secondary amount that the
    /// rest of the routine works from.
    pub fn run(
        level: &[i16; BANDS],
        sum: i16,
        noise_level: i16,
        score: i16,
        active: bool,
    ) -> ([i64; BANDS], i16) {
        if !active {
            return ([IDLE_MIX; BANDS], 0);
        }

        // Both halves of the band levels, weighted.
        let mut low_sum = 0i64;
        let mut high_sum = 0i64;
        for band in 0..BANDS / 2 {
            low_sum = sat(acc(
                low_sum + (level[band] as i64) * (HALF_SUM_WEIGHT as i64) * 2
            ));
            high_sum =
                sat(acc(high_sum
                    + (level[BANDS / 2 + band] as i64)
                        * (HALF_SUM_WEIGHT as i64)
                        * 2));
        }

        // How lopsided the spectrum is between its two halves.
        let mut change = sat(shift(sat(abs(acc(low_sum - high_sum))), 2));

        let amount;
        if (noise_level as i64) - LOUD_NOISE as i64 >= 0 {
            let bias = if (score as i64) - 1 > 0 {
                3072i64 << 16
            } else {
                sat(shift(3072i64 << 16, 1))
            };
            let scaled = shift(sat(shift(acc(((sum as i64) << 16) - bias), 2)), -16);
            amount = low(scaled.min(12288));
            change = change.min(24576i64 << 16);
        } else {
            change = change.min(8192i64 << 16);
            let x = acc(sum as i64 - 2048);
            amount = if (sum as i64) - 12288 < 0 {
                if acc(x - 1024) > 0 {
                    4096
                } else {
                    low(sat(shift(x, 4)))
                }
            } else if acc(x - 12288) > 0 {
                8192
            } else {
                hi(sat(shift((low(x) as i64) << 16, 2)))
            };
        }

        let (i, f) = split(scale32(change, EXP_SCALE));
        let first = sat(shift(pow2(i, f), 12));
        let (i, f) = split(sat((amount as i64) * (EXP_SCALE as i64) * 2));
        let second = sat(shift(pow2(i, f), 12));

        let product = sat(shift(scale32(first, hi(mul32r(second, TRIM))), 4));
        // The product replaces whichever half of the pair the band levels say
        // is the weaker one.
        let (a, b) = if acc(high_sum - low_sum) < 0 {
            (second, product)
        } else {
            (product, second)
        };

        // The pair is mixed with weights that always add up to a half.
        let mut out = [0i64; BANDS];
        for (band, mixed_out) in out.iter_mut().enumerate() {
            let c = crate::tables::MIX_COEFFS[band];
            let weight = hi(acc(((MIX_UNITY as i64) << 16) - ((c as i64) << 16)));
            let mixed = sat(acc(scale32(a, weight) + scale32(b, c)));
            *mixed_out = shift(mixed, -1);
        }
        (out, amount)
    }
}

/// Absolute value in a saturating accumulator.
fn abs(v: i64) -> i64 {
    if v < 0 { -v } else { v }
}

/// Weights the depth is carried forward with, for speech and for silence.
const CARRY: [(i16, i16); 2] = [(31130, 1638), (26214, 6554)];
/// Offset and scale of the per-band level, as in the overall level.
const LEVEL_OFFSET: i16 = 6;

/// State the depth refinement carries between frames.
#[derive(Clone, Copy)]
pub struct Depths {
    /// Last frame's smoothed depth per band.
    previous: [i64; BANDS],
}

impl Default for Depths {
    /// A fresh encoder starts out assuming the deepest suppression everywhere.
    fn default() -> Self {
        Depths {
            previous: [((DEEPEST as i64) << 16) >> 2; BANDS],
        }
    }
}

impl Depths {
    /// Latter part, up to the final weighting.
    ///
    /// The depths are carried forward from the previous frame, then clamped so
    /// that no band is ever suppressed below its own tracked noise level, and
    /// finally folded into the mix values.  Everything ends up at or below
    /// zero: these are attenuations, never gains.
    pub fn refine(
        &mut self,
        depth: &mut [i64; BANDS],
        mix: &mut [i64; BANDS],
        energy: &[i64; BANDS],
        level: &[i16; BANDS],
        weight: &[i16; BANDS],
        active: bool,
    ) {
        let (keep, fresh) = CARRY[!active as usize];
        for (current, previous) in depth.iter_mut().zip(self.previous.iter_mut()) {
            let carried = sat(acc(scale32(*current, keep) + scale32(*previous, fresh)));
            *current = carried;
            *previous = carried;
        }

        for band in 0..BANDS {
            let log = crate::bands::log2(energy[band]);
            let scaled = acc(shift(log, -1) - ((LEVEL_OFFSET as i64) << 16));
            let ceiling = sat(shift(sat(-scale32(scaled, DB_PER_OCTAVE)), 12));
            if acc(ceiling - depth[band]) >= 0 {
                depth[band] = ceiling.min(0);
            }
        }

        for band in 0..BANDS {
            let floor = (weight[band] as i64) << 16;
            let excess = hi(acc(((level[band] as i64) << 16).max(floor) - floor));
            let folded = sat(shift(scale32(mix[band], excess), 3));
            mix[band] = sat(acc(folded + depth[band])).min(0);
        }
    }
}

/// Quadratic coefficients of the attenuation curve.
/// Per-band smoothing coefficients, shared with the scale smoother.
/// Total the smoothing weights add up to.
const SMOOTH_UNITY: i16 = 16384;
/// Level at which a band starts contributing at all.
const BAND_FLOOR: i16 = 6144;
/// Full-scale value the scale is measured down from.
const FULL_SCALE: i64 = 0x7fffffff;

/// A quadratic approximation of a fractional power of two.
///
/// The argument's high word is split into an integer and a fraction; the
/// fraction goes through a three-term polynomial and the integer becomes a
/// right shift.  A non-negative integer part means the result is already at
/// full scale and is returned as it stands.
fn curve(value: i64) -> i64 {
    let scaled = sat(32 * (hi(value) as i64) * 2);
    let integer = low(acc(shift(scaled, -16) + 1));
    let fraction = low(shift(acc((scaled & 0xffff) + ((-1i64) << 16)), -1));

    // Both the square and the half added to it saturate; at the very bottom of
    // the range the sum lands one count past full scale, and letting it wrap
    // would flip the sign of the leading term.
    let square = sat(acc(
        sat((fraction as i64) * (fraction as i64) * 2) + (1 << 15)
    ));
    let mut b = sat((crate::tables::CURVE[0] as i64) * (hi(square) as i64) * 2);
    b = sat(acc(
        b + (fraction as i64) * (crate::tables::CURVE[1] as i64) * 2
    ));
    let out = acc(b + ((crate::tables::CURVE[2] as i64) << 16));
    if integer >= 0 {
        out
    } else {
        norm(out, integer)
    }
}

/// Shape of the level-to-attenuation mapping, chosen per frame.
struct Slope {
    /// Constant part of the attenuation.
    offset: i16,
    /// Slope applied to the part of the level above the floor; zero disables
    /// the level-dependent term entirely.
    gain: i16,
    /// Bias added before that slope.
    bias: i16,
}

/// Final step: the per-band weights and scales.
#[derive(Clone, Copy)]
pub struct Shaping {
    /// Last frame's weights.
    previous: [i16; BANDS],
}

impl Default for Shaping {
    /// The weights start out at the shallowest of the floor gains.
    fn default() -> Self {
        Shaping {
            previous: [5194; BANDS],
        }
    }
}

impl Shaping {
    /// Returns the per-band weight and the per-band scale.
    pub fn finish(
        &mut self,
        mix: &mut [i64; BANDS],
        level: &[i16; BANDS],
        noise_level: i16,
        sum: i16,
        active: bool,
    ) -> ([i16; BANDS], [i16; BANDS]) {
        let table = BANDS * !active as usize;
        let mut weight = [0i16; BANDS];
        for band in 0..BANDS {
            let w = curve(nonzero(sat(shift(scale32(mix[band], EXP_SCALE), 2))));
            let c = crate::tables::BAND_SMOOTH[table + band];
            let mut b = sat((c as i64) * (hi(w) as i64) * 2);
            b = sat(acc(b
                + (self.previous[band] as i64)
                    * ((SMOOTH_UNITY - c) as i64)
                    * 2));
            let rounded = hi(acc(sat(shift(b, 1)) + (1 << 15)));
            weight[band] = rounded;
            self.previous[band] = rounded;
        }

        // How the band level maps onto attenuation depends on how loud the
        // frame is: a loud one gets a flat, deep floor, a quiet one a shallower
        // floor that the level is allowed to lift.
        let slope = if (noise_level as i64) - 15360 >= 0 {
            Slope {
                offset: -10240,
                gain: 0,
                bias: 0,
            }
        } else if acc(sum as i64 - 6144) < 0 {
            Slope {
                offset: -8192,
                gain: 0,
                bias: 0,
            }
        } else if acc(sum as i64 - 12288) < 0 {
            Slope {
                offset: -8192,
                gain: 29491,
                bias: 0,
            }
        } else {
            Slope {
                offset: -8192,
                gain: 29491,
                bias: 3072,
            }
        };

        let mut scale = [0i16; BANDS];
        for band in 0..BANDS {
            let above = acc(((level[band] as i64) << 16) - ((BAND_FLOOR as i64) << 16));
            let mut a = (slope.offset as i64) << 16;
            if above > 0 {
                let mut b = shift(above, -1);
                if slope.gain != 0 {
                    b = sat(shift(b, 1));
                    let lifted = acc(b + ((slope.bias as i64) << 16));
                    b = sat((slope.gain as i64) * (hi(lifted) as i64) * 2);
                }
                a = acc(b + ((slope.offset as i64) << 16)) & !0xffff;
            }
            a = a.min(0);
            a = acc(a - sat(shift(mix[band], 2))).min(0);
            a = shift(a, -2);
            mix[band] = ((hi(a) as i64) << 16) | (a & 0xffff);
            let w = curve(nonzero(sat(shift(scale32(mix[band], EXP_SCALE), 2))));
            scale[band] = hi(acc(FULL_SCALE - w));
        }
        (weight, scale)
    }
}

/// A zero argument to the curve is nudged to the smallest negative value; the
/// curve is only ever asked for attenuations.
fn nonzero(v: i64) -> i64 {
    if v == 0 { -1 } else { v }
}
