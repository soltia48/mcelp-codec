//! Linear prediction analysis of the shaping signal.
//!
//! The autocorrelation is taken over a 400-sample window that spans both
//! half-frames plus some history, and is then lag-windowed so the eventual
//! prediction filter stays well conditioned.

use crate::fixed::{acc, exp, hi, norm_shift, sat, shift};
use crate::tables::LAG_WINDOW;

/// Samples the analysis window covers.
pub const ANALYSIS: usize = 400;
/// Order of the prediction filter, and so the number of lags kept.
pub const ORDER: usize = 10;

/// Autocorrelation of one analysis window.
///
/// Returns the eleven correlations, all normalised by the same shift so their
/// ratios survive; the energy is nudged up by one so a silent window still has
/// a defined logarithm.
pub fn autocorrelate(signal: &[i16; ANALYSIS]) -> [i64; ORDER + 1] {
    let w = &crate::tables::ANALYSIS_WINDOW[..ANALYSIS];
    let mut windowed = [0i16; ANALYSIS];
    for (o, (&x, &c)) in windowed.iter_mut().zip(signal.iter().zip(w.iter())) {
        *o = hi(acc((x as i64) * (c as i64) * 2));
    }

    let mut energy = 1i64;
    for &v in windowed.iter() {
        energy = acc(energy + (v as i64) * (v as i64) * 2);
    }
    let shift_out = exp(energy);

    let mut out = [0i64; ORDER + 1];
    out[0] = shift(energy, shift_out);
    for lag in 1..=ORDER {
        let mut sum = 0i64;
        for i in 0..ANALYSIS - lag {
            sum = acc(sum + (windowed[i] as i64) * (windowed[i + lag] as i64) * 2);
        }
        out[lag] = shift(sum, shift_out);
    }
    out
}

/// Apply the lag window to the correlations.
///
/// Only the lags are touched; the zeroth correlation is left as it is.
pub fn lag_window(correlation: &mut [i64; ORDER + 1]) {
    for (lag, c) in correlation.iter_mut().enumerate().skip(1) {
        *c = mul32(*c, pair(2 * (lag - 1)));
    }
}

/// Read a 32-bit constant stored as a high/low word pair.
fn pair(at: usize) -> i64 {
    (((LAG_WINDOW[at] as i64) << 16) | (LAG_WINDOW[at + 1] as u16 as i64)) as i32 as i64
}

/// Full 32x32 product, keeping the upper half and clipped, which is the form
/// every caller in this module wants.
#[inline]
fn mul32(x: i64, y: i64) -> i64 {
    sat(crate::fixed::mul32(x, y))
}

/// Reciprocal of a 32-bit value, to 31 bits.
///
/// A plain restoring division of `0x3fffffff` by the normalised magnitude,
/// one bit at a time.  Returns the quotient and the shift that has to be
/// applied afterwards.
fn reciprocal(value: i64) -> (i64, i16) {
    let e = exp(value);
    let normalised = shift(value, e);
    let negative = normalised < 0;
    let magnitude = sat(if normalised < 0 {
        -normalised
    } else {
        normalised
    });

    let mut quotient = 0i64;
    let mut remainder = 0x3fffffffi64;
    for _ in 0..31 {
        remainder = shift(remainder, 1);
        quotient = shift(quotient, 1);
        remainder = acc(remainder - magnitude);
        if remainder >= 0 {
            quotient |= 1;
        } else {
            remainder = acc(remainder + magnitude);
        }
    }
    (if negative { -quotient } else { quotient }, (e + 1) as i16)
}

/// `numerator / denominator` for two 32-bit values.
pub fn divide(numerator: i64, denominator: i64) -> i64 {
    let sign = if denominator >= 0 { 1 } else { -1 };
    let magnitude = if denominator < 0 {
        -denominator
    } else {
        denominator
    };
    let e = exp(magnitude);
    let (recip, adjust) = reciprocal(shift(magnitude, e));
    let product = mul32(numerator, recip);
    let signed = if sign < 0 { acc(-product) } else { product };
    shift(signed, (e as i16 + adjust) as i32)
}

/// Where the recursion starts: one, in the format the coefficients are kept in.
const UNITY: i64 = 16384 << 14;

/// Levinson-Durbin recursion over the lag-windowed correlations.
///
/// Returns the eleven prediction coefficients, rounded to sixteen bits, and the
/// ten reflection coefficients the recursion passed through.  Everything is
/// carried at 32 bits; the residual energy is renormalised at every step and
/// its exponent travels with it.
pub fn levinson(r: &[i64; ORDER + 1]) -> ([i16; ORDER + 1], [i16; ORDER]) {
    let mut a = [0i64; ORDER + 1];
    let mut reflection = [0i16; ORDER];

    let negated = acc(-divide(r[1], r[0]));
    reflection[0] = hi(negated);
    a[0] = UNITY;
    a[1] = shift(negated, -3);

    let mut energy = acc(r[0] + mul32(r[1], negated));
    let mut exponent = exp(energy);
    energy = shift(energy, exponent);

    for order in 2..=ORDER {
        let mut b = 0i64;
        for j in 0..order {
            b = acc(b + mul32(a[j], r[order - j]));
        }
        // The normalising shift comes from a signed six-bit field, so a residual
        // exponent that has grown past thirty-one wraps rather than shifting on.
        let k = shift(acc(-divide(b, energy)), norm_shift(exponent));
        let negated = shift(k, 3);
        reflection[order - 1] = hi(negated);

        let previous = a;
        for j in 1..order {
            a[j] = acc(previous[j] + mul32(negated, previous[order - j]));
        }
        // The new coefficient is kept three bits down, like every other one;
        // the reference writes it at full scale first and shifts it afterwards.
        a[order] = shift(negated, -3);

        // The correction is brought up by the residual's own exponent before it
        // is added, and then by the three bits the coefficients are kept down.
        energy = acc(energy + shift(mul32(b, negated), norm_shift(exponent) + 3));
        let e = exp(energy);
        energy = shift(energy, e);
        exponent += e;
    }

    let mut out = [0i16; ORDER + 1];
    for (o, &v) in out.iter_mut().zip(a.iter()) {
        *o = hi(acc(v + (1 << 15)));
    }
    (out, reflection)
}

/// Half of the prediction order: the number of terms in each half-polynomial.
const HALF_ORDER: usize = ORDER / 2;

/// Split the prediction filter into its symmetric and antisymmetric halves.
///
/// These are the polynomials whose roots interleave on the unit circle and
/// become the line spectral pair; each is built by folding the coefficients
/// from both ends inwards while subtracting the previous term back out.
pub fn split_polynomials(a: &[i16; ORDER + 1]) -> ([i16; HALF_ORDER + 1], [i16; HALF_ORDER + 1]) {
    let mut sum = [0i16; HALF_ORDER + 1];
    let mut difference = [0i16; HALF_ORDER + 1];
    sum[0] = hi(16384i64 << 13);
    difference[0] = sum[0];

    for k in 0..HALF_ORDER {
        let (low, high) = (a[1 + k] as i64, a[ORDER - k] as i64);
        let s = sat(shift(
            acc(((low + high) << 16) - ((sum[k] as i64) << 17)),
            -1,
        ));
        sum[k + 1] = hi(s);
        let d = sat(shift(
            acc(((low - high) << 16) + ((difference[k] as i64) << 17)),
            -1,
        ));
        difference[k + 1] = hi(d);
    }
    (sum, difference)
}

/// Terms of each half-polynomial the grid search evaluates.
const CHEB_TERMS: usize = 5;

/// Evaluate one half-polynomial at a point on the unit circle.
///
/// A Chebyshev recurrence: each step folds in the next coefficient, twice the
/// point times the running term, and the term before that.  The running terms
/// are kept three bits down so the recurrence cannot overflow.
pub fn chebyshev(poly: &[i16; CHEB_TERMS], x: i16) -> i64 {
    let mut back = 16384i64 << 10;
    let mut term = acc(((poly[0] as i64) << 16) + ((x as i64) << 13));
    for &c in poly[1..CHEB_TERMS - 1].iter() {
        let state = shift(term, -3);
        term = acc(((c as i64) << 16) - shift(back, 3));
        term = acc(term + shift(mul32(state, (x as i64) << 16), 4));
        back = state;
    }
    let state = shift(term, -3);
    let mut b = shift(mul32(state, (x as i64) << 16), 3);
    b = acc(b + ((poly[CHEB_TERMS - 1] as i64) << 15));
    b = acc(b - shift(back, 3));
    sat(shift(b, 3))
}

/// Points on the unit circle the root search steps through.
/// How many of them.
const GRID_POINTS: usize = 60;
/// Bisection steps taken once a sign change has been bracketed.
const BISECTIONS: usize = 4;

/// Find the line spectral pair.
///
/// The two half-polynomials have interleaving roots, so the search alternates
/// between them: every time a sign change is bracketed and refined, the root
/// becomes the next starting point and the other polynomial takes over.
pub fn line_spectrum(
    sum: &[i16; HALF_ORDER + 1],
    difference: &[i16; HALF_ORDER + 1],
) -> [i16; ORDER] {
    let terms = |from: &[i16; HALF_ORDER + 1]| -> [i16; CHEB_TERMS] {
        let mut a = [0i16; CHEB_TERMS];
        a.copy_from_slice(&from[1..]);
        a
    };
    let (first, second) = (terms(sum), terms(difference));

    let mut out = [0i16; ORDER];
    let mut found = 0usize;
    let mut poly = &first;
    let mut index = 0usize;

    let mut x = crate::tables::ROOT_GRID[index];
    index += 1;
    let mut value = hi(chebyshev(poly, x));

    // The grid budget is only spent on points that turn out not to bracket a
    // root; the reference winds its counter back whenever it finds one.
    let mut remaining = GRID_POINTS;
    while remaining > 0 {
        let previous = x;
        let before = value;
        x = crate::tables::ROOT_GRID[index];
        index += 1;
        value = hi(chebyshev(poly, x));
        if acc((before as i64) * (value as i64) * 2) > 0 {
            remaining -= 1;
            continue;
        }
        index -= 1;

        let (mut low, mut high) = (x, previous);
        let (mut low_value, mut high_value) = (value, before);
        for _ in 0..BISECTIONS {
            let mid = hi(shift(acc(((low as i64) << 16) + ((high as i64) << 16)), -1));
            let v = hi(chebyshev(poly, mid));
            if acc((low_value as i64) * (v as i64) * 2) > 0 {
                low = mid;
                low_value = v;
            } else {
                high = mid;
                high_value = v;
            }
        }

        // Interpolate between the two ends of the bracket.
        let span = hi(acc(((high as i64) << 16) - ((low as i64) << 16)));
        let drop = acc(((high_value as i64) << 16) - ((low_value as i64) << 16));
        let slope = shift(crate::fixed::divide(drop, span), -5);
        let root = hi(acc(-acc(shift(
            acc((low_value as i64) * (hi(slope) as i64) * 2),
            5,
        ) - ((low as i64) << 16))));

        out[found] = root;
        found += 1;
        if found == ORDER {
            break;
        }

        poly = if std::ptr::eq(poly, &first) {
            &second
        } else {
            &first
        };
        x = root;
        value = hi(chebyshev(poly, x));
    }
    out
}
