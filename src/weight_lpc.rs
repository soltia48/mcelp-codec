//! Bandwidth expansion of the prediction filter.
//!
//! Multiplying the coefficients by successive powers of a factor below one
//! pulls the filter's poles in towards the origin, which broadens its
//! resonances.  The encoder uses two different factors to build the numerator
//! and denominator of its perceptual weighting filter.

use crate::fixed::{acc, hi, mulr};
use crate::lpc::ORDER;

/// Expand `a` by successive powers of `gamma`.
pub fn expand(a: &[i16; ORDER + 1], gamma: i16) -> [i16; ORDER + 1] {
    let mut out = [0i16; ORDER + 1];
    out[0] = a[0];
    let mut power = gamma;
    for i in 1..=ORDER {
        out[i] = hi(mulr(power, a[i]));
        if i < ORDER {
            power = hi(mulr(power, gamma));
        }
    }
    out
}

pub use crate::lsp::midpoint;

use crate::fixed::{low, mul, shift};

/// Knees of the log-area approximation, in the units the magnitude is reduced
/// to, and the slope and offset of the three segments above them.
const KNEE: [i64; 3] = [1299, 1815, 1944];
const SLOPE: [i16; 3] = [4567, 11776, 27443];
const OFFSET: [i64; 3] = [0x0031_eb85, 0x00f9_9999, 0x02ca_3d71];

/// Approximate the log area ratio of one reflection coefficient.
///
/// Below the first knee the coefficient is already close enough to its own
/// logarithm to be used as it stands; above it three straight segments of
/// increasing slope stand in for the curve, which is what keeps the comparisons
/// further down meaningful as the coefficient approaches unity.
pub fn log_area(coefficient: i16) -> i16 {
    let magnitude = shift(acc((coefficient as i64) << 16).abs(), -4);
    let reduced = hi(magnitude) as i64;
    let mut value = reduced;
    for segment in 0..3 {
        if reduced <= KNEE[segment] {
            break;
        }
        if segment == 2 || reduced <= KNEE[segment + 1] {
            let scaled = acc(mul(low(reduced), SLOPE[segment]) - (OFFSET[segment] << 1));
            value = hi(shift(scaled, 4)) as i64;
            break;
        }
    }
    let value = value as i16;
    if coefficient < 0 { -value } else { value }
}

/// Thresholds the two log areas are tested against, leaving and entering the
/// periodic state respectively.
const LEAVE: [i64; 2] = [-3564, 1331];
const ENTER: [i64; 2] = [-3113, 881];

/// The two weighting factors, chosen by the state the decision settles in.
const SHARP_NUMERATOR: i16 = 32113;
const FLAT_NUMERATOR: i16 = 30802;
const FLAT_DENOMINATOR: i16 = 19661;

/// Constants of the denominator's straight line in the sharp state, and the
/// range it is held to.
const DENOMINATOR_BASE: i64 = 1024 << 16;
const DENOMINATOR_SLOPE: i16 = -24576;
const DENOMINATOR_MAX: i64 = 22938;
const DENOMINATOR_MIN: i64 = 13107;

/// The weighting filter's shape, decided once per half-frame.
#[derive(Clone, Copy, Default)]
pub struct Sharpness {
    /// The previous half-frame's log areas, for the smoothing.
    previous: [i16; 2],
    /// Set while the spectrum is being treated as flat.
    flat: i16,
}

impl Sharpness {
    /// Choose the numerator and denominator factors for the two subframes.
    ///
    /// `reflection` holds the first two reflection coefficients of the current
    /// half-frame; `lsf` holds the line spectral frequencies each subframe will
    /// be weighted with.  The first subframe is judged on the log areas
    /// smoothed against the previous half-frame and the second on this
    /// half-frame's alone, so that the filter follows a changing spectrum
    /// without chattering on a steady one.
    pub fn choose(
        &mut self,
        reflection: &[i16; 2],
        lsf: &[&[i16; ORDER]; 2],
    ) -> ([i16; 2], [i16; 2]) {
        let mut area = [[0i16; 2]; 2];
        for i in 0..2 {
            let current = log_area(reflection[i]);
            area[0][i] = hi(shift(
                acc(((current as i64) + (self.previous[i] as i64)) << 16),
                -1,
            ));
            area[1][i] = current;
            self.previous[i] = current;
        }

        let mut numerator = [0i16; 2];
        let mut denominator = [0i16; 2];
        for (subframe, pair) in area.iter().enumerate() {
            let sharp = if self.flat != 0 {
                let leaving =
                    acc((pair[0] as i64) - LEAVE[0]) < 0 && acc((pair[1] as i64) - LEAVE[1]) > 0;
                if leaving {
                    self.flat = 0;
                }
                leaving
            } else {
                let staying =
                    acc((pair[0] as i64) - ENTER[0]) <= 0 && acc((pair[1] as i64) - ENTER[1]) >= 0;
                if !staying {
                    self.flat = 1;
                }
                staying
            };
            if sharp {
                numerator[subframe] = SHARP_NUMERATOR;
                denominator[subframe] = low(spread(lsf[subframe]));
            } else {
                numerator[subframe] = FLAT_NUMERATOR;
                denominator[subframe] = FLAT_DENOMINATOR;
            }
        }
        (numerator, denominator)
    }
}

/// Denominator factor from how closely the line spectrum is bunched.
///
/// Two lines close together mean a sharp resonance, which the weighting filter
/// has to follow more closely, so the smallest gap in the set is turned into a
/// factor by a straight line and held inside a fixed range.
fn spread(lsf: &[i16; ORDER]) -> i64 {
    let mut smallest = i64::MAX;
    for pair in lsf.windows(2) {
        smallest = smallest.min(acc(((pair[1] as i64) - (pair[0] as i64)) << 16));
    }
    let line = acc(DENOMINATOR_BASE + mul(hi(smallest), DENOMINATOR_SLOPE));
    shift(line, -11).clamp(DENOMINATOR_MIN, DENOMINATOR_MAX)
}

/// Accessors used by the trace replay to line the state up with the reference.
impl Sharpness {
    pub fn restore(&mut self, previous: [i16; 2], flat: i16) {
        self.previous = previous;
        self.flat = flat;
    }
    pub fn saved(&self) -> ([i16; 2], i16) {
        (self.previous, self.flat)
    }
}
