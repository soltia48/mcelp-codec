//! Critical-band analysis of the magnitude spectrum.
//!
//! The 129 magnitude bins are grouped into twenty bands of increasing width,
//! averaged, and converted to decibels.  The band levels are what the rest of
//! the encoder's spectral machinery works on.

use crate::analysis::BINS;
use crate::fixed::{acc, exp, hi, low, mul32x16, sat, shift, square32};

/// Number of critical bands.
pub const BANDS: usize = 20;

/// Fractional part table of the base-2 logarithm.
/// `20*log10(2)` in Q12: converts a log2 value to decibels.  Shared by every
/// band-domain level in the encoder.
pub const DB_PER_OCTAVE: i16 = 24660;
/// Floor applied to every band level.
const LEVEL_FLOOR: i16 = 2560;
/// Offset subtracted from the log before scaling.
const LOG_OFFSET: i16 = 4;

/// Base-2 logarithm of a positive accumulator.
///
/// Returns `log2(x)` with the integer part in the high half and a
/// table-interpolated fraction below it; zero for a non-positive input.
pub fn log2(value: i64) -> i64 {
    if value <= 0 {
        return 0;
    }
    let e = exp(value);
    let integer = 30 - e;
    let normalised = shift(value, e);

    let shifted = shift(normalised, -9);
    let index = low(shift(shifted, -16)) as usize;
    let frac = (shift(shifted, -1) & 32767) as i16;

    let base = (crate::tables::LOG2_FRACTION[index] as i64) << 16;
    let slope = acc(base - ((crate::tables::LOG2_FRACTION[index + 1] as i64) << 16));
    let interpolated = shift(acc(base - (hi(slope) as i64) * (frac as i64) * 2), -16);
    acc(((integer as i64) << 16) + shift(interpolated, 1))
}

/// Mean level of each critical band, in decibels.
///
/// `spectrum` is the magnitude spectrum after it has been scaled back by the
/// analysis headroom.
pub fn levels(spectrum: &[i16; BINS]) -> [i64; BANDS] {
    let mut out = [0i64; BANDS];
    let mut bin = 0usize;
    for (band, level_out) in out.iter_mut().enumerate() {
        let first = crate::tables::BAND_EDGES[2 * band] as usize;
        let last = crate::tables::BAND_EDGES[2 * band + 1] as usize;
        let _ = first;

        let mut sum = 0i64;
        for _ in 0..=(last.wrapping_sub(bin)) {
            sum = acc(sum + ((spectrum[bin] as i64) << 16));
            bin += 1;
        }
        let sum = sat(sum);

        let mean = mul32x16(sum, crate::tables::BAND_INV_WIDTH_HALVED[band]);
        let log = log2(shift(mean, 1));
        let level = acc(shift(log, -1) - ((LOG_OFFSET as i64) << 16));
        let db = sat(shift(mul32x16(level, DB_PER_OCTAVE), 12));
        *level_out = db.max((LEVEL_FLOOR as i64) << 16);
    }
    out
}

/// Offset subtracted from the log of a tracked band energy.
const TRACKED_OFFSET: i16 = 6;
/// Weight each band contributes to the summary level, Q15 (halved because the
/// multiply doubles).
const SUMMARY_WEIGHT: i16 = 1638;
/// Upper clip on the level difference.
const DIFF_MAX: i16 = 20480;
/// Lower clip on the level difference.
const DIFF_MIN: i16 = -5120;

/// Convert one band energy to decibels.
///
/// The same shaping as [`levels`], but applied to an energy that is already
/// tracked at full 32-bit precision rather than summed from the spectrum.
pub fn tracked_level(energy: i64) -> i64 {
    let log = log2(energy);
    let level = acc(shift(log, -1) - ((TRACKED_OFFSET as i64) << 16));
    let db = sat(shift(mul32x16(level, DB_PER_OCTAVE), 12));
    db.max((LEVEL_FLOOR as i64) << 16)
}

/// Convert a whole set of tracked band energies.
pub fn tracked_levels(energies: &[i64; BANDS]) -> [i64; BANDS] {
    let mut out = [0i64; BANDS];
    for (o, &e) in out.iter_mut().zip(energies.iter()) {
        *o = tracked_level(e);
    }
    out
}

/// Weighted sum of a set of band levels.
pub fn summary(levels: &[i64; BANDS]) -> i64 {
    let mut total = 0i64;
    for &l in levels.iter() {
        total = acc(total + mul32x16(l, SUMMARY_WEIGHT));
    }
    total
}

/// Difference between the current and the tracked summary level, clipped.
pub fn summary_difference(current: i64, tracked: i64) -> i64 {
    let diff = acc(current - tracked);
    diff.min((DIFF_MAX as i64) << 16)
        .max((DIFF_MIN as i64) << 16)
}

/// Upper clip on a per-band excess level.
const EXCESS_MAX: i16 = 20480;

/// Per-band excess of the current level over a tracked reference, clipped to
/// `0..=20480`.
///
/// The two tracked sets give two different references, so this is run twice.
pub fn excess(levels: &[i64; BANDS], tracked: &[i64; BANDS]) -> [i64; BANDS] {
    let mut out = [0i64; BANDS];
    for i in 0..BANDS {
        let d = acc(levels[i] - tracked[i]);
        out[i] = d.min((EXCESS_MAX as i64) << 16).max(0);
    }
    out
}

/// Bands taking part in the SNR statistics (band 0 is left out).
const STAT_BANDS: usize = 18;
/// Weight applied when summing the per-band SNR, Q15.
const STAT_MEAN_WEIGHT: i16 = 8192;
/// `1/18` in Q15, used to turn the sums into a variance.
const STAT_INV_BANDS: i16 = 1820;

/// Weighted sum and variance of the per-band SNR.
///
/// Band 0 is excluded.  These are the spectral-flatness style features the
/// encoder's voice-activity logic works from.
pub fn snr_statistics(snr: &[i64; BANDS]) -> (i64, i64) {
    let mut sum = 0i64;
    let mut sum_sq = 0i64;
    for &v in snr[1..1 + STAT_BANDS].iter() {
        sum = acc(sum + mul32x16(v, STAT_MEAN_WEIGHT));
        sum_sq = acc(sum_sq + square32(v));
    }
    let scaled = shift(sum_sq, -4);
    let mean_sq = mul32x16(square32(sum), STAT_INV_BANDS);
    let variance = shift(mul32x16(acc(scaled - mean_sq), STAT_INV_BANDS), 4);
    (sum, variance)
}
