//! Weighting of the line spectral frequencies before quantisation.
//!
//! A frequency that sits close to its neighbours describes a sharp formant and
//! matters more, so it is given a larger weight; one with room around it is
//! given the flat minimum.  The two middle weights are then trimmed and the
//! whole set is normalised up to full scale.

use crate::fixed::{acc, exp, hi, low, sat, shift};
use crate::lpc::ORDER;

/// Flat weight every frequency gets before the closeness bonus.
const FLOOR: i16 = 16384;
/// How strongly closeness is rewarded.
const SLOPE: i16 = 20480;
/// Spacing at which the bonus runs out, for the first and last frequency.
const FIRST_LIMIT: i16 = 9221;
const LAST_LIMIT: i16 = 15485;
/// Spacing at which it runs out in between.
const INNER_LIMIT: i16 = 8192;
/// The two middle weights are trimmed by this.
const TRIM: i16 = 19661;
/// Shift used when a set is entirely zero.
const EMPTY_SHIFT: i16 = 15;

/// Per-frequency weights for the quantiser's distance measure.
pub fn weights(lsf: &[i16; ORDER]) -> [i16; ORDER] {
    let mut w = [0i16; ORDER];
    w[0] = bonus(acc(((lsf[1] as i64) - (FIRST_LIMIT as i64)) << 16));
    for i in 1..ORDER - 1 {
        let spacing = shift(acc(((lsf[i + 1] as i64) - (lsf[i - 1] as i64)) << 16), 1);
        w[i] = bonus(shift(acc(spacing - ((INNER_LIMIT as i64) << 17)), -1));
    }
    w[ORDER - 1] = bonus(acc(((LAST_LIMIT as i64) - (lsf[ORDER - 2] as i64)) << 16));

    for v in w[4..6].iter_mut() {
        let trimmed = sat((TRIM as i64) * (*v as i64) * 2);
        *v = hi(sat(shift(trimmed, 1)));
    }

    // Bring the set up to full scale so the search keeps its resolution.
    let mut peak = 0i64;
    for &v in w.iter() {
        peak = peak.max(sat(((v as i64) << 16).abs()));
    }
    let up = if peak == 0 {
        EMPTY_SHIFT as i32
    } else {
        exp(peak)
    };
    for v in w.iter_mut() {
        *v = low(sat(shift(acc((*v as i64) << 16), up - 16)));
    }
    w
}

/// Turn a spacing into a weight: the closer the neighbours, the bigger it gets.
fn bonus(spacing: i64) -> i16 {
    if spacing > 0 {
        return hi(shift((FLOOR as i64) << 16, -3));
    }
    let square = sat(shift(
        sat((hi(spacing) as i64) * (hi(spacing) as i64) * 2),
        2,
    ));
    let scaled = sat((SLOPE as i64) * (hi(square) as i64) * 2);
    hi(shift(acc(((FLOOR as i64) << 16) + shift(scaled, 5)), -3))
}

/// Entries in the first-stage codebook.
const STAGE1_ENTRIES: usize = 128;

/// Pick the closest first-stage entry.
///
/// A plain exhaustive search on squared distance; ties go to the earlier entry
/// because the running minimum is only replaced on a strict improvement.
pub fn search_stage1(target: &[i16; ORDER]) -> i16 {
    let mut best = 0x7fffffffi64;
    let mut winner = 0i16;
    for entry in 0..STAGE1_ENTRIES {
        let base = entry * ORDER;
        let mut distance = 0i64;
        for (j, &t) in target.iter().enumerate() {
            let d = acc(((t as i64) - (crate::tables::LSF_CB1[base + j] as i64)) << 16);
            distance = sat(acc(distance + (hi(d) as i64) * (hi(d) as i64) * 2));
        }
        if acc(distance - best) < 0 {
            winner = entry as i16;
        }
        best = best.min(distance);
    }
    winner
}

/// Smallest gap the quantiser leaves between neighbouring values.
const MIN_GAP: i16 = 10;

/// Push neighbouring values apart until they are at least `MIN_GAP` apart.
///
/// Any violation is split evenly between the two, and only pairs that actually
/// violate the gap are touched — the conditional store is what makes that
/// cheap.
pub fn enforce_spacing(values: &mut [i16]) {
    for i in 0..values.len() - 1 {
        let overlap = shift(
            acc((((values[i] as i64) - (values[i + 1] as i64)) << 16) + ((MIN_GAP as i64) << 16)),
            -1,
        );
        if overlap > 0 {
            let lowered = acc(((values[i] as i64) << 16) - overlap);
            let raised = acc(overlap + ((values[i + 1] as i64) << 16));
            values[i] = hi(lowered);
            values[i + 1] = hi(raised);
        }
    }
}

/// Weighted distortion between a candidate and the target.
///
/// The difference in each dimension is first scaled by the predictor's own
/// gain, then weighted, and the two are multiplied together.
///
/// The accumulation multiplies the scaled difference by the weight and adds
/// the upper half of the product, which is where the two shifts below come
/// from.
pub fn distortion(
    target: &[i16; ORDER],
    candidate: &[i16; ORDER],
    scale: &[i16; ORDER],
    weight: &[i16; ORDER],
) -> i64 {
    let mut total = 0i64;
    for i in 0..ORDER {
        let difference = acc(((candidate[i] as i64) - (target[i] as i64)) << 16);
        let scaled = hi(acc((hi(difference) as i64) * (scale[i] as i64) * 2));
        let weighted = shift(acc((scaled as i64) * (weight[i] as i64) * 2), 3);
        total = acc(total + (scaled as i64) * (hi(weighted) as i64) * 2);
    }
    shift(total, 1)
}

/// Entries in the second-stage codebook, and the width of one half.
const STAGE2_ENTRIES: usize = 64;
const HALF: usize = ORDER / 2;

/// Pick the closest second-stage entry for one half.
///
/// `offset` is 0 for the lower half and 5 for the upper one; the two halves are
/// quantised independently, which is what makes the second stage two six-bit
/// fields rather than one twelve-bit one.
pub fn search_stage2(residual: &[i16; HALF], weight: &[i16; HALF], offset: usize) -> i16 {
    let mut best = 0x7fffffffi64;
    let mut winner = 0i16;
    for entry in 0..STAGE2_ENTRIES {
        let base = entry * ORDER + offset;
        let mut total = 0i64;
        for k in 0..HALF {
            let d = hi(acc(((residual[k] as i64)
                - (crate::tables::LSF_CB2[base + k] as i64))
                << 16));
            let scaled = acc((d as i64) * (weight[k] as i64) * 2);
            total = acc(total + (d as i64) * (hi(scaled) as i64) * 2);
        }
        if acc(total - best) < 0 {
            winner = entry as i16;
            best = total;
        }
    }
    winner
}

/// What the quantiser decided for one frame.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Quantised {
    /// Which moving-average predictor won.
    pub mode: usize,
    /// First-stage entry.
    pub stage1: i16,
    /// Second-stage entries for the lower and upper halves.
    pub lower: i16,
    pub upper: i16,
}

impl Quantised {
    /// The two transport fields the frame carries.
    pub fn fields(&self) -> (i16, i16) {
        (
            ((self.mode as i16) << 7) | self.stage1,
            (self.lower << 6) | self.upper,
        )
    }
}

/// Quantise one line spectral vector.
///
/// Both predictors are tried in full and the one whose reconstruction sits
/// closer to the residual wins; ties go to the first, which is what the
/// reference's `>=` comparison amounts to.
pub fn quantise(
    lsf: &[i16; ORDER],
    weight: &[i16; ORDER],
    history: &[[i16; ORDER]; 3],
) -> Quantised {
    let mut best = Quantised::default();
    let mut best_distortion = 0i64;

    for mode in 0..2 {
        let residual = crate::lsp::residual(lsf, history, mode);
        let stage1 = search_stage1(&residual);
        let first = stage1 as usize * ORDER;

        let mut chosen = [0i16; ORDER];
        let mut index = [0i16; 2];
        for (half, offset) in [(0usize, 0usize), (1, HALF)] {
            let mut difference = [0i16; HALF];
            let mut w = [0i16; HALF];
            for k in 0..HALF {
                difference[k] = hi(acc(((residual[offset + k] as i64)
                    - (crate::tables::LSF_CB1[first + offset + k] as i64))
                    << 16));
                w[k] = weight[offset + k];
            }
            index[half] = search_stage2(&difference, &w, offset);
            let second = index[half] as usize * ORDER + offset;
            for k in 0..HALF {
                chosen[offset + k] = hi(acc(((crate::tables::LSF_CB1[first + offset + k] as i64)
                    + (crate::tables::LSF_CB2[second + k] as i64))
                    << 16));
            }
            // The two gap passes overlap by one: the second starts at the last
            // value of the first, so the seam between the halves is covered.
            let span = if half == 0 { 0..HALF } else { HALF - 1..ORDER };
            enforce_spacing(&mut chosen[span]);
        }
        crate::lsp::rearrange(&mut chosen, 5);

        let gain = &crate::tables::LSF_MA_GAIN[mode * ORDER..mode * ORDER + ORDER];
        let d = distortion(&residual, &chosen, gain.try_into().unwrap(), weight);
        if mode == 0 || acc(d - best_distortion) < 0 {
            best_distortion = d;
            best = Quantised {
                mode,
                stage1,
                lower: index[0],
                upper: index[1],
            };
        }
    }
    best
}
