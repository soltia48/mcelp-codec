//! The twin five-pulse codebook.
//!
//! The last of the six fixed-codebook classes is not multi-pulse at all: it is
//! a pair of 128-entry codebooks, each entry five signed unit pulses, and the
//! innovation is one entry drawn from each.  It is searched on its own terms
//! and offered to [`crate::pulses`] as one more candidate for the subframe,
//! which keeps whichever scores better.
//!
//! The search shortlists eight shapes per codebook by how well each correlates
//! with the target, builds their filtered vectors, takes out the part the
//! adaptive codebook already accounts for, and then scores every one of the
//! sixty-four cross-codebook pairs.

use crate::SUBFRAME as SPAN;
use crate::fixed::{acc, hi, low, mul, sat, shift, trunc32};
use crate::pulses::{Search, takes_lead};

/// The adaptive codebook's filtered response.
///
/// The weighted impulse response is copied out doubled, and if the pitch lag is
/// shorter than a subframe the tail is filled in by the same periodic extension
/// the codebook vector itself gets.
pub fn response(impulse: &[i16], lag: i16, fraction: i16, gain: i16) -> [i16; SPAN] {
    let mut out = [0i16; SPAN];
    for (o, &v) in out.iter_mut().zip(impulse.iter()) {
        *o = hi(shift(acc((v as i64) << 16), 1));
    }
    if acc(SPAN as i64 - (lag as i64)) > 0 {
        crate::pitch::extend(&mut out, lag, fraction, gain);
    }
    out
}

/// Correlation of the two responses at every lag.
///
/// The reference walks the two arrays forwards and backwards in turn so that
/// neither pointer ever has to be reset, which is why the lags come out in
/// pairs; the result is the plain cross-correlation all the same.
pub fn lag_products(x: &[i16], y: &[i16]) -> [i16; SPAN] {
    let mut out = [0i16; SPAN];
    for (lag, o) in out.iter_mut().enumerate() {
        let mut total = 0i64;
        for k in lag..SPAN {
            total = acc(total + mul(x[k], y[k - lag]));
        }
        *o = hi(shift(total, -1));
    }
    out
}

/// Entries in a pulse-shape codebook, and pulses in each entry.
const SHAPES: usize = 128;
const PULSES: usize = 5;
/// Shapes kept from each codebook.
pub const SHORTLIST: usize = 8;

/// Pick the most promising pulse shapes from one codebook.
///
/// Each shape is five signed positions; scoring it against the response's lag
/// correlations is a matter of adding up the correlation at each position with
/// that position's sign.  The eight best-scoring shapes are kept in a list held
/// in order by insertion, and what comes back is their numbers, worst first.
pub fn shapes(table: &[i16], products: &[i16]) -> [i16; SHORTLIST] {
    let mut score = [0i16; SHORTLIST];
    let mut number = [127i16; SHORTLIST];

    for entry in 0..SHAPES {
        let mut total = 0i16;
        for p in 0..PULSES {
            let pulse = table[entry * PULSES + p];
            let term = products[pulse.unsigned_abs() as usize - 1] as i64;
            let sum = (total as i64) + if pulse > 0 { term } else { -term };
            total = hi(acc(sum << 16));
        }

        if total <= score[SHORTLIST - 1] {
            continue;
        }
        let mut at = 0usize;
        for j in (0..SHORTLIST - 1).rev() {
            if total < score[j] {
                at = j + 1;
                break;
            }
            score[j + 1] = score[j];
            number[j + 1] = number[j];
        }
        score[at] = total;
        number[at] = (SHAPES - 1 - entry) as i16;
    }

    let mut out = [0i16; SHORTLIST];
    for k in 0..SHORTLIST {
        out[SHORTLIST - 1 - k] = low(acc((SHAPES as i64 - 1) - (number[k] as i64)));
    }
    out
}

/// Quarter the response so that five of them cannot overflow the word.
pub fn quarter(response: &mut [i16; SPAN]) {
    for v in response.iter_mut() {
        *v = hi(shift(acc((*v as i64) << 16), -2));
    }
}

/// Lay the shortlisted shapes out as filtered vectors.
///
/// Each shape is five signed positions; a positive one adds the response
/// starting there and a negative one subtracts it.  Both codebooks are laid out
/// the same way, against a response [`quarter`] has already scaled.
pub fn shape_vectors(response: &[i16; SPAN], numbers: &[i16], table: &[i16]) -> Vec<i16> {
    let mut out = vec![0i16; SHORTLIST * SPAN];
    for (slot, &number) in numbers.iter().enumerate() {
        let base = slot * SPAN;
        for p in 0..PULSES {
            let pulse = table[(number as usize) * PULSES + p] as i64;
            let at = pulse.unsigned_abs() as usize - 1;
            for k in 0..SPAN - at {
                let carried = (out[base + at + k] as i64) << 16;
                let term = (response[k] as i64) << 16;
                out[base + at + k] = hi(acc(if pulse > 0 {
                    carried + term
                } else {
                    carried - term
                }));
            }
        }
    }
    out
}

/// Take the adaptive codebook's direction out of each shape vector.
///
/// The gain search that follows assumes the two codebooks are independent, so
/// whatever each shape has in common with the adaptive contribution is
/// subtracted off first.  A zero gain leaves the shapes alone.
pub fn orthogonalise(vectors: &mut [i16], contribution: &[i16], gain: i16, gain_shift: i16) {
    if acc(gain as i64) <= 0 {
        return;
    }
    for slot in 0..SHORTLIST {
        let base = slot * SPAN;
        let mut total = 0i64;
        for k in 0..SPAN {
            total = acc(total + mul(vectors[base + k], contribution[k]));
        }
        let share = crate::fixed::mul32x16(trunc32(sat(total)), gain);
        let scale = hi(shift(shift(share, gain_shift as i32), -2));
        for k in (0..SPAN).rev() {
            let term = shift(acc(mul(scale, contribution[k])), 2);
            vectors[base + k] = hi(acc(((vectors[base + k] as i64) << 16) - term));
        }
    }
}

/// Shape vectors the two codebooks contribute between them.
const CANDIDATES: usize = 2 * SHORTLIST;

/// Score every shape vector against the target.
///
/// Each vector gets its correlation with the target and its own energy; both
/// sets are then scaled up together by whatever their largest member allows, so
/// that the comparison that follows can work in single words.
pub fn measure_shapes(target: &[i16], vectors: &[i16]) -> (Vec<i64>, i16, Vec<i64>, i16) {
    let mut correlation = Vec::with_capacity(CANDIDATES);
    let mut against = 0i64;
    for slot in 0..CANDIDATES {
        let mut total = 0i64;
        for k in 0..SPAN {
            total = acc(total + mul(target[k], vectors[slot * SPAN + k]));
        }
        correlation.push(trunc32(total));
        against = against.max(acc(total).abs());
    }
    let up = if sat(against) == 0 {
        14
    } else {
        crate::fixed::exp(against).min(15) - 1
    };
    for c in correlation.iter_mut() {
        *c = trunc32(shift(*c, up));
    }

    let mut energy = Vec::with_capacity(CANDIDATES);
    let mut largest = 0i64;
    for slot in 0..CANDIDATES {
        let mut total = 0i64;
        for k in 0..SPAN {
            let v = vectors[slot * SPAN + k];
            total = acc(total + mul(v, v));
        }
        let scaled = shift(total, -5);
        energy.push(trunc32(scaled));
        largest = largest.max(scaled);
    }
    let over = if sat(largest) == 0 {
        13
    } else {
        crate::fixed::exp(largest).min(15) - 2
    };
    for e in energy.iter_mut() {
        *e = trunc32(shift(*e, over));
    }

    (correlation, up as i16, energy, over as i16)
}

/// Choose one shape from each codebook.
///
/// The two shapes are scored together: the numerator is the squared sum of
/// their correlations with the target and the denominator the energy of the
/// pair, which is the two energies plus twice their cross term.  As everywhere
/// else the comparison is by cross-multiplication rather than division.
///
/// `up` is the energy shift plus one, which is what lines the cross term up
/// with the two energies.  A pair whose correlations do not add to anything
/// positive is passed over; if none survives, neither codebook contributes.
pub fn best_pair(
    first: &[i16],
    second: &[i16],
    correlation: &[i64],
    energy: &[i64],
    up: i16,
    numbers: (&[i16], &[i16]),
) -> ((i16, i16), (i16, i16)) {
    let mut best = (32767i64, 1i64);
    let mut chosen: Option<(usize, usize)> = None;

    for i in 0..SHORTLIST {
        for j in 0..SHORTLIST {
            // Both shapes' correlations are known before the cross term, and a
            // pair that cannot win on them cannot win at all, so the 80-tap
            // correlation below is only worth taking once that is settled.
            let sum = acc(correlation[i] + correlation[SHORTLIST + j]);
            if acc(sum) <= 0 {
                continue;
            }
            let mut cross = 0i64;
            for k in 0..SPAN {
                cross = acc(cross + mul(first[i * SPAN + k], second[j * SPAN + k]));
            }
            let total =
                acc(shift(cross, up as i32) + shift(acc(energy[i] + energy[SHORTLIST + j]), 5));
            let pair = hi(shift(total, -5)) as i64;

            if takes_lead(&mut best, sum, pair) {
                chosen = Some((i, j));
            }
        }
    }

    let score = (low(best.0), low(best.1));
    match chosen {
        Some((i, j)) => ((numbers.0[i], numbers.1[j]), score),
        None => ((0, 0), score),
    }
}

/// The shape codebook's best pair for one subframe.
pub fn shape_pair(
    input: &Search,
    scaled: &[i16; SPAN],
    residual: &[i16; SPAN],
    contribution: &[i16; SPAN],
    inverse: i16,
    left: i16,
    numbers: &mut [i16; 2 * SHORTLIST],
) -> ((i16, i16), (i16, i16), i16, i16) {
    let mut response = response(input.impulse, input.lag, input.fraction, input.extension.1);
    // The shapes are scored against the target the pulses were scored against,
    // not against the response's own autocorrelation.
    let products = lag_products(residual, &response);
    let first_numbers = shapes(&FIRST_SHAPES, &products);
    let second_numbers = shapes(&SECOND_SHAPES, &products);
    quarter(&mut response);

    let mut first = shape_vectors(&response, &first_numbers, &FIRST_SHAPES);
    let mut second = shape_vectors(&response, &second_numbers, &SECOND_SHAPES);
    orthogonalise(&mut first, contribution, inverse, left);
    orthogonalise(&mut second, contribution, inverse, left);

    let mut both = Vec::with_capacity(first.len() + second.len());
    both.extend_from_slice(&first);
    both.extend_from_slice(&second);
    let (correlation, up, energy, over) = measure_shapes(scaled, &both);
    numbers[..SHORTLIST].copy_from_slice(&first_numbers);
    numbers[SHORTLIST..].copy_from_slice(&second_numbers);
    let (chosen, score) = best_pair(
        &first,
        &second,
        &correlation,
        &energy,
        over + 1,
        (&first_numbers, &second_numbers),
    );
    (chosen, score, up, over)
}

/// The two shape tables the pair is drawn from.
use crate::tables::{FCB_ALT_A as FIRST_SHAPES, FCB_ALT_B as SECOND_SHAPES};
