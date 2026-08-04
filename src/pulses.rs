//! The fixed-codebook search, from shortlisting to a chosen
//! index.
//!
//! Every mode is tried in turn: its positions are shortlisted, the pulse
//! responses and their cross energies built, every surviving combination
//! scored, and the winner weighed against the best mode so far.  [`crate::shape_pair`] offers
//! its twin five-pulse codebook as one more alternative, and whatever wins is
//! packed into the transmitted index.
//!
//! # Shortlisting
//!
//! The fixed codebook places one pulse on each of two or three interleaved
//! tracks.  Searching every combination would be expensive, so each track is
//! first reduced to a handful of positions: the ones where the target
//! correlation is largest.  The shortlist is built by insertion — seeded with
//! the first few positions kept in order, then every remaining position on the
//! track is offered to the list and takes a place if it beats an entry already
//! there, pushing the weakest one out.
//!
//! Alongside each position the list records where on the track it came from,
//! which is what the codebook index is eventually built out of.

use crate::shape_pair::{SHORTLIST, shape_pair};
use crate::tables::{
    FCB_CLASS_BASE, FCB_POSITIONS, FCB_PULSES, FCB_RADIX, FCB_SHAPE_SEL, FCB_SHAPES, MODE_WEIGHTS,
};

/// Room each track's list is given, and the span of one track's positions.
pub const TRACK_STRIDE: usize = 20;
const POSITION_STRIDE: usize = 40;
/// Largest number of tracks any mode uses.
pub const TRACKS: usize = 3;

/// Words of [`FCB_POSITIONS`] each mode covers.
const MODE_STRIDE: usize = 120;

/// Tracks one mode uses.
pub fn tracks(mode: usize) -> usize {
    FCB_PULSES[mode] as usize
}

/// Positions kept from one track.
pub fn kept(mode: usize, track: usize) -> usize {
    crate::tables::FCB_SHORTLIST[3 * mode + track] as usize
}

/// One mode's shortlists, laid out track by track.
///
/// Only the first [`kept`] entries of each track are written; the reference
/// leaves whatever the previous call put in the rest, and nothing downstream
/// looks at it.
pub struct Shortlist {
    /// The chosen positions, best first.
    pub position: [i16; TRACKS * TRACK_STRIDE],
    /// Where each of them sits in its track's position list.
    pub rank: [i16; TRACKS * TRACK_STRIDE],
}

/// Shortlist the best positions on every track of one mode.
///
/// `value` is the correlation of the target with a pulse at each position.
pub fn shortlist(mode: usize, value: &[i16]) -> Shortlist {
    let mut out = Shortlist {
        position: [0; TRACKS * TRACK_STRIDE],
        rank: [0; TRACKS * TRACK_STRIDE],
    };

    for track in 0..tracks(mode) {
        let keep = kept(mode, track);
        let total = FCB_RADIX[3 * mode + track] as usize;
        let table = MODE_STRIDE * mode + POSITION_STRIDE * track;
        let base = TRACK_STRIDE * track;
        let better = |a: i16, b: i16| value[a as usize] < value[b as usize];

        // Seed: the first `keep` positions, sorted as they arrive.
        out.position[base] = FCB_POSITIONS[table];
        for k in 1..keep {
            let candidate = FCB_POSITIONS[table + k];
            out.position[base + k] = candidate;
            out.rank[base + k] = k as i16;
            for j in 0..k {
                if better(out.position[base + j], candidate) {
                    out.position.copy_within(base + j..base + k, base + j + 1);
                    out.rank.copy_within(base + j..base + k, base + j + 1);
                    out.position[base + j] = candidate;
                    out.rank[base + j] = k as i16;
                    break;
                }
            }
        }

        // The rest of the track, offered one at a time to a list that is now
        // full: an entry has to be pushed out for a new one to get in.
        for k in keep..total {
            let candidate = FCB_POSITIONS[table + k];
            let mut placed = false;
            for j in 0..keep - 1 {
                if better(out.position[base + j], candidate) {
                    out.position
                        .copy_within(base + j..base + keep - 1, base + j + 1);
                    out.rank
                        .copy_within(base + j..base + keep - 1, base + j + 1);
                    out.position[base + j] = candidate;
                    out.rank[base + j] = k as i16;
                    placed = true;
                    break;
                }
            }
            // Falling off the end still leaves the weakest entry to beat.
            if !placed && better(out.position[base + keep - 1], candidate) {
                out.position[base + keep - 1] = candidate;
                out.rank[base + keep - 1] = k as i16;
            }
        }
    }
    out
}

use crate::fixed::{acc, hi, shift, trunc32};

/// Which position each track ended up with, and whether the pulse is used.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Chosen {
    pub index: [i16; TRACKS],
    pub used: [i16; TRACKS],
    /// How the winning combination scored, as the reference leaves it: the
    /// pair's energy first and the squared correlation second.  The comparison
    /// against the other modes is a cross-multiplication of these two.
    pub score: (i16, i16),
}

/// The cross-multiplied ratio test every search in this module ranks by.
///
/// A candidate is described by the energy it costs and the correlation `sum` it
/// buys; the better of two is the one with the larger `sum^2 / energy`, which
/// is compared by cross multiplication so that no division is needed.  Updates
/// the running best in place and reports whether the candidate took the lead.
pub(crate) fn takes_lead(best: &mut (i64, i64), sum: i64, energy: i64) -> bool {
    if acc(sum) <= 0 {
        return false;
    }
    let squared = hi(acc((hi(sum) as i64) * (hi(sum) as i64) * 2)) as i64;
    if acc(best.0 * squared * 2 - energy * best.1 * 2) <= 0 {
        return false;
    }
    *best = (energy, squared);
    true
}

/// Search every pair of shortlisted positions.
///
/// The two pulses are scored together: the numerator is the squared sum of the
/// two correlations and the denominator the energy of the pair, which is the
/// two self energies plus the cross term.  Rather than divide, each candidate is
/// compared against the best so far by cross-multiplying.
///
/// `correlation` and `energy` hold both tracks' shortlists end to end, `cross`
/// the pair energies in row-major order, and `sign` says whether a pulse at the
/// chosen position would help at all.
pub fn search_pairs(
    keep: [usize; 2],
    correlation: &[i16],
    energy: &[i16],
    cross: &[i16],
    sign: &[i16],
) -> Chosen {
    // A ratio so small that any real candidate beats it.
    let mut best = (32767i64, 1i64);
    let mut out = Chosen {
        index: [0; TRACKS],
        used: [1; TRACKS],
        score: (0, 0),
    };

    for i in 0..keep[0] {
        // The first track's energy is rounded down to a multiple of four on the
        // way through a sixteen-bit slot.
        let first = shift(acc((energy[i] as i64) << 16), -2);
        for j in 0..keep[1] {
            let total = acc(((energy[keep[0] + j] as i64) << 16)
                + shift(acc((hi(first) as i64) << 16), 2)
                + ((cross[i * keep[1] + j] as i64) << 16));
            let sum = acc(((correlation[keep[0] + j] as i64) + (correlation[i] as i64)) << 16);
            if takes_lead(&mut best, sum, hi(shift(total, -2)) as i64) {
                out.index[0] = i as i16;
                out.index[1] = j as i16;
            }
        }
    }

    out.score = (low(best.0), low(best.1));
    drop_unhelpful(&mut out, &keep, sign);
    out
}

/// A pulse whose position turned out not to help is dropped.
fn drop_unhelpful(out: &mut Chosen, keep: &[usize], sign: &[i16]) {
    let mut at = 0usize;
    for (track, &n) in keep.iter().enumerate() {
        if acc(sign[at + out.index[track] as usize] as i64) <= 0 {
            out.used[track] = 0;
        }
        at += n;
    }
}

/// Search every triple of shortlisted positions.
///
/// The same cross-multiplied comparison as [`search_pairs`], one nesting level
/// deeper.  Each of the three energies is passed through a sixteen-bit slot on
/// its way down, so the running total is truncated to a multiple of four twice
/// before the final comparison.
///
/// `cross` is the pair-energy table as [`energies`] leaves it: for each position
/// on the first track, its energies against the second track then against the
/// third, followed by the second-against-third block.  The reference shuffles
/// that into three contiguous blocks first; indexing it in place comes to the
/// same thing.
pub fn search_triples(
    keep: [usize; 3],
    correlation: &[i16],
    energy: &[i16],
    cross: &[i16],
    sign: &[i16],
) -> Chosen {
    let mut best = (32767i64, 1i64);
    let mut out = Chosen {
        index: [0; TRACKS],
        used: [1; TRACKS],
        score: (0, 0),
    };

    let row = keep[1] + keep[2];
    let last = keep[0] * row;
    let quarter = |v: i64| shift(acc(v), -2);
    let restore = |v: i16| shift(acc((v as i64) << 16), 2);

    for i in 0..keep[0] {
        let first = hi(quarter(acc((energy[i] as i64) << 16)));
        for j in 0..keep[1] {
            let pair = acc(((energy[keep[0] + j] as i64) << 16)
                + restore(first)
                + ((cross[i * row + j] as i64) << 16));
            let second = hi(quarter(pair));
            let two = hi(acc(((correlation[i] as i64)
                + (correlation[keep[0] + j] as i64))
                << 16));

            for k in 0..keep[2] {
                let total = acc(((energy[keep[0] + keep[1] + k] as i64) << 16)
                    + restore(second)
                    + ((cross[i * row + keep[1] + k] as i64) << 16)
                    + ((cross[last + j * keep[2] + k] as i64) << 16));
                let sum = acc(((two as i64) + (correlation[keep[0] + keep[1] + k] as i64)) << 16);
                if takes_lead(&mut best, sum, hi(quarter(total)) as i64) {
                    out.index = [i as i16, j as i16, k as i16];
                }
            }
        }
    }

    out.score = (low(best.0), low(best.1));
    drop_unhelpful(&mut out, &keep, sign);
    out
}

use crate::fixed::{low, mul, sat};

/// Samples one pulse's filtered response spans.
const SPAN: usize = 80;

/// Correlations between the filtered pulse responses.
///
/// Each track has its own copy of the weighted impulse response, filtered
/// through the synthesis filter; a pulse at position `p` on that track
/// contributes the response shifted right by `p`.  What the search needs is the
/// energy of each such contribution and the cross energy of every pair drawn
/// from two different tracks, which is what this builds.
///
/// The responses are first scaled up to use the full word, and the shift that
/// took is reported so that the corrections applied afterwards can be lined up
/// with it.  The pair energies come out in the order the search reads them:
/// for each position on each track, its energies against every later track.
pub fn energies(
    keep: &[usize],
    filtered: &mut [i16],
    position: &[i16],
) -> (i16, Vec<i16>, Vec<i16>) {
    let tracks = keep.len();
    let span = SPAN * tracks;

    // Bring the responses up to the top of the word.
    let mut peak = 0i64;
    for &v in filtered[..span].iter() {
        peak = peak.max(acc((v as i64) << 16).abs());
    }
    let up = if sat(peak) == 0 {
        15
    } else {
        crate::fixed::exp(peak)
    };
    for v in filtered[..span].iter_mut() {
        *v = hi(shift(acc((*v as i64) << 16), up));
    }

    let at = |track: usize, i: usize| position[TRACK_STRIDE * track + i] as usize;
    let correlate = |x: &[i16], y: &[i16], n: usize| {
        let mut total = 0i64;
        for k in 0..n {
            total = acc(total + mul(x[k], y[k]));
        }
        hi(shift(total, -6))
    };

    // Every pair of positions on two different tracks, upper triangle.
    let pair_count: usize = (0..tracks)
        .map(|first| keep[first] * keep[first + 1..].iter().sum::<usize>())
        .sum();
    let mut pairs = Vec::with_capacity(pair_count);
    for first in 0..tracks - 1 {
        for i in 0..keep[first] {
            let p = at(first, i);
            for second in first + 1..tracks {
                for j in 0..keep[second] {
                    let q = at(second, j);
                    let overlap = p.max(q);
                    pairs.push(correlate(
                        &filtered[SPAN * first + overlap - p..],
                        &filtered[SPAN * second + overlap - q..],
                        SPAN - overlap,
                    ));
                }
            }
        }
    }

    // Then each position's own energy.
    let mut alone = Vec::with_capacity(keep.iter().sum());
    for (track, &kept_here) in keep.iter().enumerate() {
        for i in 0..kept_here {
            let base = SPAN * track;
            alone.push(correlate(
                &filtered[base..],
                &filtered[base..],
                SPAN - at(track, i),
            ));
        }
    }

    (low(shift(acc(up as i64), 1)), pairs, alone)
}

/// Remove the adaptive codebook's share and fold in the pulse signs.
///
/// The energies built by [`energies`] describe the pulse responses on their own.
/// What the search actually needs is their energy against the part of the target
/// the adaptive codebook has not already accounted for, so the adaptive
/// contribution is subtracted from each: `adaptive[m]` is that codebook's
/// correlation with the response at position `m` and `gain` the gain it was
/// given.  A zero gain leaves the energies alone.
///
/// Afterwards every pair energy picks up the signs of the two pulses it belongs
/// to, and the self energies are halved.
#[allow(clippy::too_many_arguments)]
pub fn correct(
    keep: &[usize],
    pairs: &mut [i16],
    alone: &mut [i16],
    adaptive: &[i16],
    sign: &[i16],
    gain: i16,
    gain_shift: i16,
    headroom: i16,
) {
    let total: usize = keep.iter().sum();
    let line = |carried: i16, term: i64| {
        let lifted = shift(acc((carried as i64) << 16), 4);
        hi(shift(
            acc(lifted - shift(shift(term, gain_shift as i32), headroom as i32)),
            -4,
        ))
    };

    if acc(gain as i64) > 0 {
        debug_assert!(total <= TRACKS * TRACK_STRIDE);
        let mut share = [0i16; TRACKS * TRACK_STRIDE];
        for m in 0..total {
            share[m] = hi(acc(mul(gain, adaptive[m])));
            alone[m] = line(alone[m], acc(mul(adaptive[m], share[m])));
        }
        let mut at = 0usize;
        let mut base = 0usize;
        for &kept_here in &keep[..keep.len() - 1] {
            for i in 0..kept_here {
                for &later in &adaptive[base + kept_here..total] {
                    pairs[at] = line(pairs[at], acc(mul(later, share[base + i])));
                    at += 1;
                }
            }
            base += kept_here;
        }
    }

    // Each pair energy carries the signs of both its pulses.
    let mut at = 0usize;
    let mut base = 0usize;
    for &kept_here in &keep[..keep.len() - 1] {
        for i in 0..kept_here {
            let flip = acc(sign[base + i] as i64) <= 0;
            for &later in &sign[base + kept_here..total] {
                let product = acc(mul(pairs[at], later));
                pairs[at] = hi(shift(if flip { acc(-product) } else { product }, 1));
                at += 1;
            }
        }
        base += kept_here;
    }

    for e in alone[..total].iter_mut() {
        *e = hi(shift(acc((*e as i64) << 16), -1));
    }
}

/// Positions each two-track mode covers.
/// Which track each position in turn belongs to.
const THREE_TRACK_ORDER: [usize; 5] = [0, 1, 1, 2, 2];
const TWO_TRACK_ORDER: [usize; 2] = [0, 1];
/// Positions a three-track mode covers.
const THREE_TRACK_POSITIONS: usize = 40;
/// Positions in a subframe.
const POSITIONS: usize = 80;

/// How far apart one mode's candidate positions sit, and how many there are.
pub fn grid(mode: usize) -> (usize, usize) {
    if tracks(mode) == 3 {
        (mode + 1, THREE_TRACK_POSITIONS)
    } else {
        (1, crate::tables::FCB_POSITION_TOTAL[mode] as usize)
    }
}

/// What the target correlates with at every candidate pulse position.
///
/// A pulse at position `p` on track `t` contributes that track's filtered
/// response shifted right by `p`, so its correlation with the target is the
/// target from `p` onwards against the response from the start.  Positions are
/// dealt out to the tracks in turn, and how far apart they sit depends on the
/// mode.
///
/// What comes back is the magnitude of each correlation, normalised to use the
/// full word, the sign it had, the correlation of the same response with the
/// adaptive codebook's contribution, and the shift the normalisation took.
///
/// `adaptive` is written in place and only where `active` says the adaptive
/// codebook contributed anything, and then only on the mode's own grid — the
/// rest keeps what the previous call left there, which is what the reference
/// does and what the gather downstream goes on to read.
pub fn correlations(
    mode: usize,
    target: &[i16],
    contribution: &[i16],
    responses: &[i16],
    gain: i16,
    active: i16,
    adaptive: &mut [i16; POSITIONS],
) -> ([i16; POSITIONS], [i16; POSITIONS], i16) {
    let order: &[usize] = if tracks(mode) == 3 {
        &THREE_TRACK_ORDER
    } else {
        &TWO_TRACK_ORDER
    };
    let (stride, count) = grid(mode);

    let against = |source: &[i16], p: usize, track: usize, down: i32| {
        let mut total = 0i64;
        for n in 0..POSITIONS - p {
            total = acc(total + mul(source[p + n], responses[SPAN * track + n]));
        }
        shift(total, down)
    };

    let mut correlation = [0i64; POSITIONS];
    for k in 0..count {
        let (p, track) = (k * stride, order[k % order.len()]);
        correlation[p] = trunc32(against(target, p, track, -2));
    }

    // The adaptive codebook has already taken its share of the target.
    if acc(active as i64) > 0 {
        for k in 0..count {
            let (p, track) = (k * stride, order[k % order.len()]);
            let share = sat(against(contribution, p, track, -1));
            adaptive[p] = hi(share);
            correlation[p] = trunc32(acc(correlation[p] - shift(acc(mul(gain, hi(share))), 1)));
        }
    }

    // Signs off, then everything scaled up together.
    let mut sign = [0i16; POSITIONS];
    let mut peak = 0i64;
    for p in 0..POSITIONS {
        sign[p] = if correlation[p] < 0 { -16384 } else { 16384 };
        correlation[p] = correlation[p].abs();
        peak = peak.max(correlation[p]);
    }
    let up = if sat(peak) == 0 {
        13
    } else {
        crate::fixed::exp(peak).min(15) - 2
    };

    let mut value = [0i16; POSITIONS];
    for p in 0..POSITIONS {
        value[p] = hi(shift(correlation[p], up));
    }
    (value, sign, up as i16)
}

/// Pack the chosen positions and signs into a codebook index.
///
/// The tracks have different numbers of positions, so the indices are packed
/// mixed-radix rather than by bit fields, with each track's position count as
/// the next radix.  The sign bits go in underneath, one per track, and the
/// mode's own base is added on top so that every mode occupies its own stretch
/// of the codebook.
pub fn codebook_index(mode: usize, index: &[i16], used: &[i16]) -> i16 {
    let n = tracks(mode);

    let mut signs = 0i64;
    for t in (0..n).rev() {
        signs = acc(shift(signs, 1) + (used[t] as i64));
    }

    // The multiply runs with the fractional bit off, so it is a plain product.
    let mut packed = index[0] as i64;
    for (t, &position) in index.iter().enumerate().take(n).skip(1) {
        let radix = FCB_RADIX[3 * mode + t] as i64;
        packed = acc(shift(acc(packed * radix), 16) + ((position as i64) << 16));
        packed = shift(packed, -16);
    }

    low(acc(shift(packed, n as i32)
        + signs
        + (FCB_CLASS_BASE[mode] as i64)))
}

/// Scale the target up to use the full word.
///
/// Three bits of headroom are left so that the correlations built on top of it
/// cannot overflow; a silent subframe gets a fixed shift instead.
pub fn normalise_target(target: &[i16]) -> [i16; SPAN] {
    let mut peak = 0i64;
    for &v in target.iter().take(SPAN) {
        peak = peak.max(acc((v as i64) << 16).abs());
    }
    let peak = sat(peak);
    let up = if peak == 0 {
        12
    } else {
        crate::fixed::exp(peak) - 3
    };

    let mut out = [0i16; SPAN];
    for (o, &v) in out.iter_mut().zip(target.iter()) {
        *o = hi(shift(acc((v as i64) << 16), up));
    }
    out
}

/// What the adaptive codebook already accounts for, taken out of the target.
///
/// The contribution is scaled up by half its own reciprocal's exponent, the
/// gain that best fits it to the target is worked out, and that multiple of it
/// is subtracted, leaving the residual the fixed codebook has to cover.
///
/// Returns the reciprocal's mantissa, the shift left over, and the gain — the
/// three words the rest of the search reads.  A silent contribution leaves the
/// residual untouched and reports zeros.
pub fn remove_adaptive(
    target: &[i16],
    contribution: &mut [i16; SPAN],
    residual: &mut [i16; SPAN],
) -> (i16, i16, i16) {
    let mut energy = 0i64;
    for &v in contribution.iter() {
        energy = acc(energy + mul(v, v));
    }
    if acc(energy) <= 0 {
        return (0, 0, 0);
    }

    let (quotient, exponent) = crate::fixed::normalised_reciprocal(energy);
    let inverse = hi(quotient);
    let up = shift(acc(exponent as i64), -1) - 1;
    for v in contribution.iter_mut() {
        *v = hi(shift(acc((*v as i64) << 16), up as i32));
    }
    let left = low(acc((exponent as i64) - shift(up, 1)));

    let mut cross = 0i64;
    for k in (0..SPAN).rev() {
        cross = acc(cross + mul(target[k], contribution[k]));
    }
    let scaled = crate::fixed::mul32x16(trunc32(sat(cross)), inverse);
    let gain = hi(shift(shift(scaled, left as i32), -2));

    for k in 0..SPAN {
        let term = shift(acc(mul(contribution[k], gain)), 2);
        residual[k] = hi(acc(((target[k] as i64) << 16) - term));
    }
    (inverse, left, gain)
}

/// Which pulse shape each track of each mode uses.
/// The pulse shapes themselves, twenty taps each.
const SHAPE_TAPS: usize = 20;

/// Filter one shape per track through the weighted synthesis filter.
///
/// Each track has its own pulse shape; filtering it once gives the response a
/// pulse anywhere on that track will produce, which is what the whole search is
/// built on.  The result is extended periodically for the same reason the
/// adaptive codebook vector is.
pub fn track_responses(
    mode: usize,
    impulse: &[i16],
    lag: i16,
    fraction: i16,
    gain: i16,
) -> Vec<i16> {
    let n = tracks(mode);
    let mut out = vec![0i16; n * SPAN];
    for track in 0..n {
        let shape = FCB_SHAPE_SEL[3 * mode + track] as usize;
        let taps = &FCB_SHAPES[SHAPE_TAPS * shape..][..SHAPE_TAPS];
        let mut filtered = crate::convolve::filter(taps, impulse);
        crate::pitch::extend(&mut filtered, lag, fraction, gain);
        out[track * SPAN..(track + 1) * SPAN].copy_from_slice(&filtered);
    }
    out
}

/// Gather what the search needs about each shortlisted position.
///
/// The shortlist holds positions; everything downstream wants the three
/// quantities worked out for them, laid out consecutively track
/// after track rather than indexed by position.
pub fn gather(
    mode: usize,
    positions: &[i16],
    value: &[i16],
    adaptive: &[i16],
    sign: &[i16],
) -> (Vec<i16>, Vec<i16>, Vec<i16>) {
    let total: usize = (0..tracks(mode)).map(|t| kept(mode, t)).sum();
    let (mut values, mut shares, mut signs) = (
        Vec::with_capacity(total),
        Vec::with_capacity(total),
        Vec::with_capacity(total),
    );
    for track in 0..tracks(mode) {
        for i in 0..kept(mode, track) {
            let at = positions[TRACK_STRIDE * track + i] as usize;
            values.push(value[at]);
            shares.push(adaptive[at]);
            signs.push(sign[at]);
        }
    }
    (values, shares, signs)
}

/// Turn the chosen shortlist slots into positions on their tracks.
///
/// The search picks entries of the shortlist; what the index packing wants is
/// where those entries sit on their own track, which the shortlist recorded
/// alongside them.
pub fn ranks(mode: usize, chosen: &mut [i16], list: &[i16]) {
    for track in 0..tracks(mode) {
        chosen[track] = list[TRACK_STRIDE * track + chosen[track] as usize];
    }
}

/// Keep the better of this mode and the best so far.
///
/// Each mode leaves a pair standing for how well it did; comparing two of them
/// is again a cross-multiplication.  The running best starts at a pair nothing
/// can lose to.  When a mode wins, its pulse positions and sign flags are
/// copied aside, because the next mode is about to overwrite the working ones.
#[allow(clippy::too_many_arguments)]
pub fn keep_best(
    mode: usize,
    pair: (i16, i16),
    best: &mut (i16, i16, i16),
    chosen: &[i16],
    used: &[i16],
    positions: &mut [i16],
    flags: &mut [i16],
) -> bool {
    let test = acc(mul(pair.0, best.0) - mul(pair.1, best.1));
    if test >= 0 {
        return false;
    }
    // The pair is stored the other way round from how it is compared.
    *best = (pair.1, pair.0, mode as i16);
    let n = tracks(mode);
    positions[..n].copy_from_slice(&chosen[..n]);
    flags[..n].copy_from_slice(&used[..n]);
    true
}

/// Weigh one mode's result and line the pair up.
///
/// Each mode is worth a different amount — a mode that spends more bits has to
/// do better to be chosen — so its score is multiplied by a weight read from a
/// table that advances one entry per mode.  The pair is then brought onto one
/// exponent so that the comparison against the running best is a plain
/// cross-multiplication; whichever of the two needs shifting down gets it, and
/// a gap of more than sixteen bits drops the second word altogether.
pub fn weigh_mode(pair: &mut (i16, i16), weight: i16, correlation: i16, energy: i16) {
    let weighted = acc(mul(pair.0, weight));
    let e = crate::fixed::exp(weighted);
    pair.0 = hi(shift(weighted, e));
    align_pair(pair, e, correlation, energy);
}

/// Bring the two words of a scored pair onto one exponent.
fn align_pair(pair: &mut (i16, i16), e: i32, correlation: i16, energy: i16) {
    let gap = acc((e as i64) - shift(correlation as i64, 1) + (energy as i64));
    if gap > 0 {
        pair.0 = hi(shift(acc((pair.0 as i64) << 16), -gap as i32));
    } else if acc(gap + 16) < 0 {
        pair.1 = 0;
    } else {
        pair.1 = hi(shift(acc((pair.1 as i64) << 16), gap as i32));
    }
}

/// Weight of the shape-pair codebook, and the factor its first mode's weight is
/// kept at for the comparison that follows.
const PAIR_KEPT: i16 = 27306;

/// Weigh the shape pair's result the same way.
///
/// The weight table has one more entry than there are modes; that last one
/// belongs to the shape pair.  On the way past, the table's first entry is
/// scaled and left in the slot the pointer occupied, ready for the comparison
/// mode 2 makes.
pub fn weigh_pair(
    pair: &mut (i16, i16),
    weight: i16,
    first: i16,
    correlation: i16,
    energy: i16,
) -> i16 {
    weigh_mode(pair, weight, correlation, energy);
    hi(acc(mul(first, PAIR_KEPT)))
}

/// How far the innovation's centre may drift from the target's before the mode
/// is penalised, and how much of its score it keeps when it does.
const DRIFT: i64 = 8;
const KEPT_SCORE: i16 = 22938;

/// Penalise a mode whose pulses land in the wrong place.
///
/// A mode can score well by putting its energy in the wrong part of the
/// subframe, which sounds worse than the numbers suggest.  Comparing where the
/// innovation's rectified area balances against where the target's does catches
/// that; drifting more than eight samples costs the mode three tenths of its
/// score.
pub fn centroid_penalty(score: &mut i16, headroom: &mut i16, innovation: i64, target: i16) {
    let drift = acc(innovation - (target as i64)).abs();
    if acc(drift - DRIFT) <= 0 {
        return;
    }
    *score = hi(acc(mul(*score, KEPT_SCORE)));
    *headroom = low(acc((*headroom as i64) - 1));
}

/// Bits the first shape number is shifted up by.
const PAIR_SHIFT: i32 = 7;

/// Choose between the pulse codebook and the shape-pair codebook.
///
/// The two are scored the same way and compared by cross-multiplication one
/// last time.  If the shape pair wins, the index is rebuilt out of the two
/// shape numbers and the base its part of the codebook starts at.
pub fn choose_innovation(
    pulses: (i16, i16),
    pair: (i16, i16),
    shapes: (i16, i16),
    index: &mut i16,
) -> bool {
    let test = acc(mul(pulses.0, pair.0) - mul(pulses.1, pair.1));
    if test >= 0 {
        return false;
    }
    *index = low(acc(shift(shapes.0 as i64, PAIR_SHIFT)
        + (shapes.1 as i64)
        + (crate::tables::FCB_CLASS_BASE[crate::tables::FCB_PAIR_CLASS]
            as i64)));
    true
}

/// Try the residual on its own against the shape pair.
///
/// Only for subframes the mode decision called generic.  Coding nothing but the
/// residual's own energy is the fallback the shape pair has to beat; if it does
/// not, the pair is replaced by it.
pub fn plain_alternative(pair: &mut (i16, i16), kept: i16, residual: &[i16]) -> bool {
    let mut energy = 0i64;
    for &v in residual.iter().take(SPAN) {
        energy = acc(energy + mul(v, v));
    }
    let e = crate::fixed::exp(energy);
    let mut normalised = shift(energy, e);

    let gap = acc((e as i64) - 3);
    let mut kept = kept;
    if gap > 0 {
        normalised = shift(normalised, -gap as i32);
    } else if acc(gap + 16) < 0 {
        kept = 0;
    } else {
        kept = hi(shift(acc((kept as i64) << 16), gap as i32));
    }

    if acc(mul(pair.0, hi(normalised)) - mul(kept, pair.1)) <= 0 {
        return false;
    }
    *pair = (kept, hi(normalised));
    true
}

/// Modes the fixed codebook search tries.
pub const MODES: usize = 5;

/// The gain each mode's periodic extension uses.
///
/// Two gains are kept between subframes; which of them a mode gets depends on
/// how that mode places its pulses.  A lag as long as the subframe leaves
/// nothing to extend, so every mode gets nothing.
pub fn mode_gains(lag: i16, first: i16, second: i16) -> [i16; MODES] {
    if acc((lag as i64) - SPAN as i64) >= 0 {
        [0; MODES]
    } else {
        [first, second, second, first, first]
    }
}

/// Lay the chosen pulses into an innovation vector.
///
/// Each track contributes its own response, shifted to the chosen position and
/// added or subtracted according to the pulse's sign.  The decoder reaches the
/// same routine with the codebook's shape tables instead of filtered responses.
pub fn place_pulses(mode: usize, chosen: &[i16], used: &[i16], responses: &[i16]) -> [i16; SPAN] {
    let mut out = [0i16; SPAN];
    for track in 0..tracks(mode) {
        let table = MODE_STRIDE * mode + POSITION_STRIDE * track;
        let at = FCB_POSITIONS[table + chosen[track] as usize] as usize;
        let base = SPAN * track;
        for k in 0..SPAN - at {
            let carried = (out[at + k] as i64) << 16;
            let term = (responses[base + k] as i64) << 16;
            out[at + k] = hi(acc(if used[track] != 0 {
                carried + term
            } else {
                carried - term
            }));
        }
    }
    out
}

/// Weights per quantiser mode in [`MODE_WEIGHTS`]: one for each of the five
/// pulse modes and one for the shape pair.
const WEIGHTS_PER_ROW: usize = 6;

/// What one subframe's fixed codebook search settles on.
///
/// Only the index comes out: the vector itself is built from it by the same
/// routine the decoder uses, so that the two cannot drift apart.
pub struct Innovation {
    /// The codebook index that goes into the bitstream.
    pub index: i16,
    /// How each mode scored as it went into the comparison, which is what the
    /// choice between them is made on.
    pub scores: [(i16, i16); MODES],
    /// How the shape pair scored against them.
    pub pair: (i16, i16),
    /// The two shape codebooks' shortlists.
    pub numbers: [i16; 2 * SHORTLIST],
}

/// Everything the search needs from the rest of the subframe.
pub struct Search<'a> {
    /// The weighted target, and the adaptive codebook's filtered contribution.
    pub target: &'a [i16],
    pub contribution: &'a [i16],
    /// The weighted synthesis filter's impulse response.
    pub impulse: &'a [i16],
    /// The lag the adaptive codebook settled on, and its fraction.
    pub lag: i16,
    pub fraction: i16,
    /// The two extension gains carried between subframes.
    pub extension: (i16, i16),
    /// How the subframe is being coded, and which row of the weight table the
    /// line spectrum's own mode selects.
    pub coding: i16,
    pub weights: usize,
}

/// Search the fixed codebook.
///
/// Five pulse modes are tried in turn, each laying one pulse on each of its two
/// or three tracks, and the best of them is compared against a pair drawn from
/// the shape codebook.  Every comparison is a cross-multiplication of a score
/// and an energy rather than a division, so the two are carried side by side
/// the whole way down.
///
/// `adaptive` carries the per-position correlations against the adaptive
/// codebook between calls: each mode only refreshes the positions on its own
/// grid, and the rest keeps what the last mode left there.
pub fn search(input: &Search, adaptive: &mut [i16; POSITIONS]) -> Innovation {
    let centre = centroid(input.target);
    let gains = mode_gains(input.lag, input.extension.0, input.extension.1);

    let scaled = normalise_target(input.target);
    let mut contribution = [0i16; SPAN];
    contribution.copy_from_slice(&input.contribution[..SPAN]);
    let mut residual = scaled;
    let (inverse, left, gain) = remove_adaptive(&scaled, &mut contribution, &mut residual);

    // A pair nothing can lose to, and where the winner's pulses are kept while
    // the next mode overwrites the working ones.
    let mut best = (-32768i16, 0i16, -1i16);
    let (mut positions, mut flags) = ([0i16; TRACKS], [0i16; TRACKS]);
    let mut scores = [(0i16, 0i16); MODES];

    for mode in 0..MODES {
        let mut filtered =
            track_responses(mode, input.impulse, input.lag, input.fraction, gains[mode]);
        let (value, sign, up) = correlations(
            mode,
            &scaled,
            &contribution,
            &filtered,
            gain,
            inverse,
            adaptive,
        );
        let list = shortlist(mode, &value);
        let (values, shares, signs) = gather(mode, &list.position, &value, adaptive, &sign);

        let mut keep = [0usize; TRACKS];
        for (t, k) in keep.iter_mut().enumerate().take(tracks(mode)) {
            *k = kept(mode, t);
        }
        let keep = &keep[..tracks(mode)];
        let (mut headroom, mut pairs, mut alone) = energies(keep, &mut filtered, &list.position);
        correct(
            keep, &mut pairs, &mut alone, &shares, &signs, inverse, left, headroom,
        );

        let mut chosen = if tracks(mode) == 3 {
            search_triples([keep[0], keep[1], keep[2]], &values, &alone, &pairs, &signs)
        } else {
            search_pairs([keep[0], keep[1]], &values, &alone, &pairs, &signs)
        };
        ranks(mode, &mut chosen.index, &list.rank);

        // A mode that puts its energy in the wrong part of the subframe is
        // worth less than its numbers say.
        let mut pair = chosen.score;
        if acc(input.coding as i64) > 0 {
            // The responses come out of the energy pass scaled up to the top of
            // the word; five of them added together would overflow, so they are
            // quartered before the innovation is laid out.
            for v in filtered.iter_mut() {
                *v = hi(shift(acc((*v as i64) << 16), -2));
            }
            let placed = place_pulses(mode, &chosen.index, &chosen.used, &filtered);
            let drift = centroid(&placed);
            centroid_penalty(&mut pair.0, &mut headroom, drift as i64, centre);
        }

        let weight = crate::tables::MODE_WEIGHTS[WEIGHTS_PER_ROW * input.weights + mode];
        weigh_mode(&mut pair, weight, up, headroom);
        scores[mode] = pair;
        keep_best(
            mode,
            pair,
            &mut best,
            &chosen.index,
            &chosen.used,
            &mut positions,
            &mut flags,
        );
    }

    let mut index = codebook_index(best.2 as usize, &positions, &flags);

    // The shape codebook's own best pair, which the pulses have to beat.
    let mut numbers = [0i16; 2 * SHORTLIST];
    let (shapes, mut pair, up, over) = shape_pair(
        input,
        &scaled,
        &residual,
        &contribution,
        inverse,
        left,
        &mut numbers,
    );
    let row = WEIGHTS_PER_ROW * input.weights;
    let first = MODE_WEIGHTS[row];
    let kept = weigh_pair(&mut pair, MODE_WEIGHTS[row + MODES], first, up, over);
    if acc(input.coding as i64) == 2 {
        plain_alternative(&mut pair, kept, &residual);
    }
    choose_innovation((best.0, best.1), pair, shapes, &mut index);
    Innovation {
        index,
        scores,
        pair,
        numbers,
    }
}

/// Where the first half of the subframe's absolute sum has gone by.
///
/// This is the point that splits the subframe into two halves of equal
/// rectified area — a cheap stand-in for where the energy sits, used to place
/// the codebook search's window.  A silent subframe reports the middle.
pub fn centroid(x: &[i16]) -> i16 {
    let magnitude = |v: i16| acc((v as i64) << 16).abs();
    let mut total = 0i64;
    for &v in x.iter().take(SPAN) {
        total = acc(total + magnitude(v));
    }
    let mut half = shift(total, -1);
    if acc(half) == 0 {
        return (SPAN / 2) as i16;
    }
    for (i, &v) in x.iter().enumerate().take(SPAN) {
        half = acc(half - magnitude(v));
        if half <= 0 {
            return i as i16;
        }
    }
    (SPAN - 1) as i16
}
