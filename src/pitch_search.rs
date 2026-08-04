//! The encoder's pitch search, open loop and closed.
//!
//! The open loop runs first, over the whole weighted signal: three passes at
//! decreasing lag ranges, each keeping the lag whose correlation is largest
//! once normalised by its own energy.  What it
//! settles on brackets the closed loop, which is not run over the whole range
//! every time: the first subframe of a pair looks in a window
//! around the open-loop estimate and the second in a window around the first
//! one's choice, which is what makes the second lag fit in five bits.
//!
//! The closed loop's own arithmetic sits at the end of the module: the
//! correlations it scores lags on and
//! the adaptive gain that falls out of the winner.  Everything that happens
//! *after* a lag has been chosen is in [`crate::excitation`].
use crate::SUBFRAME;
use crate::fixed::{acc, divide, exp, hi, inverse_sqrt, low, mul, mul32, shift};
use crate::tables::{
    LAG_INTERP, LAG_INTERP_BACKWARD as INTERP_BACKWARD, LAG_INTERP_FORWARD as INTERP_FORWARD,
    LAG_SLOTS, LAG_SLOTS_HIGH, LAG_SLOTS_LOW,
};

/// Narrowest and widest lag the search will consider.
const SHORTEST: i16 = 20;
const LONGEST: i16 = 143;
/// Lag above which the window stops being stretched.
const SPLIT: i16 = 85;
/// How far below the previous lag the window starts.
const BACK: i16 = 5;
/// Width of the window, in the two regimes.
const WIDE: i16 = 117;
const NARROW: i16 = 43;
/// Half-width of the second subframe's window.
const RELATIVE: i16 = 2;

/// The window one subframe searches in.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Window {
    /// Inclusive bounds, as the reference stores them.
    pub low: i16,
    pub high: i16,
    /// Where the coded value is measured from.
    pub origin: i16,
}

/// Window for the first subframe of a pair, from the previous frame's lag.
pub fn first(previous: i16, base: i16) -> Window {
    let back = acc(previous as i64 - BACK as i64);
    let origin = if acc(previous as i64 - SPLIT as i64) > 0 {
        low(acc(back + WIDE as i64))
    } else {
        low(acc(acc(back + (back << 1)) - NARROW as i64 + base as i64))
    };
    let high = low(acc(back.max(SHORTEST as i64) + 9).min(LONGEST as i64));
    Window {
        low: low(acc(high as i64 - 9)),
        high,
        origin,
    }
}

/// Window for the second subframe, from the first one's.
pub fn second(low_bound: i16, previous: i16, base: i16) -> i16 {
    let d = acc(previous as i64 - low_bound as i64);
    low(acc(acc(d + (d << 1)) + RELATIVE as i64 + base as i64))
}

/// Level a peak has to clear.  The correlations themselves are passed in.
const THRESHOLD: i64 = 30000 << 16;
/// How far below the lag the window starts.
const LOOKBACK: i16 = 90;

/// Decide whether the band around the current lag is strongly periodic.
///
/// The lag picks out a run of open-loop correlations; if twice the largest of
/// them clears the threshold the caller treats the subframe as voiced.
pub fn periodic(lag: i16, fraction: i16, correlations: &[i64]) -> bool {
    // The `add #1` sits in the branch's two delay-slot words, so it happens
    // either way; only the `sub #1` that takes it back is skipped.
    let adjusted = if fraction > 0 {
        acc(lag as i64 + 1)
    } else {
        lag as i64
    };
    let high = LAG_SLOTS[LAG_SLOTS_HIGH + low(acc(adjusted)) as u16 as usize];
    let start = low(acc(adjusted - LOOKBACK as i64).max(0));
    let low_slot = LAG_SLOTS[LAG_SLOTS_LOW + start as u16 as usize];

    let mut peak = 0i64;
    for k in 0..=(high - low_slot) {
        peak = peak.max(correlations[(low_slot + k) as usize]);
    }
    acc(shift(peak, 1) - THRESHOLD) > 0
}

/// Samples the open-loop search correlates over.
pub const SEARCH_SPAN: usize = 343;

/// Scale the signal up before the open-loop correlation.
///
/// The correlation squares its input, so the signal is only lifted by half the
/// headroom it has - one bit less, to leave room for the accumulation.
pub fn normalise(signal: &[i16; SEARCH_SPAN]) -> [i16; SEARCH_SPAN] {
    let mut energy = 0i64;
    for &x in signal.iter() {
        energy = acc(energy + (x as i64) * (x as i64) * 2);
    }
    let up = shift(acc(crate::fixed::exp(energy) as i64 - 1), -1) as i32;
    let mut out = [0i16; SEARCH_SPAN];
    for (o, &x) in out.iter_mut().zip(signal.iter()) {
        *o = low(shift(acc((x as i64) << 16), up - 16));
    }
    out
}

/// Where the correlation window starts inside the normalised block, and how
/// long it is.
const WINDOW_START: usize = 143;
const WINDOW: usize = 200;

/// Find the lag with the largest correlation over a range.
///
/// `longest` is the first lag tried and `count` how many are tried below it.
/// The scan runs downwards and the winner is recorded whenever the correlation
/// merely equals the running maximum, so ties go to the shortest lag.
pub fn best_lag(signal: &[i16], longest: i16, count: i16) -> i16 {
    let mut best = acc((-32768i64) << 16);
    let mut winner = longest;
    for step in 0..=count {
        let lag = (longest - step) as usize;
        let mut total = 0i64;
        for i in 0..WINDOW {
            let a = signal[WINDOW_START + i] as i64;
            let b = signal[WINDOW_START + i - lag] as i64;
            total = acc(total + a * b * 2);
        }
        best = best.max(total);
        if acc(total - best) >= 0 {
            winner = longest - step;
        }
    }
    winner
}

/// Correlation of the best lag in a range, normalised by its own energy.
///
/// Returns the lag and the normalised peak, which is what the three passes are
/// compared on.
pub fn normalised_peak(signal: &[i16], longest: i16, count: i16) -> (i16, i64) {
    let lag = best_lag(signal, longest, count);

    let mut correlation = 0i64;
    let mut energy = 0i64;
    for i in 0..WINDOW {
        let a = signal[WINDOW_START + i] as i64;
        let b = signal[WINDOW_START + i - lag as usize] as i64;
        correlation = acc(correlation + a * b * 2);
        energy = acc(energy + b * b * 2);
    }
    (lag, shift(mul32(correlation, inverse_sqrt(energy)), 15))
}

/// The three ranges the open-loop search covers, longest first.
const PASSES: [(i16, i16); 3] = [(143, 63), (79, 39), (39, 19)];
/// How much better a shorter lag has to be to displace a longer one.
const BIAS: i16 = 27853;

/// Pick the open-loop pitch lag.
///
/// Three ranges are searched separately and the shorter one only takes over if
/// its normalised peak still reaches most of the longer one's — which is what
/// stops the search settling on a multiple of the true period.
pub fn open_loop_lag(signal: &[i16; SEARCH_SPAN]) -> i16 {
    let scaled = normalise(signal);
    let mut best = 0i16;
    let mut winner = 0i16;
    for (n, &(longest, count)) in PASSES.iter().enumerate() {
        let (lag, value) = normalised_peak(&scaled, longest, count);
        let candidate = crate::fixed::hi(value);
        if n == 0 || acc(((candidate as i64) << 16) - (BIAS as i64) * (best as i64) * 2) > 0 {
            best = candidate;
            winner = lag;
        }
    }
    winner
}

/// Bias added to every scaled peak, and how many entries the history shifts.
const PEAK_BIAS: i16 = 16384;
const HISTORY: usize = 3;

/// Track the largest recent correlation peak.
///
/// Short lags look only at the newest entry and its own rescaling; longer ones
/// sweep the slot range the lag maps onto, the same range `periodic` uses.  The
/// winner is pushed onto the front of the history.
pub fn track_peak(peaks: &mut [i64], lag: i16, gain: i16) -> i64 {
    let scale = |v: i64| -> i64 {
        crate::fixed::sat(acc(
            shift(crate::fixed::mul32(v, (gain as i64) << 16), 1) + PEAK_BIAS as i64
        ))
    };
    let mut best = -1i64;
    if acc(lag as i64 - 80) < 0 {
        let first = scale(peaks[0]);
        best = best.max(first);
        best = best.max(scale(first));
    } else {
        let low_slot = LAG_SLOTS[LAG_SLOTS_LOW + low(acc(lag as i64 - 80)) as u16 as usize];
        let high = LAG_SLOTS[LAG_SLOTS_LOW + low(acc(lag as i64 - 1)) as u16 as usize];
        for slot in low_slot..=high {
            best = best.max(scale(peaks[slot as usize]));
        }
    }
    for i in (0..HISTORY).rev() {
        peaks[i + 1] = peaks[i];
    }
    peaks[0] = best;
    best
}

/// Pick the best integer lag from the closed-loop correlations.
///
/// The running maximum is compared with `>=`, and the position is only recorded
/// when the new value wins outright — so among equals the *last* one is kept.
pub fn best_closed_lag(correlations: &[i16], low_bound: i16, high_bound: i16) -> i16 {
    let mut best = 0usize;
    let count = acc(high_bound as i64 - low_bound as i64 - 1);
    let mut running = 0i64;
    for i in 0..=count as usize {
        let here = (correlations[i] as i64) << 16;
        let next = (correlations[i + 1] as i64) << 16;
        if i == 0 {
            running = here;
        }
        let top = running.max(next);
        if acc(next - top) == 0 {
            best = i + 1;
        }
        running = top;
    }
    low_bound + best as i16
}

/// Interpolate the correlation between two integer lags.
///
/// `fraction` selects which of the sub-sample positions to evaluate; a negative
/// one is folded onto the next lag down, which is why the caller adjusts the
/// integer lag afterwards.  The caller tries all five of -2..2.
pub fn interpolate(correlations: &[i16], centre: usize, fraction: i16) -> i64 {
    let (mut at, mut phase) = (centre, fraction);
    if acc(fraction as i64) < 0 {
        phase += 3;
        at -= 1;
    }
    let mut total = 0i64;
    for tap in 0..3 {
        let back = LAG_INTERP[INTERP_FORWARD + phase as usize + 3 * tap];
        let forward = LAG_INTERP[INTERP_BACKWARD - phase as usize + 3 * tap];
        total = acc(total + (correlations[at - tap] as i64) * (back as i64) * 2);
        total = acc(total + (correlations[at + 1 + tap] as i64) * (forward as i64) * 2);
    }
    let last_back = LAG_INTERP[INTERP_FORWARD + phase as usize + 9];
    total = acc(total + (correlations[at - 3] as i64) * (last_back as i64) * 2);
    let last_forward = LAG_INTERP[INTERP_BACKWARD - phase as usize + 9];
    acc(acc(total + (correlations[at + 4] as i64) * (last_forward as i64) * 2) + 0x8000) & !0xffff
}

/// Shortest and longest lag the closed-loop search will consider.
const CLOSED_SHORTEST: i64 = 20;
const CLOSED_LONGEST: i64 = 143;
/// Lags either side of the open-loop estimate.
const SPREAD: i64 = 6;

/// Bracket the closed-loop search around the open-loop lag.
///
/// The open-loop estimate is only approximate, so the closed-loop search looks
/// at a short stretch of lags around it — three below and three above, before
/// the ends of the range are clamped.  Returns the first and last lag.
pub fn search_range(open_loop: i16) -> (i16, i16) {
    let start = acc((open_loop as i64) - 3).max(CLOSED_SHORTEST);
    let last = acc(start + SPREAD).min(CLOSED_LONGEST);
    (
        crate::fixed::low(acc(last - SPREAD)),
        crate::fixed::low(last),
    )
}

/// A value kept as a mantissa and the shift that was applied to it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Scaled {
    pub mantissa: i16,
    pub exponent: i16,
}

/// Bias added to each exponent; the first is doubled because it describes a
/// sum of squares.
const ENERGY_BIAS: i16 = 26;
const CORRELATION_BIAS: i16 = 13;

/// Energy and two correlations, normalised.
///
/// `before` and `after` are the two neighbouring stretches of the reference
/// signal; a zero result is replaced by the smallest non-zero value so that the
/// divisions downstream never see one.
pub fn measures(x: &[i16], before: &[i16], after: &[i16]) -> [Scaled; 3] {
    let mut energy = 0i64;
    let mut back = 0i64;
    let mut forward = 0i64;
    for i in 0..SUBFRAME {
        energy = acc(energy + (x[i] as i64) * (x[i] as i64) * 2);
        back = acc(back + (before[i] as i64) * (x[i] as i64) * 2);
        forward = acc(forward + (after[i] as i64) * (x[i] as i64) * 2);
    }
    [
        normalise_pair(energy, ENERGY_BIAS),
        normalise_pair(acc(-shift(back, 1)), CORRELATION_BIAS),
        normalise_pair(shift(forward, 1), CORRELATION_BIAS),
    ]
}

/// Split a value into mantissa and exponent, or the placeholder if it is zero.
fn normalise_pair(value: i64, bias: i16) -> Scaled {
    if acc(value) == 0 {
        return Scaled {
            mantissa: 1,
            exponent: 15,
        };
    }
    let e = exp(value);
    Scaled {
        mantissa: hi(shift(value, e)),
        exponent: bias + e as i16,
    }
}

/// Headroom above which the target is left alone, and the two shifts applied
/// below it.
const KEEP_ABOVE: i32 = 5;
const DEEP_BELOW: i32 = -4;

/// Scale the target down before the closed-loop correlations.
///
/// A loud subframe would overflow the correlation accumulator, so it is shifted
/// down by two bits, or four when it is louder still.  A quiet one is left as
/// it is.
pub fn prescale(target: &mut [i16; SUBFRAME]) -> i32 {
    let mut energy = 0i64;
    for &x in target.iter() {
        energy = acc(energy + (x as i64) * (x as i64) * 2);
    }
    let headroom = exp(energy);
    if headroom > KEEP_ABOVE {
        return 0;
    }
    let down = if headroom > DEEP_BELOW { 2 } else { 4 };
    for x in target.iter_mut() {
        *x = low(shift(acc((*x as i64) << 16), -16 - down));
    }
    -down
}

/// Cross-correlation of the two reference stretches, and the negated double of
/// it as a mantissa/exponent pair.
///
/// Returns the rounded correlation, the pair, and the headroom the second
/// block lines its own energy up against - which is the scaling the
/// correlation was taken at, not the pair's exponent.
pub fn cross(x: &[i16], y: &[i16]) -> (i16, Scaled, i16) {
    let mut total = 0i64;
    for i in 0..SUBFRAME {
        total = acc(total + (x[i] as i64) * (y[i] as i64) * 2);
    }
    let headroom = exp(total);
    let scaled = shift(shift(total, headroom), -2);
    let rounded = hi(acc(scaled + (1 << 15)));

    let against = (headroom - 2) as i16;
    let doubled = acc(-shift(scaled, 1));
    (rounded, normalise_pair(doubled, against), against)
}

/// Largest gain the adaptive codebook is allowed.
const GAIN_LIMIT: i64 = 19661 << 17;

/// Gain of the adaptive codebook and the energy it was divided by.
///
/// `correlation` and `correlation_exponent` come from the first block; the
/// quotient is lined up with them and clipped to just above unity, then halved
/// on the way out.
pub fn adaptive_gain(y: &[i16], correlation: i16, correlation_exponent: i16) -> (i16, Scaled) {
    let mut energy = 0i64;
    for &v in y.iter().take(SUBFRAME) {
        energy = acc(energy + (v as i64) * (v as i64) * 2);
    }
    if acc(energy) == 0 {
        return (
            0,
            Scaled {
                mantissa: 1,
                exponent: 15,
            },
        );
    }
    let e = exp(energy) as i16;
    let normalised = shift(energy, e as i32);
    let pair = Scaled {
        mantissa: hi(normalised),
        exponent: e,
    };

    let quotient = divide(normalised, correlation);
    if acc(quotient) <= 0 {
        return (0, pair);
    }

    // Line the quotient up with the correlation's own exponent, dropping a byte
    // or a word first when the gap is too wide for the shift field.
    let gap = acc(e as i64 - correlation_exponent as i64);
    let (shifted, amount) = if acc(gap - 8) > 0 {
        (shift(quotient, -8), acc(gap - 8))
    } else if acc(gap + 16) >= 0 {
        (quotient, gap)
    } else {
        (shift(quotient, -16), acc(gap + 16))
    };
    let limited = shift(acc(shifted), amount as i32).min(GAIN_LIMIT);
    (hi(shift(limited, -1)), pair)
}

/// Taps of the impulse response the recursive update convolves with.
const RESPONSE: usize = 40;

/// Normalised correlation of the target against each candidate lag.
///
/// The candidates are consecutive lags, so rather than filtering each one from
/// scratch the routine slides the filtered vector along by a sample and folds
/// in the one excitation sample that has just come into range:
///
/// ```text
/// y[n] = y[n-1] + x[-lag] h[n]
/// ```
///
/// Only the first forty terms need the correction; beyond that the response has
/// run out and the slide is all there is to it.  `past[k]` is the excitation
/// sample the k-th slide brings in, and `down` is the shift [`prescale`] chose,
/// applied here too so that the correlations stay in the target's units.
pub fn lag_correlations(
    target: &[i16],
    filtered: &mut [i16; SUBFRAME],
    past: &[i16],
    response: &[i16],
    down: i32,
    lags: usize,
) -> Vec<i16> {
    let mut out = Vec::with_capacity(lags);
    let mut history = past.iter();
    for lag in 0..lags {
        let mut energy = 0i64;
        for &v in filtered.iter() {
            energy = acc(energy + mul(v, v));
        }
        let inverse = inverse_sqrt(energy);

        let mut correlation = 0i64;
        for (&x, &y) in target.iter().zip(filtered.iter()) {
            correlation = acc(correlation + mul(x, y));
        }
        correlation = shift(correlation, down);
        out.push(hi(shift(mul32(inverse, correlation), 15)));

        if lag + 1 == lags {
            break;
        }
        let excitation = *history
            .next()
            .expect("past must cover every lag but the last");
        for n in (RESPONSE..SUBFRAME).rev() {
            filtered[n] = filtered[n - 1];
        }
        for n in (1..RESPONSE).rev() {
            let tap = shift(acc(mul(excitation, response[n])), down + 3);
            filtered[n] = hi(acc(((filtered[n - 1] as i64) << 16) + tap));
        }
        filtered[0] = hi(shift(acc((excitation as i64) << 16), down));
    }
    out
}
