//! Voice-activity detection.
//!
//! The band analysis leaves behind two parallel descriptions of the frame: how
//! far the current spectrum sits above each of two long-term noise references,
//! summarised both per band and as a whole.  The detector turns each of those
//! into a small integer score by comparing three features against a set of
//! thresholds, then runs the score through a counter/hangover pair so that a
//! single quiet frame in the middle of speech does not switch the encoder over.
//!
//! Both references are scored by exactly the same code with exactly the same
//! feature layout; only the thresholds and the hangover lengths differ, and
//! those are chosen per frame from how loud the frame is.  The second score
//! additionally takes the open-loop voicing measure into account, and it is the
//! second one that the rest of the encoder consumes.

use crate::bands::BANDS;
use crate::fixed::{acc, hi, low, shift};

/// Bands taking part in the excess sums; band 0 is left out.
const SCORED_BANDS: usize = 18;

/// Threshold on the tracked summary level that selects the loud parameter set.
const LOUD_TRACKED: i64 = 15360 << 16;
/// Threshold on the current summary level, used both to pick a parameter set
/// and as the gate on the hangover counters.
const LOUD_SUMMARY: i64 = 7680 << 16;

/// Thresholds voicing is compared against when adjusting the second score.
const VOICING_HIGH: i16 = 26214;
const VOICING_MID: i16 = 19661;
const VOICING_LOW: i16 = 13107;

/// Range the second score is clamped to before the hangover logic sees it.
const SCORE_LIMIT: i16 = 4;

/// One scoring configuration: an upper and a lower threshold for each of the
/// three features, plus the hangover arming delay and length.
#[derive(Clone, Copy)]
struct Thresholds {
    /// Upper/lower pairs for the excess sum, the summary difference and the
    /// SNR variance, in that order.  The last two are compared against the
    /// high half of a 32-bit feature.
    limit: [i16; 6],
    /// Frames of continuous activity before the hangover is armed.
    delay: i16,
    /// Frames the hangover keeps the decision alive once armed.
    length: i16,
}

/// Pick the thresholds for one reference from how loud the frame is.
///
/// Quiet frames get a harder test on the excess sum and the summary difference
/// but a much longer hangover, which is what keeps low-level speech from being
/// chopped up.
fn thresholds(tracked: i64, summary: i64) -> Thresholds {
    if acc(tracked - LOUD_TRACKED) > 0 {
        Thresholds {
            limit: [50, 15, 768, 128, 40, 6],
            delay: 4,
            length: 4,
        }
    } else if acc(summary - LOUD_SUMMARY) > 0 {
        Thresholds {
            limit: [50, 15, 768, 128, 40, 6],
            delay: 6,
            length: 2,
        }
    } else {
        Thresholds {
            limit: [80, 15, 1280, 128, 40, 6],
            delay: 8,
            length: 1,
        }
    }
}

/// The three features of one reference, in the order the detector tests them.
#[derive(Clone, Copy, Default)]
pub struct Features {
    /// Summary level of the long-term reference.
    pub tracked: i64,
    /// Sum of the per-band excess over the reference, coarsely quantised.
    pub excess_sum: i16,
    /// Clipped difference between the current and the tracked summary level.
    pub difference: i64,
    /// Variance of the per-band signal-to-noise ratios.
    pub variance: i64,
}

/// Everything the detector looks at for one frame.
#[derive(Clone, Copy, Default)]
pub struct Frame {
    /// Summary level of the current spectrum.
    pub summary: i64,
    /// Features of the fast and the slow noise reference.
    pub reference: [Features; 2],
    /// Open-loop voicing measure, Q15.
    pub voicing: i16,
}

/// Score one reference: +1 for every feature above its upper threshold, -1 for
/// every feature below its lower one, and nothing in between.
fn score(t: &Thresholds, f: &Features) -> i16 {
    let mut total = 0i16;
    let mut vote = |value: i64, upper: i64, lower: i64| {
        if acc(value - upper) > 0 {
            total += 1;
        } else if acc(value - lower) < 0 {
            total -= 1;
        }
    };
    vote(f.excess_sum as i64, t.limit[0] as i64, t.limit[1] as i64);
    vote(
        f.difference,
        (t.limit[2] as i64) << 16,
        (t.limit[3] as i64) << 16,
    );
    vote(
        f.variance,
        (t.limit[4] as i64) << 16,
        (t.limit[5] as i64) << 16,
    );
    total
}

/// Nudge a score by how periodic the frame looks.
///
/// Strong periodicity is evidence of voicing even when the spectral features
/// are unconvincing; the absence of it counts double against.
fn voicing_bias(voicing: i16) -> i16 {
    if voicing > VOICING_HIGH {
        2
    } else if voicing > VOICING_MID {
        1
    } else if voicing > VOICING_LOW {
        0
    } else {
        -2
    }
}

/// Hangover state of one reference.
#[derive(Clone, Copy, Default)]
struct Hangover {
    /// Consecutive frames scored active.
    active: i16,
    /// Frames left to hold the decision after activity stops.
    remaining: i16,
}

/// Outcome of one hangover step.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Held {
    /// The frame scored active on its own.
    Active,
    /// The frame did not, but the hangover is still holding the decision up.
    Coasting,
    /// Neither.
    Idle,
}

impl Hangover {
    /// Run one frame through the counter/hangover pair and return the decision.
    fn step(&mut self, loud: bool, score: i16, delay: i16, length: i16) -> Held {
        if loud && score >= 0 {
            self.active = self.active.wrapping_add(1);
            if self.active > delay {
                self.remaining = length;
            }
            Held::Active
        } else {
            self.active = 0;
            self.remaining = self.remaining.wrapping_sub(1);
            if self.remaining > 0 {
                Held::Coasting
            } else {
                self.remaining = 0;
                Held::Idle
            }
        }
    }
}

/// Persistent state of the detector.
#[derive(Clone, Copy, Default)]
pub struct Vad {
    hangover: [Hangover; 2],
    /// Last frame's reported score, decayed while the hangover runs.
    previous: i16,
}

/// What the detector reports to the rest of the encoder.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Decision {
    /// Set while the frame is treated as speech.
    pub active: bool,
    /// The same decision taken against the fast reference.  Nothing downstream
    /// reads it except the fast noise estimate, which uses it to decide whether
    /// this frame may update it.
    pub fast: bool,
    /// Confidence score, clamped to `-4..=4`.
    pub score: i16,
}

impl Vad {
    /// Run the detector over one frame's features.
    pub fn run(&mut self, frame: &Frame) -> Decision {
        let loud = acc(frame.summary - LOUD_SUMMARY) > 0;

        let t0 = thresholds(frame.reference[0].tracked, frame.summary);
        let s0 = score(&t0, &frame.reference[0]);

        let t1 = thresholds(frame.reference[1].tracked, frame.summary);
        let mut s1 = score(&t1, &frame.reference[1]);
        s1 += voicing_bias(frame.voicing);
        s1 = s1.clamp(-SCORE_LIMIT, SCORE_LIMIT);

        let fast = self.hangover[0].step(loud, s0, t0.delay, t0.length) != Held::Idle;

        let held = self.hangover[1].step(loud, s1, t1.delay, t1.length);

        // Once the second reference is only coasting on its hangover the
        // reported score stops coming from this frame and decays instead.
        if held == Held::Coasting {
            s1 = if self.previous > 0 {
                self.previous - 1
            } else {
                0
            };
        }
        self.previous = s1;
        Decision {
            active: held != Held::Idle,
            fast,
            score: s1,
        }
    }
}

/// Coarse sum of the per-band excess levels.
///
/// Each band contributes the top half of its 32-bit excess shifted down by
/// eight, so the sum stays well inside a word.  Band 0 is left out.
pub fn excess_sum(excess: &[i64; BANDS]) -> i16 {
    let mut total = 0i64;
    for &e in excess[1..=SCORED_BANDS].iter() {
        total = acc(total + shift(hi(e) as i64, -8));
    }
    low(total)
}
