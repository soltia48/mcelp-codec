//! Deciding how a subframe is to be coded.
//!
//! Before the fixed codebook is searched the encoder settles on one of three
//! modes by asking two questions.  The first is whether the adaptive codebook
//! alone already accounts for most of the target: if the predicted energy
//! covers roughly two thirds of it, the subframe is periodic enough that mode 0
//! is used.  Otherwise the second question is whether the target has just grown
//! sharply against the previous subframe, which is what an onset looks like;
//! onsets get mode 1, and so do the three subframes after one, so that the
//! attack is not cut short.  Everything else is mode 2.

use crate::fixed::{acc, mul, mul32x16 as mul32, sat, shift, trunc32};

/// Fraction of the target energy the prediction has to reach for mode 0.
const COVERAGE: i16 = 20972;
/// Halved because the comparison shifts back up by two afterwards: an onset is
/// a subframe carrying twice the energy of the one before it.
const ONSET_RATIO: i16 = 16384;
/// Subframes the onset decision is held for.
const HOLD: i16 = 3;

/// Headroom taken out of both energies before they are compared.
const HEADROOM: i32 = -7;

/// What the rest of the search does with the subframe.
pub const PERIODIC: i16 = 0;
pub const ONSET: i16 = 1;
pub const GENERIC: i16 = 2;

/// State the decision carries between subframes.
#[derive(Clone, Copy, Default)]
pub struct Mode {
    /// Energy of the previous subframe's target, at the same headroom.
    previous: i64,
    /// Subframes still to be held at [`ONSET`].
    hold: i16,
}

impl Mode {
    /// Choose the mode for one subframe.
    ///
    /// `target` is the weighted speech and `filtered` the adaptive codebook's
    /// contribution to it; `gain` is the gain chosen for that contribution.
    pub fn choose(&mut self, target: &[i16], filtered: &[i16], gain: i16) -> i16 {
        let predicted = trunc32(shift(energy(filtered), HEADROOM));
        // The gain is applied twice because it is the energy that is scaled.
        let predicted = trunc32(mul32(predicted, gain));
        let predicted = mul32(predicted, gain);
        let predicted = shift(sat(shift(predicted, 2)), 2);

        let current = trunc32(shift(energy(target), HEADROOM));
        let wanted = mul32(current, COVERAGE);

        if acc(predicted - wanted) >= 0 {
            self.previous = current;
            self.decay();
            return PERIODIC;
        }

        let threshold = shift(mul32(self.previous, ONSET_RATIO), 2);
        self.previous = current;
        if acc(current - threshold) > 0 {
            self.hold = HOLD;
            ONSET
        } else if self.hold > 0 {
            self.hold -= 1;
            ONSET
        } else {
            GENERIC
        }
    }

    /// Count one subframe off the onset hold.
    fn decay(&mut self) {
        if self.hold > 0 {
            self.hold -= 1;
        }
    }
}

/// Energy of one subframe.
fn energy(block: &[i16]) -> i64 {
    let mut total = 0i64;
    for &v in block.iter() {
        total = acc(total + mul(v, v));
    }
    total
}

/// Accessors used by the trace replay to line the state up with the reference.
impl Mode {
    pub fn restore(&mut self, previous: i64, hold: i16) {
        self.previous = previous;
        self.hold = hold;
    }
    pub fn energy(&self) -> i64 {
        self.previous
    }
    pub fn hold(&self) -> i16 {
        self.hold
    }
}
