//! Encoder input conditioning.
//!
//! The encoder takes 320 mu-law samples, expands them to linear and runs them
//! through a second-order high-pass that removes DC and the lowest part of the
//! band before any analysis happens.

use crate::fixed::{acc, hi, sat, shift, trunc32};

/// Samples per frame.
pub use crate::FRAME;
/// Samples per half-frame, which is the unit the analysis works on.
pub use crate::HALF;

/// Expand one mu-law byte to linear.
///
/// A plain table lookup, against the same table the reference uses.
pub fn ulaw_to_linear(sample: u8) -> i16 {
    crate::tables::ULAW_TO_LINEAR[sample as usize]
}

/// State of the input high-pass filter.
#[derive(Clone, Copy, Default)]
pub struct Highpass {
    /// The two previous outputs, kept at full 32-bit precision.
    output: [i64; 2],
    /// The two previous inputs.
    input: [i16; 2],
}

impl Highpass {
    /// Filter one half-frame in place.
    ///
    /// `y[n] = 8 * (a1 y[n-1] + a2 y[n-2] + b0 x[n] + b1 x[n-1] + b2 x[n-2])`,
    /// where the recursive part is evaluated against the undamaged 32-bit
    /// history rather than the rounded output.
    pub fn run(&mut self, block: &mut [i16]) {
        let c = &crate::tables::INPUT_HIGHPASS[..5];
        for sample in block.iter_mut() {
            let x = *sample;
            let mut a = acc(recursive(self.output[0], c[0]) + recursive(self.output[1], c[1]));
            a = acc(a + (x as i64) * (c[2] as i64) * 2);
            a = acc(a + (self.input[0] as i64) * (c[3] as i64) * 2);
            a = acc(a + (self.input[1] as i64) * (c[4] as i64) * 2);
            a = sat(shift(a, 3));

            self.input[1] = self.input[0];
            self.input[0] = x;
            self.output[1] = self.output[0];
            self.output[0] = trunc32(a);

            // The rounding add saturates, so a sample that overshoots full scale
            // clips there instead of wrapping round to the opposite rail.
            *sample = hi(sat(acc(a + 32768)));
        }
    }
}

/// One recursive tap: a 32x16 product built from the halves of the state, the
/// low half being pre-shifted so the two partial products line up.
fn recursive(state: i64, coef: i16) -> i64 {
    let low = ((state as i16 as u16) >> 1) as i64;
    let partial = shift(acc(low * (coef as i64) * 2), -16);
    acc(shift(partial, 1) + (hi(state) as i64) * (coef as i64) * 2)
}
