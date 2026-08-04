//! Filtering of the shaping signal.
//!
//! What comes out of the noise suppression is still a raw spectral estimate;
//! before the search can use it, it goes through three filters in a row.  The
//! first of them is the same second-order high-pass the encoder puts its input
//! through, so the shaping signal is band-limited the same way the speech is.

use crate::fixed::{acc, hi, mul32x16, round, sat, shift};

/// State of the shaping signal's high-pass.
///
/// The two previous outputs are kept at 32-bit precision, as is the case
/// throughout: rounding them to the emitted word would let the recursion drift.
#[derive(Clone, Copy, Default)]
pub struct Highpass {
    output: [i64; 2],
    input: [i16; 2],
}

impl Highpass {
    /// Filter one half-frame in place.
    pub fn run(&mut self, block: &mut [i16]) {
        let ff = &crate::tables::HIGHPASS_FF[..3];
        let fb = &crate::tables::HIGHPASS_FB[..2];
        for sample in block.iter_mut() {
            let x = *sample;
            let mut a = mul32x16(self.output[0], fb[0]);
            a = acc(a + mul32x16(self.output[1], fb[1]));
            a = shift(a, 1);
            a = acc(a + (x as i64) * (ff[0] as i64) * 2);
            let taps = acc((self.input[1] as i64) * (ff[2] as i64) * 2
                + (self.input[0] as i64) * (ff[1] as i64) * 2);
            a = acc(a + taps);

            self.input[1] = self.input[0];
            self.input[0] = x;
            self.output[1] = self.output[0];
            self.output[0] = sat(shift(a, 2));

            *sample = hi(sat(acc(shift(a, 3) + 32768)));
        }
    }
}

/// Symmetric six-tap smoothing filter.
/// Feedback coefficients of the fifth-order recursion that follows it.
/// Taps of each.
const FIR_TAPS: usize = 6;
const IIR_TAPS: usize = 5;
/// Coefficient of the final differencing pass.
const TILT: i16 = 9830;

/// State of the shaping signal's smoothing filter.
#[derive(Clone, Copy, Default)]
pub struct Smooth {
    /// The five previous inputs of the symmetric part.
    input: [i16; FIR_TAPS - 1],
    /// The five previous outputs, at 32-bit precision.
    output: [i64; IIR_TAPS],
}

impl Smooth {
    /// Filter one half-frame in place.
    pub fn run(&mut self, block: &mut [i16]) {
        let ff = &crate::tables::SMOOTH_FIR[..FIR_TAPS];
        let fb = &crate::tables::SMOOTH_IIR[..IIR_TAPS];
        for sample in block.iter_mut() {
            let x = *sample;
            // Overflow mode is off through the whole shaping path, so only the
            // explicit saturations below clip; the accumulations run in the
            // full forty bits.
            let mut b = acc((x as i64) * (ff[0] as i64) * 2);
            for (tap, &past) in self.input.iter().enumerate() {
                b = acc(b + (past as i64) * (ff[tap + 1] as i64) * 2);
            }
            self.input.rotate_right(1);
            self.input[0] = x;

            for (tap, &past) in self.output.iter().enumerate() {
                b = acc(b - shift(mul32x16(past, fb[tap]), 1));
            }
            self.output.rotate_right(1);
            self.output[0] = sat(shift(b, 2));

            *sample = hi(sat(acc(shift(b, 3) + 32768)));
        }
    }
}

/// State of the final differencing pass.
#[derive(Clone, Copy, Default)]
pub struct Tilt {
    /// The last sample of the previous half-frame.
    previous: i16,
}

impl Tilt {
    /// Apply `y[n] = x[n] - 0.3 x[n-1]` in place.
    ///
    /// The reference runs it backwards so that it can work in place without a
    /// copy; the result is the same as a forward pass over the untouched input,
    /// which is how it is written here.
    pub fn run(&mut self, block: &mut [i16]) {
        let mut previous = self.previous;
        self.previous = block[block.len() - 1];
        for sample in block.iter_mut() {
            let x = *sample;
            *sample = hi(sat(round(acc(
                ((x as i64) << 16) - (TILT as i64) * (previous as i64) * 2
            ))));
            previous = x;
        }
    }
}
