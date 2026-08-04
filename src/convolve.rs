//! Filtering the search's shapes through the impulse response.
//!
//! Both encoder searches work against a signal that has already been through
//! the weighted synthesis filter, so rather than filtering every candidate they
//! convolve each shape with the impulse response once and combine the results.
//! This is that convolution; the searches themselves are in
//! [`crate::pitch_search`] and [`crate::pulses`].

use crate::fixed::{acc, hi, sat, shift};

/// Samples in a subframe, and taps in one pulse shape.
pub use crate::SUBFRAME;
pub const SHAPE: usize = 20;

/// Convolve one shape with the impulse response.
///
/// The shape's length is whatever the caller passes: the fixed codebook's
/// twenty-tap pulses, or the eighty samples of excitation history the closed
/// loop searches against.
pub fn filter(shape: &[i16], response: &[i16]) -> [i16; SUBFRAME] {
    let mut out = [0i16; SUBFRAME];
    for k in 0..SUBFRAME {
        let mut total = 0i64;
        for i in 0..=k.min(shape.len() - 1) {
            total = acc(total + (shape[i] as i64) * (response[k - i] as i64) * 2);
        }
        out[k] = hi(sat(shift(total, 3)));
    }
    out
}
