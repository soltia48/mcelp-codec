//! Short-term synthesis filter.

use crate::fixed::{acc, hi, mul, sat, shift};
use crate::lsp::ORDER;

/// Longest block any caller filters in one go, which is one subframe.
const BLOCK: usize = 80;

/// Run `1/A(z)` over `input`, writing `output` and carrying `mem` forward.
///
/// `a` is the direct-form predictor in Q12 (`a[0]` is unused here); the signal
/// is Q15, so the accumulated prediction is shifted left by three to line the
/// two formats up.  The reference has two entry points into it — the same loop, called
/// with different lengths.
pub fn synthesis(a: &[i16; ORDER + 1], input: &[i16], output: &mut [i16], mem: &mut [i16; ORDER]) {
    // A sliding window holding the ten previous outputs followed by the new
    // ones, which is how the filter memory is addressed.  No caller runs
    // the filter over more than one subframe at a time.
    debug_assert!(input.len() <= BLOCK);
    let mut storage = [0i16; ORDER + BLOCK];
    let window = &mut storage[..ORDER + input.len()];
    window[..ORDER].copy_from_slice(mem);

    for n in 0..input.len() {
        let mut sum = 0i64;
        for i in 1..=ORDER {
            sum = acc(sum - mul(a[i], window[ORDER + n - i]));
        }
        let y = hi(sat(acc(shift(sum, 3)
            + ((input[n] as i64) << 16)
            + (1 << 15))));
        window[ORDER + n] = y;
        output[n] = y;
    }

    mem.copy_from_slice(&window[window.len() - ORDER..]);
}
