//! Fixed (innovation) codebook.
//!
//! The 16-bit codebook index is first classified: six index ranges select six
//! different codebook structures.  Classes 0..=4 are multi-pulse — the index
//! within the class is a mixed-radix number giving two or three pulse positions
//! plus their signs, and each pulse contributes a 20-sample waveform taken from
//! a shape table.  Class 5 is a twin five-pulse codebook where each half of the
//! index selects five signed unit pulses.

use crate::fixed::{acc, hi, shift};
use crate::pitch::SUBFRAME;

/// Number of index classes.
const CLASSES: usize = 6;
/// The class that uses the alternative twin-pulse codebook.
const ALT_CLASS: usize = 5;
/// Pulses per half of the alternative codebook.
const ALT_PULSES: usize = 5;
/// Amplitude of an alternative-codebook pulse, Q15.
const ALT_AMPLITUDE: i16 = 8192;
/// Length of one pulse shape.
const SHAPE_LEN: usize = 20;
/// Candidate positions per pulse slot.
const POSITIONS: usize = 40;
/// Maximum pulses in a multi-pulse class.
const MAX_PULSES: usize = 3;

/// A decoded innovation vector plus the class it came from.
pub struct Innovation {
    /// The innovation itself.
    pub code: [i16; SUBFRAME],
    /// Index class, which also selects which pitch-sharpening gain applies.
    pub class: usize,
}

/// Pick the structure a received index selects: the largest class whose lower
/// bound the index clears.
fn classify(index: u16) -> (usize, u16) {
    let index = index.min(crate::tables::FCB_INDEX_LIMIT[0] - 1);
    let class = (0..CLASSES)
        .rev()
        .find(|&k| index >= crate::tables::FCB_CLASS_BASE[k])
        .unwrap_or(0);
    (class, index)
}

/// Split a class index into pulse positions and signs.
///
/// The low bits carry one sign per pulse; what remains is a mixed-radix number
/// whose digits are the pulse position indices, most significant digit first.
fn decompose(class: usize, index: u16, pulses: usize) -> ([i16; MAX_PULSES], [i16; MAX_PULSES]) {
    let mut rest = (index as i64) - (crate::tables::FCB_CLASS_BASE[class] as i64);

    let mut signs = [0i16; MAX_PULSES];
    for s in signs.iter_mut().take(pulses) {
        *s = (rest & 1) as i16;
        rest = shift(rest, -1);
    }

    let mut positions = [0i16; MAX_PULSES];
    let mut value = rest << 16;
    for k in (0..pulses).rev() {
        let recip = crate::tables::FCB_RADIX_INV[MAX_PULSES * class + k];
        let radix = crate::tables::FCB_RADIX[MAX_PULSES * class + k];
        let quotient = hi(shift(acc((hi(value) as i64) * (recip as i64) * 2), -2));
        value = acc(value - (((quotient as i64) * (radix as i64) * 2) << 15));
        positions[k] = hi(value);
        value = (quotient as i64) << 16;
    }
    (positions, signs)
}

/// Lay the pulse shapes into the innovation vector.
///
/// A zero sign bit subtracts the shape, a set one adds it.
fn place_shapes(class: usize, positions: &[i16], signs: &[i16], pulses: usize) -> [i16; SUBFRAME] {
    let first_shape = crate::tables::FCB_SHAPE_SEL[MAX_PULSES * class] as usize;
    let mut code = [0i16; SUBFRAME];
    for k in 0..pulses {
        let slot = positions[k] as usize;
        let pos =
            crate::tables::FCB_POSITIONS[POSITIONS * (MAX_PULSES * class + k) + slot] as usize;
        let shape = &crate::tables::FCB_SHAPES
            [SHAPE_LEN * (first_shape + k)..SHAPE_LEN * (first_shape + k) + SHAPE_LEN];
        for n in pos..SUBFRAME.min(pos + SHAPE_LEN) {
            let x = (code[n] as i64) << 16;
            let y = (shape[n - pos] as i64) << 16;
            code[n] = hi(if signs[k] != 0 {
                acc(x + y)
            } else {
                acc(x - y)
            });
        }
    }
    code
}

/// The twin five-pulse codebook used by the top index class.
///
/// Each table entry is a signed 1-based position: the magnitude places the
/// pulse, the sign chooses its polarity.
fn alternative(index: u16) -> [i16; SUBFRAME] {
    let rest = (index as i64) - (crate::tables::FCB_CLASS_BASE[ALT_CLASS] as i64);
    let upper = shift(rest, -7);
    let lower = rest - shift(upper, 7);

    let mut code = [0i16; SUBFRAME];
    for (base, offset) in [
        (&crate::tables::FCB_ALT_A, upper as usize),
        (&crate::tables::FCB_ALT_B, lower as usize),
    ] {
        for k in 0..ALT_PULSES {
            let entry = base[ALT_PULSES * offset + k];
            // Positions are 1-based; the sign of the entry is the pulse's sign.
            let pos = entry.unsigned_abs() as usize - 1;
            let x = (code[pos] as i64) << 16;
            let amp = (ALT_AMPLITUDE as i64) << 16;
            code[pos] = hi(if entry <= 0 {
                acc(x - amp)
            } else {
                acc(x + amp)
            });
        }
    }
    code
}

/// Decode one subframe's innovation vector.
pub fn decode(index: u16) -> Innovation {
    let (class, index) = classify(index);
    if class == ALT_CLASS {
        return Innovation {
            code: alternative(index),
            class,
        };
    }
    let pulses = crate::tables::FCB_PULSES[class] as usize;
    let (positions, signs) = decompose(class, index, pulses);
    Innovation {
        code: place_shapes(class, &positions, &signs, pulses),
        class,
    }
}
