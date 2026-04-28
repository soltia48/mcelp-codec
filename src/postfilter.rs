//! Postfilter (`func_c86e`).
//!
//! Operates on half-frame units (160 samples). Two IIR stages — forward and
//! reverse — perform formant emphasis + tilt + output-gain correction.

use crate::arith::{acc32_16, mpya_from_acc32_16, sat32};
use crate::tables::postfilter::{FEEDBACK_COEFF, INPUT_COEFF, OUTPUT_GAIN};

/// `round_acc(b)` — `(b<<2) sat32 + 0x8000 sat32`.
fn round_acc(acc: i64) -> i64 {
    let scaled = sat32(acc << 2);
    sat32(scaled + 0x8000)
}

/// 4-tap MAC + 3-cell history shift register.
fn input_step(input_acc: i64, history: &mut [i16; 3], coeffs: &[i16; 4]) -> i64 {
    let mut acc = mpya_from_acc32_16(coeffs[0], input_acc, 2);
    acc = sat32(acc + i64::from(history[0]) * i64::from(coeffs[1]) * 2);
    acc = sat32(acc + i64::from(history[1]) * i64::from(coeffs[2]) * 2);
    acc = sat32(acc + i64::from(history[2]) * i64::from(coeffs[3]) * 2);

    history[2] = history[1];
    history[1] = history[0];
    history[0] = ((input_acc >> 16) & 0xffff) as i16;
    let _ = acc32_16; // helper used inside mpya_from_acc32_16
    acc
}

/// 3 IIR feedback sections.
fn feedback_correction(mut acc: i64, delay: &[i16; 6], coeffs: &[i16; 3]) -> i64 {
    for s in 0..3 {
        let hi = i64::from(delay[2 * s]);
        let lo_u = delay[2 * s + 1] as u16 as i64;
        let dword = hi * 65536 + lo_u;
        acc -= ((dword * i64::from(coeffs[s]) * 2) >> 16) << 1;
    }
    acc
}

/// 3-element shift register on dword pairs.
fn shift_delay_line(delay: &mut [i16; 6], acc: i64) {
    let filtered = sat32(acc << 1);
    let f_hi = ((filtered >> 16) & 0xffff) as i16;
    let f_lo = (filtered & 0xffff) as i16;
    delay[4] = delay[2];
    delay[5] = delay[3];
    delay[2] = delay[0];
    delay[3] = delay[1];
    delay[0] = f_hi;
    delay[1] = f_lo;
}

fn output_sample(acc: i64, output_gain: i16) -> i64 {
    let rounded = round_acc(acc);
    let filtered = mpya_from_acc32_16(output_gain, rounded, 2);
    let scaled = sat32(filtered << 1);
    sat32(scaled + 0x8000)
}

/// 160-sample half-frame postfilter.
pub fn postfilter_apply(
    pre_samples: &[i16; 160],
    history_6: &mut [i16; 6],
    delay_12: &mut [i16; 12],
) -> [i16; 160] {
    let mut fwd_hist: [i16; 3] = [history_6[0], history_6[1], history_6[2]];
    let mut rev_hist: [i16; 3] = [history_6[3], history_6[4], history_6[5]];
    let mut delay = *delay_12;

    let mut output: [i16; 160] = *pre_samples;
    for sample in output.iter_mut() {
        let input_acc: i64 = (*sample as i64) << 16;

        // Forward stage
        let b1 = input_step(input_acc, &mut fwd_hist, &INPUT_COEFF);
        let fwd_fb: [i16; 6] = std::array::from_fn(|k| delay[k]);
        let b2 = feedback_correction(b1, &fwd_fb, &FEEDBACK_COEFF);
        let mut fwd_shift: [i16; 6] = std::array::from_fn(|k| delay[k]);
        shift_delay_line(&mut fwd_shift, b2);
        delay[..6].copy_from_slice(&fwd_shift);

        let rounded = round_acc(b2);

        // Reverse stage
        let b3 = input_step(rounded, &mut rev_hist, &INPUT_COEFF);
        let rev_fb: [i16; 6] = std::array::from_fn(|k| delay[6 + k]);
        let b4 = feedback_correction(b3, &rev_fb, &FEEDBACK_COEFF);
        let mut rev_shift: [i16; 6] = std::array::from_fn(|k| delay[6 + k]);
        shift_delay_line(&mut rev_shift, b4);
        delay[6..].copy_from_slice(&rev_shift);

        // Output gain
        let pf = output_sample(b4, OUTPUT_GAIN);
        *sample = ((pf >> 16) & 0xffff) as i16;
    }

    history_6[0..3].copy_from_slice(&fwd_hist);
    history_6[3..6].copy_from_slice(&rev_hist);
    *delay_12 = delay;
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_input_produces_zero_output_with_zero_state() {
        let zero = [0i16; 160];
        let mut hist = [0i16; 6];
        let mut delay = [0i16; 12];
        let out = postfilter_apply(&zero, &mut hist, &mut delay);
        assert!(out.iter().all(|&v| v == 0));
        assert!(hist.iter().all(|&v| v == 0));
        assert!(delay.iter().all(|&v| v == 0));
    }
}
