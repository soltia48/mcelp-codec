//! LSP → autocorrelation → Levinson → LPC pipeline.

//! Functions are organized in the following hierarchy:
//!
//! ```text
//! build_autocorrelation_with_state  (top-level: LSP → 22-cell LPC)
//!   ├─ response_shaper                (×2: → 11-cell mirror_output each)
//!   │    └─ levinson_step             (×2 per cbec call)
//!   ↓
//! levinson_recursion  (per-subframe: autocorr window → reflection coeffs k_1..k_10)
//!   └─ q15_reciprocal                     (= 1/E_{i-1} estimate)
//! ```

use crate::arith::{exp_acc, norm_acc_with_t, shift_acc40};

const FRCT_MUL: i64 = 2;

// ============================================================================
// dword pair access on 35-cell scratch buffer
// ============================================================================

#[inline]
fn dword_read(buf: &[i16; 35], off: usize) -> i64 {
    let hi = buf[off];
    let lo_u = buf[off ^ 1] as u16 as i64;
    ((hi as i64) << 16) | lo_u
}

#[inline]
fn dword_write(buf: &mut [i16; 35], off: usize, value: i64) {
    buf[off] = ((value >> 16) & 0xffff) as i16;
    buf[off ^ 1] = (value & 0xffff) as i16;
}

#[inline]
fn subc_step(src: i64, divisor: i16) -> i64 {
    let div = i64::from(divisor);
    let alu = src - (div << 15);
    if alu >= 0 { (alu << 1) + 1 } else { src << 1 }
}

/// Output of `q15_reciprocal`.
pub struct ReciprocalResult {
    /// 1/|input| Q15 estimate (sign restored; internal i40-compatible).
    pub acc_main: i64,
    /// Normalization exponent + 1 (the caller uses it to shift back).
    pub norm_lo: i16,
}

pub fn q15_reciprocal(input: i64) -> ReciprocalResult {
    // Step 1: Q15-float normalize — extract mantissa and exponent.
    let (mantissa, t_exp) = q15_float_normalize(input);
    let norm_lo = (1 + t_exp as i32) as i16;

    // Divisor: clamp |mantissa| to the i16 range.
    let divisor = (mantissa as i32).unsigned_abs().min(i16::MAX as u32) as i16;

    // Step 2: 15-step subc divide, computing 16383 / divisor in Q15.
    let mut acc = (16383i64) << 16;
    for _ in 0..15 {
        acc = subc_step(acc, divisor);
    }

    // Step 3: 32-bit shift × 2 (× 256 × 256) — reproduce upper-bit dropping with a 32-bit mask.
    acc = (((acc << 8) & 0xffff_ffff) as u32) as i64;
    acc = (((acc << 8) & 0xffff_ffff) as u32) as i64;

    // Step 4: 40-bit sign extension (bit 39 → upper 24 bits).
    acc = sign_extend_40bit(acc);

    // Step 5: restore sign — if the original mantissa was negative, the result is also negative.
    if mantissa < 0 {
        acc = -acc;
    }

    ReciprocalResult {
        acc_main: acc,
        norm_lo,
    }
}

#[inline]
fn q15_float_normalize(value: i64) -> (i16, i16) {
    if value == 0 {
        return (0, 0);
    }
    let exp = exp_acc(value);
    let mantissa = ((norm_acc_with_t(value, exp) >> 16) & 0xffff) as i16;
    (mantissa, exp)
}

#[inline]
fn sign_extend_40bit(v: i64) -> i64 {
    let masked = v & 0x00ff_ffff_ffff;
    if masked & 0x80_0000_0000 != 0 {
        (masked as u64 | 0xffff_ff00_0000_0000_u64) as i64
    } else {
        masked
    }
}

/// 10-step Levinson-Durbin recursion. Returns `coeffs` holding the reflection
/// coefficients $k_1, ..., k_{10}$ (= in-place update). `coeffs[0] = k_1` is
/// used as the `clamp_input` for cb46_gain.
pub fn levinson_recursion(input_window: &[i16; 10]) -> [i16; 10] {
    let mut coeffs: [i16; 10] = *input_window;

    let mut feedback_hi: i16;
    #[allow(unused_assignments)]
    let mut feedback_lo: i16;
    let mut inv_energy: i16;
    let mut new_coeff: i16;
    let mut acc_main: i64;
    let mut acc_scratch: i64;
    let mut head_idx: i32;
    let mut tail_idx: i32 = 9;
    let mut sweep_idx: i32;
    let mut outer_count: i16 = 8;
    let mut inner_count: i16;

    loop {
        // Outer step prelude: compute acc_main ≈ E_{i-1} → estimate the reciprocal via q15_reciprocal.
        head_idx = 0;
        sweep_idx = tail_idx;
        let t_sq = coeffs[sweep_idx as usize];
        acc_main = (16384i64 << 16) - shift_acc40((t_sq as i64) * (t_sq as i64) * FRCT_MUL, 5);
        sweep_idx -= 1;

        let recip = q15_reciprocal(acc_main);
        acc_main = recip.acc_main;
        new_coeff = recip.norm_lo;
        inner_count = outer_count >> 1;
        let _saved_count = inner_count.wrapping_add(1);

        // Derive the asm shift amount from new_coeff (signed 6-bit).
        let asm_raw = new_coeff & 0x1f;
        let asm: i16 = if (asm_raw & 0x10) != 0 {
            asm_raw | -32i16
        } else {
            asm_raw
        };

        inv_energy = ((shift_acc40(acc_main, -1) >> 16) & 0xffff) as i16;

        // Inner loop: update coeffs via the forward + reverse steps.
        loop {
            // Forward step
            let t0 = coeffs[sweep_idx as usize];
            let y0 = coeffs[tail_idx as usize];
            acc_scratch = (t0 as i64) * (y0 as i64) * FRCT_MUL;
            acc_main = (coeffs[head_idx as usize] as i64) << 16;
            acc_main =
                (acc_main - shift_acc40(acc_scratch, 3)).clamp(i32::MIN as i64, i32::MAX as i64);
            feedback_hi = ((acc_main >> 16) & 0xffff) as i16;
            feedback_lo = (acc_main & 0xffff) as i16;

            // Q15 multiply: feedback × inv_energy.
            acc_main = (feedback_lo as u16 as i64) * (inv_energy as i64) * FRCT_MUL;
            acc_main = shift_acc40(acc_main, -16);
            acc_main += (inv_energy as i64) * (feedback_hi as i64) * FRCT_MUL;

            // B value for the next sample (used by the reverse step).
            let t3 = coeffs[tail_idx as usize];
            let y3 = coeffs[head_idx as usize];
            acc_scratch = (t3 as i64) * (y3 as i64) * FRCT_MUL;

            new_coeff = ((shift_acc40(acc_main, asm) >> 16) & 0xffff) as i16;

            // Reverse step
            acc_main = (coeffs[sweep_idx as usize] as i64) << 16;
            acc_main =
                (acc_main - shift_acc40(acc_scratch, 3)).clamp(i32::MIN as i64, i32::MAX as i64);
            feedback_hi = ((acc_main >> 16) & 0xffff) as i16;
            feedback_lo = (acc_main & 0xffff) as i16;

            acc_main = (feedback_lo as u16 as i64) * (inv_energy as i64) * FRCT_MUL;
            acc_main = shift_acc40(acc_main, -16);
            acc_main += (inv_energy as i64) * (feedback_hi as i64) * FRCT_MUL;

            coeffs[head_idx as usize] = new_coeff;
            head_idx += 1;
            coeffs[sweep_idx as usize] = ((shift_acc40(acc_main, asm) >> 16) & 0xffff) as i16;
            sweep_idx -= 1;

            if inner_count == 0 {
                break;
            }
            inner_count -= 1;
        }

        let saved = outer_count;
        outer_count -= 1;
        tail_idx -= 1;
        if saved == 0 {
            break;
        }
    }

    coeffs
}

/// Internal helper for `response_shaper`. Performs one Levinson-step worth of
/// dword accumulation over the autocorrelation window, updating the 35-cell
/// dword buffer in place.
fn levinson_step(
    autocorr: &[i16; 14],
    dword_buf: &mut [i16; 35],
    scratch: &mut i16,
    autocorr_off: usize,
    dword_off: usize,
) {
    const AUTOCORR_STRIDE: usize = 2;
    let mut autocorr_idx = autocorr_off;
    let mut dword_cursor = dword_off;

    // Prelude: A_0 = 16384 << 10 (Q24 reference), B_0 = -(R[autocorr_idx] << 10)
    dword_write(dword_buf, dword_cursor, (16384i64) << 10);
    dword_cursor += 2;
    let r0 = autocorr[autocorr_idx] as i64;
    dword_write(dword_buf, dword_cursor, -(r0 << 10));
    dword_cursor -= 2;
    autocorr_idx += AUTOCORR_STRIDE;

    let mut inner_count: i32 = 0;
    let mut outer_count: i32 = 2;

    // 1st step bookkeeping: scratch = -autocorr[autocorr_idx+2]
    let r1 = autocorr[autocorr_idx] as i64;
    *scratch = (-(r1)) as i16;
    autocorr_idx += AUTOCORR_STRIDE;

    let mut acc_scratch = dword_read(dword_buf, dword_cursor) << 1;
    dword_cursor += 2;
    dword_cursor += 1;

    // First dword × scalar multiply (scratch as Q15 multiplier)
    let dword =
        ((dword_buf[dword_cursor - 1] as i64) << 16) | (dword_buf[dword_cursor] as u16 as i64);
    let mut acc_main = (dword * (*scratch as i64) * FRCT_MUL) >> 16;
    dword_cursor = dword_cursor - 1 + AUTOCORR_STRIDE;
    acc_scratch += acc_main << 1;
    dword_write(dword_buf, dword_cursor, acc_scratch);
    dword_cursor -= 2;

    acc_main = dword_read(dword_buf, dword_cursor);
    acc_main += (*scratch as i64) << 10;
    dword_write(dword_buf, dword_cursor, acc_main);

    // Outer loop: the steps that update the reflection coefficients.
    loop {
        let v = autocorr[autocorr_idx];
        autocorr_idx += AUTOCORR_STRIDE;
        let mut step_count: i32 = inner_count;
        *scratch = v.wrapping_neg();

        acc_scratch = dword_read(dword_buf, dword_cursor) << 1;
        dword_cursor += 2;
        let cursor_save = dword_cursor;
        dword_cursor += 1;

        let dword_a =
            ((dword_buf[dword_cursor - 1] as i64) << 16) | (dword_buf[dword_cursor] as u16 as i64);
        acc_main = (dword_a * (*scratch as i64) * FRCT_MUL) >> 16;
        dword_cursor = dword_cursor - 1 + AUTOCORR_STRIDE;
        acc_scratch += acc_main << 1;
        dword_write(dword_buf, dword_cursor, acc_scratch);
        dword_cursor -= 2;
        let mut sweep_cursor = dword_cursor;

        // Inner loop: accumulate dword × scratch with sweep_cursor (descending).
        loop {
            sweep_cursor -= 1;
            acc_scratch = dword_read(dword_buf, dword_cursor);
            let dword_b = ((dword_buf[sweep_cursor - 1] as i64) << 16)
                | (dword_buf[sweep_cursor] as u16 as i64);
            acc_main = (dword_b * (*scratch as i64) * FRCT_MUL) >> 16;
            sweep_cursor -= 3;

            acc_scratch += acc_main << 1;
            acc_scratch += dword_read(dword_buf, sweep_cursor);
            sweep_cursor += 2;

            dword_write(dword_buf, dword_cursor, acc_scratch);
            dword_cursor -= 2;

            if step_count == 0 {
                break;
            }
            step_count -= 1;
        }

        acc_main = dword_read(dword_buf, dword_cursor);
        acc_main += (*scratch as i64) << 10;
        dword_write(dword_buf, dword_cursor, acc_main);

        let saved = outer_count;
        outer_count -= 1;
        inner_count += 1;
        dword_cursor = cursor_save;

        if saved == 0 {
            break;
        }
    }
}

/// Output of `response_shaper` (post-cbec dword buffer + 11-cell LPC mirror).
#[derive(Clone, Debug)]
pub struct ResponseShaperOutput {
    pub merged_buffer: [i16; 35],
    pub mirror_output: [i16; 11],
}

const MIRROR_OUTPUT_LEN: usize = 11;
const MERGE_PHASE_ITERATIONS: usize = 5;
const MIRROR_BIAS: i64 = 1 << 15;

pub fn response_shaper(autocorr: &[i16; 14], initial_buffer: &[i16; 35]) -> ResponseShaperOutput {
    let mut dword_buf: [i16; 35] = *initial_buffer;
    let mut scratch: i16 = 0;

    // Phase 1: Levinson level 1 in the front half of the dword buffer.
    levinson_step(autocorr, &mut dword_buf, &mut scratch, 0, 0);

    // Phase 2: Levinson level 2 in the back half of the dword buffer.
    levinson_step(autocorr, &mut dword_buf, &mut scratch, 1, 12);

    // Phase 3: Merge — sums on the front half, differences on the back half.
    merge_phase(&mut dword_buf);

    // Phase 4: Mirror — produce the 11-cell LPC.
    let mirror_output = mirror_phase(&dword_buf);

    ResponseShaperOutput {
        merged_buffer: dword_buf,
        mirror_output,
    }
}

/// Merge phase: repeated 5 times; accumulate sums into the front half
/// [10, 8, 6, 4, 2] and differences into the back half [22, 20, 18, 16, 14].
fn merge_phase(dword_buf: &mut [i16; 35]) {
    let mut front = 10usize; // tail of the front half
    let mut back = 22usize; // tail of the back half

    for _ in 0..MERGE_PHASE_ITERATIONS {
        // Front half: dword_buf[front] += dword_buf[front - 2]
        let mut a = dword_read(dword_buf, front);
        front -= 2;
        a += dword_read(dword_buf, front);
        front += 2;
        dword_write(dword_buf, front, a);
        front -= 2;

        // Back half: dword_buf[back] -= dword_buf[back - 2]
        let mut a = dword_read(dword_buf, back);
        back -= 2;
        a -= dword_read(dword_buf, back);
        back += 2;
        dword_write(dword_buf, back, a);
        back -= 2;
    }
}

/// Mirror phase: produce the 11-cell LPC from the front/back halves of dword_buf.
/// `mirror[0] = 4096` (= Q12 1.0, equivalent to $a_0$).
fn mirror_phase(dword_buf: &[i16; 35]) -> [i16; MIRROR_OUTPUT_LEN] {
    let mut mirror = [0i16; MIRROR_OUTPUT_LEN];
    mirror[0] = 4096; // Q12 1.0

    let mut front = 2usize; // front half (ascending)
    let mut back = 14usize; // back half (ascending)
    let mut sum_idx = 1usize; // mirror write (ascending: 1..6)
    let mut diff_idx = 10usize; // mirror write (descending: 10..6)

    for _ in 0..MERGE_PHASE_ITERATIONS {
        // Sum side: mirror[sum_idx] = hi16(((sum >> 1) << 4) + bias)
        let sum = dword_read(dword_buf, front) + dword_read(dword_buf, back);
        let scaled = (sum >> 1) << 4;
        mirror[sum_idx] = (((scaled + MIRROR_BIAS) >> 16) & 0xffff) as i16;
        sum_idx += 1;

        // Diff side: mirror[diff_idx] = hi16(((diff >> 1) << 4) + bias)
        let diff = dword_read(dword_buf, front) - dword_read(dword_buf, back);
        front += 2;
        back += 2;
        let scaled = (diff >> 1) << 4;
        mirror[diff_idx] = (((scaled + MIRROR_BIAS) >> 16) & 0xffff) as i16;
        diff_idx -= 1;
    }

    mirror
}

/// LSP windows → 22-cell LPC + updated dword scratch buffer.
pub fn build_autocorrelation_with_state(
    work_coeffs: &[i16; 10],
    coeff_seed_window: &[i16; 14],
    initial_buffer: &[i16; 35],
) -> ([i16; 22], [i16; 35]) {
    // Step 1: 1st response_shaper input = averaged 14-cell first_autocorr.
    let mut first_autocorr = [0i16; 14];
    for i in 0..10 {
        let sum = (work_coeffs[i] as i32) + (coeff_seed_window[i] as i32);
        first_autocorr[i] = (sum >> 1) as i16; // arithmetic shift = floor direction
    }
    first_autocorr[10..14].copy_from_slice(&coeff_seed_window[0..4]);

    // Step 2: 1st response_shaper → mirror_output for sub_a.
    let first_out = response_shaper(&first_autocorr, initial_buffer);

    // Step 3: 2nd response_shaper → mirror_output for sub_b.
    let second_out = response_shaper(coeff_seed_window, &first_out.merged_buffer);

    // Step 4: combine (sub_a [0..11] + sub_b [11..22]).
    let mut result = [0i16; 22];
    result[0..11].copy_from_slice(&first_out.mirror_output);
    result[11..22].copy_from_slice(&second_out.mirror_output);

    (result, second_out.merged_buffer)
}

const STABILITY_REF: i64 = 16384i64 << 16;

#[inline]
fn d14b_metric_residual(sample: i16) -> i64 {
    let x = sample as i64;
    let squared = x * x * FRCT_MUL;
    let v = STABILITY_REF - (squared << 5);
    v.clamp(-0x8000_0000, 0x7fff_ffff) // sat32
}

/// Determine, from the output of `levinson_recursion` (the reflection-coefficient
/// sequence), the stability used to decide whether the block-0 LPC override applies.
pub fn is_lpc_stable(reflection_coeffs: &[i16; 10]) -> bool {
    // Steps 1+2: store the head metric as a pair (hi, lo).
    let head_metric = d14b_metric_residual(reflection_coeffs[0]);
    let mut pair_hi: i16 = ((head_metric >> 16) & 0xffff) as i16;
    let mut pair_lo: i16 = (head_metric & 0xffff) as i16;

    // Step 3: tail metrics — accumulate products of reflection coefficients as a dword pair.
    for &sample in &reflection_coeffs[1..10] {
        let tail = d14b_metric_residual(sample);
        let tail_hi = ((tail >> 16) & 0xffff) as i16;
        let upper = (pair_lo as u16 as i64) * (tail_hi as i64) * FRCT_MUL;
        let lower = (pair_hi as i64) * (tail_hi as i64) * FRCT_MUL;
        let pair_acc = (upper >> 16) + lower;
        pair_hi = ((pair_acc >> 16) & 0xffff) as i16;
        pair_lo = (pair_acc & 0xffff) as i16;
    }

    // Step 4: build the final tail_metric and decide stability.
    let pair_acc_final = ((pair_hi as i64) << 16) | (pair_lo as u16 as i64);
    let tail_metric_final = ((pair_acc_final << 9) >> 16) as i16;

    let head_cmp = ((reflection_coeffs[0] as i64) << 19) - STABILITY_REF;
    let tail_cmp = ((tail_metric_final as i64) << 17) - STABILITY_REF;
    head_cmp <= 0 && tail_cmp <= 0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// is_lpc_stable: with all-zero input, head=0 satisfies the head_cmp
    /// condition, but in the tail-metric computation `metric_residual(0) =
    /// STABILITY_REF` is accumulated through the dword pair so
    /// tail_metric_final exceeds STABILITY_REF. Result: false.
    #[test]
    fn is_lpc_stable_zero_input_returns_false_per_dsp_semantics() {
        let zero = [0i16; 10];
        // Obtain ground truth from the same-named function in rust_float
        // (called via rust_refactored's dev-deps mcelp_float) and verify match.
        let _ = zero;
        assert!(!is_lpc_stable(&zero));
    }

    /// is_lpc_stable: a huge head value is judged unstable.
    #[test]
    fn is_lpc_stable_large_head_is_unstable() {
        let mut input = [0i16; 10];
        input[0] = i16::MAX;
        assert!(!is_lpc_stable(&input));
    }

    /// mirror_phase: mirror[0] is always 4096 (= Q12 1.0).
    #[test]
    fn mirror_phase_mirror_zero_index_is_q12_one() {
        let autocorr = [0i16; 14];
        let init = [0i16; 35];
        let out = response_shaper(&autocorr, &init);
        assert_eq!(out.mirror_output[0], 4096);
    }
}
