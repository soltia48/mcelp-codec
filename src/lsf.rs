//! LSF (Line Spectral Frequency) decoding.

use crate::LPC_ORDER;
use crate::arith::sat16;
use crate::bitstream::ControlFrame;
use crate::tables::lsf::{INIT_LSF_TEMPLATE, PREDICTOR, SMOOTH, STAGE1_CODEBOOK, STAGE2_CODEBOOK};
use crate::tables::lsp::{COS_LUT_SLOPE, COS_LUT_VALUE};

/// Four indices extracted from $F_0$ / $F_1$.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LsfIndices {
    pub mode: bool,
    pub seed: u8,
    pub upper: u8,
    pub lower: u8,
}

impl LsfIndices {
    pub fn from_control(c: &ControlFrame) -> Self {
        let f0 = c.fields[0];
        let f1 = c.fields[1];
        Self {
            mode: (f0 & 0x80) != 0,
            seed: (f0 & 0x7f) as u8,
            upper: ((f1 >> 6) & 0x3f) as u8,
            lower: (f1 & 0x3f) as u8,
        }
    }
}

/// Rolling history holding the scratch for the past 3 frames.
/// `slots[0]` = most recent, `slots[2]` = oldest. `new()` initializes all 3
/// slots with `INIT_LSF_TEMPLATE`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LsfHistory {
    pub slots: [[i16; LPC_ORDER]; 3],
}

impl Default for LsfHistory {
    fn default() -> Self {
        Self::new()
    }
}

impl LsfHistory {
    pub fn new() -> Self {
        Self {
            slots: [INIT_LSF_TEMPLATE, INIT_LSF_TEMPLATE, INIT_LSF_TEMPLATE],
        }
    }

    pub fn reset(&mut self) {
        *self = Self::new();
    }

    /// Shift slots back and insert the new vector at slot 0.
    pub fn splice_new_vector(&mut self, new_vector: &[i16; LPC_ORDER]) {
        self.slots[2] = self.slots[1];
        self.slots[1] = self.slots[0];
        self.slots[0] = *new_vector;
    }
}

/// Split-VQ composition of stage1[seed] + stage2[upper/lower].
pub fn combine_codebook(idx: &LsfIndices) -> [i16; LPC_ORDER] {
    debug_assert!((idx.seed as usize) < STAGE1_CODEBOOK.len());
    debug_assert!((idx.upper as usize) < STAGE2_CODEBOOK.len());
    debug_assert!((idx.lower as usize) < STAGE2_CODEBOOK.len());

    let s1 = &STAGE1_CODEBOOK[idx.seed as usize];
    let s2_first = &STAGE2_CODEBOOK[idx.upper as usize];
    let s2_second = &STAGE2_CODEBOOK[idx.lower as usize];

    let mut out = [0i16; LPC_ORDER];
    for k in 0..5 {
        out[k] = s1[k].wrapping_add(s2_first[k]);
    }
    for k in 5..LPC_ORDER {
        out[k] = s1[k].wrapping_add(s2_second[k]);
    }
    out
}

/// Smooth adjacent pairs.
/// The caller passes 10 (first call) or 5 (second call) as `min_gap`.
pub fn enforce_min_gap(lsf: &mut [i16; LPC_ORDER], min_gap: i16) {
    for i in 0..LPC_ORDER - 1 {
        let left = lsf[i] as i32;
        let right = lsf[i + 1] as i32;
        let diff = left - right + min_gap as i32;
        if diff > 0 {
            let delta_right = diff >> 1;
            let delta_left = (diff + 1) >> 1;
            lsf[i] = (left - delta_left) as i16;
            lsf[i + 1] = (right + delta_right) as i16;
        }
    }
}

/// Linear-predictive combination of scratch + 3 frames of history.
pub fn predictive_combine(
    scratch: &[i16; LPC_ORDER],
    history: &LsfHistory,
    mode: bool,
) -> [i16; LPC_ORDER] {
    let mode_idx = usize::from(mode);
    let pred = &PREDICTOR[mode_idx];
    let smooth = &SMOOTH[mode_idx];

    let mut out = [0i16; LPC_ORDER];
    for i in 0..LPC_ORDER {
        let mut acc: i64 = (scratch[i] as i64) * (smooth[i] as i64) * 2;
        acc += (history.slots[0][i] as i64) * (pred[i] as i64) * 2;
        acc += (history.slots[1][i] as i64) * (pred[i + 10] as i64) * 2;
        acc += (history.slots[2][i] as i64) * (pred[i + 20] as i64) * 2;
        out[i] = ((acc >> 16) & 0xffff) as i16;
    }
    out
}

/// Four-stage alignment (bubble + clamp + min gap).
pub fn stabilize_and_finalize(lsf: &mut [i16; LPC_ORDER]) {
    // Phase 1: bubble 1 pass
    for i in 0..LPC_ORDER - 1 {
        if lsf[i + 1] < lsf[i] {
            lsf.swap(i, i + 1);
        }
    }
    // Phase 2: lower clamp (41 / 32768 ≈ 0.00125)
    if lsf[0] < 41 {
        lsf[0] = 41;
    }
    // Phase 3: forward min-gap = 321
    for i in 0..LPC_ORDER - 1 {
        let floor = (lsf[i] as i32).wrapping_add(321) as i16;
        if lsf[i + 1] < floor {
            lsf[i + 1] = floor;
        }
    }
    // Phase 4: upper clamp (25682 / 32768 ≈ 0.7838)
    if lsf[LPC_ORDER - 1] > 25682 {
        lsf[LPC_ORDER - 1] = 25682;
    }
}

/// Full LSF decoding pipeline.
pub fn decode_lsf_direct_mode(idx: &LsfIndices, history: &mut LsfHistory) -> [i16; LPC_ORDER] {
    let mut scratch = combine_codebook(idx);
    enforce_min_gap(&mut scratch, 10);
    enforce_min_gap(&mut scratch, 5);
    let mut output = predictive_combine(&scratch, history, idx.mode);
    history.splice_new_vector(&scratch);
    stabilize_and_finalize(&mut output);
    output
}

/// LSF (Q15) → LSP (Q15) cosine transform.
pub fn lsf_to_lsp(lsf: &[i16; LPC_ORDER]) -> [i16; LPC_ORDER] {
    let mut out = [0i16; LPC_ORDER];
    for k in 0..LPC_ORDER {
        // scaled_index = lsf[k] * 20861 * 2 (i64), Q15 * 17-bit factor * frct_mul
        let scaled = (lsf[k] as i64) * 20861 * 2;
        // index = clamp((scaled >> 24), 0, 63), fraction_bits = (scaled >> 16) & 0xff
        let index = ((scaled >> 24) as i32).clamp(0, 63) as usize;
        let fraction_bits = ((scaled >> 16) & 0xff) as i32;
        let value = COS_LUT_VALUE[index] as i32;
        let slope = COS_LUT_SLOPE[index] as i32;
        // out = value + slope * fraction_bits / 4096
        let interp = (slope * fraction_bits) >> 12;
        out[k] = value.wrapping_add(interp) as i16;
    }
    out
}

/// Element-wise average of two 10-word vectors.
/// `out[i] = (prev[i] + curr[i]) >> 1` (= round toward floor).
pub fn block_average_lsp(prev: &[i16; LPC_ORDER], curr: &[i16; LPC_ORDER]) -> [i16; LPC_ORDER] {
    std::array::from_fn(|i| {
        let sum = (prev[i] as i32) + (curr[i] as i32);
        (sum >> 1) as i16
    })
}

/// Derive the LPC coefficients for the 2 subframes within a block.
/// Returns: `(sub_a, sub_b)`. Side effect: `prev_block_coeffs ← block_coeffs`.
pub fn subframe_lsp_interp(
    prev_block_coeffs: &mut [i16; LPC_ORDER],
    block_coeffs: &[i16; LPC_ORDER],
) -> ([i16; LPC_ORDER], [i16; LPC_ORDER]) {
    let sub_a = block_average_lsp(prev_block_coeffs, block_coeffs);
    let sub_b = *block_coeffs;
    *prev_block_coeffs = *block_coeffs;
    (sub_a, sub_b)
}

/// Convert LSP (Q15) to LPC (Q12) (P/Q polynomial decomposition).
///
/// $A(z) = \frac{1}{2}\left((1+z^{-1})F_1(z) + (1-z^{-1})F_2(z)\right)$
/// where $F_1$ uses even-indexed LSPs, $F_2$ uses odd-indexed LSPs.
///
/// Output: `[i16; 11]` Q12, with `a[0] = 4096` (= 1.0 in Q12).
pub fn lsp_to_lpc(lsp: &[i16; LPC_ORDER]) -> [i16; 11] {
    let f1 = poly_from_5_lsps([lsp[0], lsp[2], lsp[4], lsp[6], lsp[8]]);
    let f2 = poly_from_5_lsps([lsp[1], lsp[3], lsp[5], lsp[7], lsp[9]]);

    let mut a = [0i16; 11];
    a[0] = 1 << 12; // Q12 1.0
    for m in 1..=10 {
        // p_tilde[m] = F1[m] + F1[m-1], q_tilde[m] = F2[m] - F2[m-1]
        let f1_m = f1[m];
        let f1_prev = f1[m - 1];
        let f2_m = f2[m];
        let f2_prev = f2[m - 1];
        let sum_q14 = (f1_m + f1_prev) + (f2_m - f2_prev); // Q14
        // a[m] = sum_q14 / 2 in Q14, then Q14 → Q12 (>> 2). Combined: (sum + 4) >> 3.
        let q12 = (sum_q14 + 4) >> 3;
        a[m] = sat16(q12);
    }
    a
}

/// Extract the **first reflection coefficient $k_1$** from LPC (Q12, length 11) via inverse Levinson.
///
/// Output is Q15. Used as `clamp_input` for `cb46_gain`.
/// Returns 0 if the LPC is unstable (`|k_step| >= 1`) or denom becomes ≤ 0.
///
/// Algorithm: Schur recursion (= inverse Levinson)
/// ```text
/// for step = p..=2:
///     k_step = a[step]
///     if |k_step| >= 1.0: stop (unstable)
///     denom = 1 - k_step^2
///     for m = 1..step:
///         a'[m] = (a[m] - k_step * a[step - m]) / denom
///     a = a'
/// k_1 = a[1]
/// ```
///
/// Q-format: computed internally in Q31 (i64); output is Q15 (i16). LPC input is Q12.
pub fn lpc_to_first_reflection(lpc_q12: &[i16; 11]) -> i16 {
    let one_q31: i64 = 1i64 << 31;
    let mut a: [i64; 11] = std::array::from_fn(|i| (lpc_q12[i] as i64) << 19); // Q12 → Q31
    a[0] = one_q31; // a_0 = 1.0 in Q31

    for step in (2..=10).rev() {
        let k = a[step];
        if k >= one_q31 || k <= -one_q31 {
            return 0; // unstable
        }
        // k^2 in Q31: (k * k) >> 31
        let k_sq = ((k as i128 * k as i128) >> 31) as i64;
        let denom = one_q31 - k_sq;
        if denom <= 0 {
            return 0;
        }
        let mut new_a = [0i64; 11];
        for m in 1..step {
            // (a[m] - k * a[step - m]) / denom, all in Q31
            let kx = ((k as i128 * a[step - m] as i128) >> 31) as i64;
            let num = a[m] - kx;
            // divide num by denom: result in Q31
            // (num as i128 << 31) / (denom as i128) = Q62 / Q31 = Q31
            let result = ((num as i128) << 31) / (denom as i128);
            // saturate
            new_a[m] = result.clamp(i64::MIN as i128, i64::MAX as i128) as i64;
        }
        a[1..step].copy_from_slice(&new_a[1..step]);
    }
    // k_1 = a[1] in Q31, convert to Q15
    let k1_q15 = a[1] >> 16; // Q31 >> 16 = Q15
    k1_q15.clamp(i16::MIN as i64, i16::MAX as i64) as i16
}

/// From 5 LSPs (Q15), compute the full coefficients [i32; 11] of
/// $F(z) = \prod_{k=0}^{4}(1 - 2 \omega_k z^{-1} + z^{-2})$
/// (= an 11th-order palindromic polynomial) in Q14.
fn poly_from_5_lsps(lsps: [i16; 5]) -> [i32; 11] {
    let mut p = [0i64; 11];
    p[0] = 1 << 14; // Q14 1.0
    let mut len: usize = 1;
    for &omega in &lsps {
        let mut next = [0i64; 11];
        for j in 0..(len + 2) {
            let mut acc: i64 = 0;
            if j < len {
                acc += p[j];
            }
            if j >= 1 && (j - 1) < len {
                // -2ω · p[j-1] in Q14: ((omega_q15 * p_q14) >> 14)
                let prod = (omega as i64) * p[j - 1];
                acc -= prod >> 14;
            }
            if j >= 2 && (j - 2) < len {
                acc += p[j - 2];
            }
            next[j] = acc;
        }
        p = next;
        len += 2;
    }
    let mut out = [0i32; 11];
    for i in 0..11 {
        out[i] = p[i] as i32; // Q14; bounds are guaranteed by LSP stability
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lsf_indices_from_control_extracts_4_fields() {
        let mut c = ControlFrame::default();
        c.fields[0] = 0b1_0010_101;
        c.fields[1] = 0b1011_01_001110;
        let idx = LsfIndices::from_control(&c);
        assert!(idx.mode);
        assert_eq!(idx.seed, 21);
        assert_eq!(idx.upper, 45);
        assert_eq!(idx.lower, 14);
    }

    #[test]
    fn history_init_uses_template() {
        let h = LsfHistory::new();
        for slot in &h.slots {
            assert_eq!(slot, &INIT_LSF_TEMPLATE);
        }
    }

    #[test]
    fn history_splice_shifts_slots() {
        let mut h = LsfHistory::new();
        let v0 = [1i16; LPC_ORDER];
        let v1 = [2i16; LPC_ORDER];
        h.splice_new_vector(&v0);
        h.splice_new_vector(&v1);
        assert_eq!(h.slots[0], v1);
        assert_eq!(h.slots[1], v0);
        assert_eq!(h.slots[2], INIT_LSF_TEMPLATE);
    }

    #[test]
    fn combine_codebook_rows_monotone_for_sample_indices() {
        // By table design, stage1[seed] + stage2[*] is monotonically increasing.
        for &(seed, upper, lower) in &[(0, 0, 0), (10, 7, 31), (127, 63, 63)] {
            let idx = LsfIndices {
                mode: false,
                seed,
                upper,
                lower,
            };
            let v = combine_codebook(&idx);
            for i in 1..LPC_ORDER {
                assert!(
                    v[i] > v[i - 1],
                    "non-monotone seed={seed} u={upper} l={lower}"
                );
            }
        }
    }

    #[test]
    fn stabilize_and_finalize_satisfies_invariants() {
        // Even with initially narrow pairs, after stabilize all are separated by ≥ 321.
        let mut x = [40, 100, 200, 300, 400, 500, 600, 700, 800, 26000];
        stabilize_and_finalize(&mut x);
        assert!(x[0] >= 41);
        assert!(x[LPC_ORDER - 1] <= 25682);
        for i in 1..LPC_ORDER {
            assert!(x[i] - x[i - 1] >= 321, "gap violation at {i}");
        }
    }
}
