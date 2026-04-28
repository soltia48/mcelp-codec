//! Gain-quantization pipeline.

use crate::arith::{acc32_16, exp_acc, norm_acc_with_t, sat32, shift_acc40, to_i40};
use crate::tables::gain::{
    COEFF_TABLE, PHASE_TABLE, PREDICT_NORMALIZE_TABLE, PREDICT_SECONDARY_TABLE,
};

#[derive(Clone, Copy, Debug)]
pub struct GainPhaseSetup {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// post-incremented address pointing to the next phase-table cell (= upper half)
    pub phase_first_addr: i16,
    /// post-incremented address pointing to the next phase-table cell (= lower half)
    pub phase_second_addr: i16,
    /// initial value of the gain update base
    pub update_base: i16,
    pub saved_ah: i16,
}

pub fn gain_phase_setup(
    phase_word: i16,
    phase_first_upper: i16,
    phase_first_lower: i16,
    phase_second_upper: i16,
    phase_second_lower: i16,
) -> GainPhaseSetup {
    let phi = phase_word as i64;
    let mut acc_scratch = shift_acc40(phi, -4);
    let mut acc_main: i64 = phi & 0xf;
    acc_scratch = shift_acc40(acc_scratch, 1);
    acc_main = shift_acc40(acc_main, 1);

    acc_scratch += -8248i64;
    let phase_first_addr = (acc_scratch & 0xffff) as u16 as i16;
    acc_main += -8232i64;
    let phase_second_addr = (acc_main & 0xffff) as u16 as i16;

    let first_sum: i64 = ((phase_first_upper as i64) + (phase_first_lower as i64)) << 16;
    let update_base = ((first_sum >> 16) & 0xffff) as i16;
    let second_sum: i64 = ((phase_second_upper as i64) + (phase_second_lower as i64)) << 16;
    let saved_ah = ((second_sum >> 16) & 0xffff) as i16;

    GainPhaseSetup {
        acc_main: second_sum,
        acc_scratch,
        phase_first_addr: phase_first_addr.wrapping_add(1),
        phase_second_addr: phase_second_addr.wrapping_add(1),
        update_base,
        saved_ah,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GainTailDecayInput {
    pub threshold_in: i16,
    pub counter_in: i16,
    /// gain update base before tail-decay scaling (= phase setup output)
    pub update_base: i16,
    /// initial gain before tail-decay scaling (= history setup output)
    pub initial_gain: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct GainTailDecayOutput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    pub threshold_out: i16,
    pub counter_out: i16,
    /// gain update base after tail-decay (or pass-through if no decay)
    pub update_base_out: i16,
    /// initial gain after tail-decay (or pass-through if no decay)
    pub initial_gain_out: i16,
    pub decay_applied: bool,
}

pub fn gain_tail_decay(input: &GainTailDecayInput) -> GainTailDecayOutput {
    let mut acc_main: i64 = input.threshold_in as i64;
    acc_main = to_i40(acc_main - 8); // sub TAIL_THRESHOLD_STEP

    let mut acc_scratch: i64 = input.counter_in as i64;
    if acc_main >= 0 {
        acc_scratch = 4; // TAIL_COUNTER_RELOAD
    }

    let threshold_out: i16 = 0;
    let counter_after_load: i16 = (acc_scratch & 0xffff) as u16 as i16;

    if acc_scratch <= 0 {
        return GainTailDecayOutput {
            acc_main,
            acc_scratch,
            threshold_out,
            counter_out: counter_after_load,
            update_base_out: input.update_base,
            initial_gain_out: input.initial_gain,
            decay_applied: false,
        };
    }

    acc_scratch -= 1;
    acc_main = -acc_scratch;
    acc_main += 3; // TAIL_EXP_BIAS
    let t_value: i16 = (acc_main & 0xffff) as u16 as i16;
    let counter_out: i16 = (acc_scratch & 0xffff) as u16 as i16;

    let mut acc_main: i64 = (3277i64) << 16; // MIN_SYNTH_STATE << 16
    acc_main = shift_acc40(acc_main, t_value);

    let mut acc_scratch: i64 = (input.update_base as i64) * (acc32_16(acc_main) as i64) * 2;
    let update_base_out = ((acc_scratch >> 16) & 0xffff) as i16;
    acc_scratch = (input.initial_gain as i64) * (acc32_16(acc_main) as i64) * 2;
    let initial_gain_out = ((acc_scratch >> 16) & 0xffff) as i16;

    GainTailDecayOutput {
        acc_main,
        acc_scratch,
        threshold_out,
        counter_out,
        update_base_out,
        initial_gain_out,
        decay_applied: true,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GainHistorySetup {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// shift count to apply during history history update (= norm backoff result)
    pub shift_count: i16,
    pub initial_gain: i16,
}

pub fn gain_history_setup(
    predict_acc_main: i64,
    saved_ah: i16,
    sample: i16,
    history_shift: i16,
) -> GainHistorySetup {
    // Restore bits 16..32 of predict_acc_main using saved_ah.
    let restored_bits = (predict_acc_main as u64) & 0x00ff_ffff_ffff_u64;
    let restored_bits = (restored_bits & !(0xffff_u64 << 16)) | ((saved_ah as u16 as u64) << 16);
    let restored_acc = to_i40(restored_bits as i64);

    // norm backoff (no-CMPT path): shift -= 7
    let shift_count_after_backoff = history_shift.wrapping_sub(7);

    // build_initial_gain
    let mut acc_scratch: i64 = (sample as i64) * (acc32_16(restored_acc) as i64) * 2;
    acc_scratch = shift_acc40(acc_scratch, shift_count_after_backoff);
    acc_scratch = shift_acc40(acc_scratch, -1);
    acc_scratch = sat32(acc_scratch);

    let initial_gain = ((acc_scratch >> 16) & 0xffff) as i16;

    GainHistorySetup {
        acc_main: restored_acc,
        acc_scratch,
        shift_count: shift_count_after_backoff,
        initial_gain,
    }
}

pub fn gain_normalize_primary_lookup_index(a_acc_in: i64) -> i16 {
    if a_acc_in <= 0 {
        return 0;
    }
    let mut acc_main = shift_acc40(a_acc_in, 5);
    let t_exp = exp_acc(acc_main);
    acc_main = norm_acc_with_t(acc_main, t_exp);
    acc_main -= 16384i64 << 16;
    acc_main = shift_acc40(acc_main, 5);
    ((shift_acc40(acc_main, -14) >> 16) & 0xffff) as i16
}

pub fn gain_normalize_secondary_lookup_index(a_acc_in: i64) -> i16 {
    let acc_main = shift_acc40(a_acc_in, 5);
    ((shift_acc40(acc_main, -15) >> 16) & 0xffff) as i16
}

#[derive(Clone, Copy, Debug)]
pub struct NormalizePrimaryInput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// initial value for the negated-exponent register (typically 0 from caller)
    pub neg_t_exp_in: i16,
    /// initial value for the lookup-table address pointer (typically 0 from caller)
    pub next_table_addr_in: i16,
    pub table0: i16,
    pub table1: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct NormalizePrimaryOutput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// negated exponent of the normalized accumulator (= -t_exp)
    pub neg_t_exp: i16,
    /// post-incremented address pointing to the next gain_normalize_primary table cell
    pub next_table_addr: i16,
    /// integer index used for `table[idx]` / `table[idx + 1]` lookup
    pub lookup_index: i16,
    /// fractional weight (Q15) used for linear interpolation between the two table cells
    pub interp_frac: i16,
    pub early_return: bool,
}

/// "Normalize correction" leaf helper used within the gain pipeline.
pub fn gain_normalize_primary(input: &NormalizePrimaryInput) -> NormalizePrimaryOutput {
    if input.acc_main <= 0 {
        return NormalizePrimaryOutput {
            acc_main: 0,
            acc_scratch: input.acc_scratch,
            neg_t_exp: input.neg_t_exp_in,
            next_table_addr: input.next_table_addr_in,
            lookup_index: 0,
            interp_frac: 0,
            early_return: true,
        };
    }

    let mut acc_main = shift_acc40(input.acc_main, 5);
    let t_exp = exp_acc(acc_main);
    acc_main = norm_acc_with_t(acc_main, t_exp);

    let mut acc_scratch: i64 = -i64::from(t_exp);
    let neg_t_exp: i16 = (acc_scratch & 0xffff) as u16 as i16;

    acc_main -= 16384i64 << 16;
    acc_main = shift_acc40(acc_main, 5);
    let lookup_index: i16 = ((shift_acc40(acc_main, -14) >> 16) & 0xffff) as i16;

    acc_scratch = (lookup_index as i64) << 16;
    acc_main -= shift_acc40(acc_scratch, 14);
    let interp_frac: i16 = ((shift_acc40(acc_main, -1) >> 16) & 0xffff) as i16;

    let mut next_table_addr = (-5592i16).wrapping_add(lookup_index);
    acc_scratch = (input.table0 as i64) << 16;
    next_table_addr = next_table_addr.wrapping_add(1);

    let y_ext: i64 = i64::from(input.table1);
    acc_main = acc_scratch - (y_ext << 16);
    acc_main = shift_acc40(acc_main, 2);

    acc_scratch -= (interp_frac as i64) * (acc32_16(acc_main) as i64) * 2;
    let hi_word = ((acc_scratch >> 16) & 0xffff) as i16;
    let acc_main = (hi_word as i64) << 16;

    NormalizePrimaryOutput {
        acc_main,
        acc_scratch,
        neg_t_exp,
        next_table_addr,
        lookup_index,
        interp_frac,
        early_return: false,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct NormalizeSecondaryInput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    pub table0: i16,
    pub table1: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct NormalizeSecondaryOutput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// post-incremented address pointing to the next gain_normalize_secondary table cell
    pub next_table_addr: i16,
    /// integer index used for `table[idx]` / `table[idx + 1]` lookup
    pub lookup_index: i16,
    /// fractional weight (Q15) used for linear interpolation between the two table cells
    pub interp_frac: i16,
}

/// Alternative normalization for `gain_normalize_primary` (different coefficients).
/// Called once for the final normalization in predict.
pub fn gain_normalize_secondary(input: &NormalizeSecondaryInput) -> NormalizeSecondaryOutput {
    let mut acc_main = shift_acc40(input.acc_main, 5);
    let lookup_index: i16 = ((shift_acc40(acc_main, -15) >> 16) & 0xffff) as i16;

    let mut acc_scratch: i64 = (lookup_index as i64) << 16;
    acc_main -= shift_acc40(acc_scratch, 15);
    let interp_frac: i16 = ((shift_acc40(acc_main, -2) >> 16) & 0xffff) as i16;

    let mut next_table_addr = (-5625i16).wrapping_add(lookup_index);
    acc_scratch = (input.table0 as i64) << 16;
    next_table_addr = next_table_addr.wrapping_add(1);

    let y_ext: i64 = i64::from(input.table1);
    acc_main = acc_scratch - (y_ext << 16);
    acc_main = shift_acc40(acc_main, 3);
    acc_scratch = shift_acc40(acc_scratch, 1);

    acc_scratch -= (interp_frac as i64) * (acc32_16(acc_main) as i64) * 2;
    acc_scratch = shift_acc40(acc_scratch, -1);

    let hi_word = ((acc_scratch >> 16) & 0xffff) as i16;
    let acc_main = (hi_word as i64) << 16;

    NormalizeSecondaryOutput {
        acc_main,
        acc_scratch,
        next_table_addr,
        lookup_index,
        interp_frac,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GainHistoryInput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// shift count forwarded from history_setup (currently unused inside the body)
    pub shift_count: i16,
    /// negated exponent forwarded from gain_normalize_primary (typically 0 for fresh entry)
    pub neg_t_exp_in: i16,
    pub history: [i16; 3],
    pub normalize_table0: i16,
    pub normalize_table1: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct GainHistoryOutput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    pub history_out: [i16; 4],
    pub scratch_29852: [i16; 4],
}

pub fn gain_history_update(input: &GainHistoryInput) -> GainHistoryOutput {
    let normalize_in = NormalizePrimaryInput {
        acc_main: input.acc_main,
        acc_scratch: input.acc_scratch,
        neg_t_exp_in: input.neg_t_exp_in,
        next_table_addr_in: 0,
        table0: input.normalize_table0,
        table1: input.normalize_table1,
    };
    let normalize_out = gain_normalize_primary(&normalize_in);

    let mut acc_main = normalize_out.acc_main;
    let acc_scratch: i64 = (normalize_out.neg_t_exp as i64) << 16;

    acc_main += shift_acc40(acc_scratch, 15);
    acc_main = shift_acc40(acc_main, -3);
    let saturated_acc = sat32(acc_main);

    let scratch_hi = ((saturated_acc >> 16) & 0xffff) as i16;
    let scratch_lo = (saturated_acc & 0xffff) as i16;
    let scratch_29854: i16 = 24660;

    let mut acc_main: i64 = 0;
    acc_main += (scratch_lo as u16 as i64) * (scratch_29854 as i64) * 2;
    acc_main = shift_acc40(acc_main, -16);
    acc_main += (scratch_hi as i64) * (scratch_29854 as i64) * 2;

    let new_gain = ((acc_main >> 16) & 0xffff) as i16;

    let history_out = [
        new_gain,
        input.history[0],
        input.history[1],
        input.history[2],
    ];

    GainHistoryOutput {
        acc_main,
        acc_scratch,
        history_out,
        scratch_29852: [scratch_hi, scratch_lo, scratch_29854, 0],
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GainPredictiveInput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    pub candidate: [i16; 80],
    pub coeff: [i16; 4],
    pub history: [i16; 5],
    pub normalize_table0: i16,
    pub normalize_table1: i16,
    pub secondary_table0: i16,
    pub secondary_table1: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct GainPredictiveOutput {
    pub acc_main: i64,
    pub acc_scratch: i64,
    /// shift count to feed into the downstream history setup
    pub history_shift_init: i16,
    /// audio sample to feed into the downstream history setup
    pub history_sample_init: i16,
}

pub fn gain_predictive_refresh(input: &GainPredictiveInput) -> GainPredictiveOutput {
    // Phase 1: Energy sum
    let mut acc_main: i64 = (10i64 << 16) << 1;
    for &t in input.candidate.iter() {
        acc_main += (t as i64) * (t as i64) * 2;
    }
    let energy_acc_hi = ((shift_acc40(acc_main, -2) >> 16) & 0xffff) as i16;
    acc_main = (energy_acc_hi as i64) << 16;

    // Phase 2: gain_normalize_primary
    let normalize_in = NormalizePrimaryInput {
        acc_main,
        acc_scratch: input.acc_scratch,
        neg_t_exp_in: 0,
        next_table_addr_in: 0,
        table0: input.normalize_table0,
        table1: input.normalize_table1,
    };
    let normalize_out = gain_normalize_primary(&normalize_in);

    // Phase 3: Linear correction + scratch save
    let mut acc_main = normalize_out.acc_main;
    let scratch_post_normalize: i64 = (normalize_out.neg_t_exp as i64) << 16;
    acc_main += shift_acc40(scratch_post_normalize, 15);
    acc_main = shift_acc40(acc_main, -3);
    let saturated_acc = sat32(acc_main);
    let scratch_29852 = ((saturated_acc >> 16) & 0xffff) as i16;
    let scratch_29853 = (saturated_acc & 0xffff) as i16;
    let const_24660: i16 = 24660;

    // Phase 4: 24660 dword × hi/lo multiply
    let mut acc_main: i64 = 0;
    acc_main += (scratch_29853 as u16 as i64) * (const_24660 as i64) * 2;
    acc_main = shift_acc40(acc_main, -16);
    acc_main += (scratch_29852 as i64) * (const_24660 as i64) * 2;

    // Phase 5: subtract 19488<<19; acc_scratch setup
    acc_main = shift_acc40(acc_main, 3);
    let mut acc_scratch: i64 = (19488i64) << 16;
    acc_main -= shift_acc40(acc_scratch, 3);
    acc_scratch = (18432i64) << 16;
    acc_scratch = shift_acc40(acc_scratch, 4);
    acc_scratch -= acc_main;

    // Phase 6: coeff × history inner product
    let mut acc_main: i64 = (input.coeff[0] as i64) * (input.history[0] as i64) * 2;
    for k in 1..4 {
        acc_main += (input.coeff[k] as i64) * (input.history[k] as i64) * 2;
    }

    // Phase 7: combine + sat32 + dword save
    acc_scratch += shift_acc40(acc_main, 6);
    acc_scratch = shift_acc40(acc_scratch, -5);
    let saturated_scratch = sat32(acc_scratch);
    let scratch_29852_b = ((saturated_scratch >> 16) & 0xffff) as i16;
    let scratch_29853_b = (saturated_scratch & 0xffff) as i16;
    let const_5439: i16 = 5439;

    // Phase 8: 5439 dword × hi/lo multiply → history_shift_init
    let mut acc_main: i64 = 0;
    acc_main += (scratch_29853_b as u16 as i64) * (const_5439 as i64) * 2;
    acc_main = shift_acc40(acc_main, -16);
    acc_main += (scratch_29852_b as i64) * (const_5439 as i64) * 2;
    acc_main = shift_acc40(acc_main, 5);
    let history_shift_init = ((shift_acc40(acc_main, -13) >> 16) & 0xffff) as i16;

    // Phase 9: gain_normalize_secondary
    let b_acc_q15: i64 = (history_shift_init as i64) << 16;
    acc_main -= shift_acc40(b_acc_q15, 13);
    let secondary_input_hi = ((shift_acc40(acc_main, 2) >> 16) & 0xffff) as i16;
    let secondary_acc_in = (secondary_input_hi as i64) << 16;

    let secondary_in = NormalizeSecondaryInput {
        acc_main: secondary_acc_in,
        acc_scratch: b_acc_q15,
        table0: input.secondary_table0,
        table1: input.secondary_table1,
    };
    let secondary_out = gain_normalize_secondary(&secondary_in);

    // Phase 10: history_sample_init = hi16 of post-gain_normalize_secondary acc_main
    let history_sample_init = ((secondary_out.acc_main >> 16) & 0xffff) as i16;

    GainPredictiveOutput {
        acc_main: secondary_out.acc_main,
        acc_scratch: secondary_out.acc_scratch,
        history_shift_init,
        history_sample_init,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GainOrchestrateInput {
    pub phase_word: i16,
    pub threshold_in: i16,
    pub counter_in: i16,
    pub candidate: [i16; 80],
    pub coeff: [i16; 4],
    pub history: [i16; 5],
    pub phase_first_upper: i16,
    pub phase_first_lower: i16,
    pub phase_second_upper: i16,
    pub phase_second_lower: i16,
    pub predict_normalize_table0: i16,
    pub predict_normalize_table1: i16,
    pub predict_secondary_table0: i16,
    pub predict_secondary_table1: i16,
    pub history_normalize_table0: i16,
    pub history_normalize_table1: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct GainOrchestrateOutput {
    pub pitch_gain: i16,
    pub fcb_gain: i16,
    pub history_out: [i16; 5],
    pub threshold_out: i16,
    pub counter_out: i16,
}

/// Orchestrator for the suppress=0 path.
pub fn gain_orchestrate(input: &GainOrchestrateInput) -> GainOrchestrateOutput {
    let phase = gain_phase_setup(
        input.phase_word,
        input.phase_first_upper,
        input.phase_first_lower,
        input.phase_second_upper,
        input.phase_second_lower,
    );

    let predict_in = GainPredictiveInput {
        acc_main: phase.acc_main,
        acc_scratch: phase.acc_scratch,
        candidate: input.candidate,
        coeff: input.coeff,
        history: input.history,
        normalize_table0: input.predict_normalize_table0,
        normalize_table1: input.predict_normalize_table1,
        secondary_table0: input.predict_secondary_table0,
        secondary_table1: input.predict_secondary_table1,
    };
    let predict = gain_predictive_refresh(&predict_in);

    let history_setup = gain_history_setup(
        predict.acc_main,
        phase.saved_ah,
        predict.history_sample_init,
        predict.history_shift_init,
    );

    let history_in = GainHistoryInput {
        acc_main: history_setup.acc_main,
        acc_scratch: history_setup.acc_scratch,
        shift_count: history_setup.shift_count,
        neg_t_exp_in: 0,
        history: [input.history[0], input.history[1], input.history[2]],
        normalize_table0: input.history_normalize_table0,
        normalize_table1: input.history_normalize_table1,
    };
    let history = gain_history_update(&history_in);

    let tail_in = GainTailDecayInput {
        threshold_in: input.threshold_in,
        counter_in: input.counter_in,
        update_base: phase.update_base,
        initial_gain: history_setup.initial_gain,
    };
    let tail = gain_tail_decay(&tail_in);

    let history_out = [
        history.history_out[0],
        input.history[0],
        input.history[1],
        input.history[2],
        input.history[4],
    ];

    GainOrchestrateOutput {
        pitch_gain: tail.update_base_out,
        fcb_gain: tail.initial_gain_out,
        history_out,
        threshold_out: tail.threshold_out,
        counter_out: tail.counter_out,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GainOrchestrateCodecInput {
    pub phase_word: i16,
    pub threshold_in: i16,
    pub counter_in: i16,
    pub candidate: [i16; 80],
    pub history: [i16; 5],
}

pub fn gain_orchestrate_codec(input: &GainOrchestrateCodecInput) -> GainOrchestrateOutput {
    // PHASE_TABLE indices: upper4(phase_word)*2, lower4(phase_word)*2 + 16
    // (phase_word is a 7-bit ctrl field = [0, 127], so idx falls in [0, 14] / [16, 46])
    let phase_first_idx = ((input.phase_word as usize) >> 4) << 1;
    let phase_second_idx = (((input.phase_word as usize) & 0xf) << 1) + 16;

    let phase_first_upper = PHASE_TABLE[phase_first_idx];
    let phase_first_lower = PHASE_TABLE[phase_second_idx];
    let phase_second_upper = PHASE_TABLE[phase_first_idx + 1];
    let phase_second_lower = PHASE_TABLE[phase_second_idx + 1];

    let coeff = COEFF_TABLE;

    // predict path: primary normalize lookup index prediction
    let mut energy_acc: i64 = (10i64 << 16) << 1;
    for &t in input.candidate.iter() {
        energy_acc += (t as i64) * (t as i64) * 2;
    }
    let energy_acc_hi = ((shift_acc40(energy_acc, -2) >> 16) & 0xffff) as i16;
    let predict_normalize_acc_in = (energy_acc_hi as i64) << 16;
    let predict_normalize_idx =
        gain_normalize_primary_lookup_index(predict_normalize_acc_in) as usize;
    let predict_normalize_table0 = PREDICT_NORMALIZE_TABLE
        .get(predict_normalize_idx)
        .copied()
        .unwrap_or(0);
    let predict_normalize_table1 = PREDICT_NORMALIZE_TABLE
        .get(predict_normalize_idx.wrapping_add(1))
        .copied()
        .unwrap_or(0);

    // Phase setup (full)
    let phase = gain_phase_setup(
        input.phase_word,
        phase_first_upper,
        phase_first_lower,
        phase_second_upper,
        phase_second_lower,
    );

    // Replay phases 2-9 of predict to get secondary normalize input acc_main (need partial predict)
    // Replicate predict phases 2-9 here:
    let normalize_in = NormalizePrimaryInput {
        acc_main: predict_normalize_acc_in,
        acc_scratch: phase.acc_scratch,
        neg_t_exp_in: 0,
        next_table_addr_in: 0,
        table0: predict_normalize_table0,
        table1: predict_normalize_table1,
    };
    let normalize_out = gain_normalize_primary(&normalize_in);

    let mut acc_main = normalize_out.acc_main;
    let scratch_post_normalize: i64 = (normalize_out.neg_t_exp as i64) << 16;
    acc_main += shift_acc40(scratch_post_normalize, 15);
    acc_main = shift_acc40(acc_main, -3);
    let saturated_acc = sat32(acc_main);
    let scratch_29852 = ((saturated_acc >> 16) & 0xffff) as i16;
    let scratch_29853 = (saturated_acc & 0xffff) as i16;

    let mut acc_main: i64 = 0;
    acc_main += (scratch_29853 as u16 as i64) * (24660i64) * 2;
    acc_main = shift_acc40(acc_main, -16);
    acc_main += (scratch_29852 as i64) * (24660i64) * 2;
    acc_main = shift_acc40(acc_main, 3);
    let mut acc_scratch: i64 = (19488i64) << 16;
    acc_main -= shift_acc40(acc_scratch, 3);
    acc_scratch = (18432i64) << 16;
    acc_scratch = shift_acc40(acc_scratch, 4);
    acc_scratch -= acc_main;

    let mut acc_main: i64 = (coeff[0] as i64) * (input.history[0] as i64) * 2;
    for k in 1..4 {
        acc_main += (coeff[k] as i64) * (input.history[k] as i64) * 2;
    }

    acc_scratch += shift_acc40(acc_main, 6);
    acc_scratch = shift_acc40(acc_scratch, -5);
    let saturated_scratch = sat32(acc_scratch);
    let scratch_29852_b = ((saturated_scratch >> 16) & 0xffff) as i16;
    let scratch_29853_b = (saturated_scratch & 0xffff) as i16;

    let mut acc_main: i64 = 0;
    acc_main += (scratch_29853_b as u16 as i64) * (5439i64) * 2;
    acc_main = shift_acc40(acc_main, -16);
    acc_main += (scratch_29852_b as i64) * (5439i64) * 2;
    acc_main = shift_acc40(acc_main, 5);
    let history_shift_init = ((shift_acc40(acc_main, -13) >> 16) & 0xffff) as i16;

    let b_acc_q15: i64 = (history_shift_init as i64) << 16;
    acc_main -= shift_acc40(b_acc_q15, 13);
    let secondary_input_hi = ((shift_acc40(acc_main, 2) >> 16) & 0xffff) as i16;
    let secondary_acc_in = (secondary_input_hi as i64) << 16;

    // secondary normalize lookup index prediction
    let predict_secondary_idx = gain_normalize_secondary_lookup_index(secondary_acc_in) as usize;
    let predict_secondary_table0 = PREDICT_SECONDARY_TABLE
        .get(predict_secondary_idx)
        .copied()
        .unwrap_or(0);
    let predict_secondary_table1 = PREDICT_SECONDARY_TABLE
        .get(predict_secondary_idx.wrapping_add(1))
        .copied()
        .unwrap_or(0);

    // Now run predict with predicted tables
    let predict_in = GainPredictiveInput {
        acc_main: phase.acc_main,
        acc_scratch: phase.acc_scratch,
        candidate: input.candidate,
        coeff,
        history: input.history,
        normalize_table0: predict_normalize_table0,
        normalize_table1: predict_normalize_table1,
        secondary_table0: predict_secondary_table0,
        secondary_table1: predict_secondary_table1,
    };
    let predict = gain_predictive_refresh(&predict_in);

    let history_setup = gain_history_setup(
        predict.acc_main,
        phase.saved_ah,
        predict.history_sample_init,
        predict.history_shift_init,
    );

    // history path: primary normalize lookup index prediction
    let history_normalize_idx =
        gain_normalize_primary_lookup_index(history_setup.acc_main) as usize;
    let history_normalize_table0 = PREDICT_NORMALIZE_TABLE
        .get(history_normalize_idx)
        .copied()
        .unwrap_or(0);
    let history_normalize_table1 = PREDICT_NORMALIZE_TABLE
        .get(history_normalize_idx.wrapping_add(1))
        .copied()
        .unwrap_or(0);

    let history_in = GainHistoryInput {
        acc_main: history_setup.acc_main,
        acc_scratch: history_setup.acc_scratch,
        shift_count: history_setup.shift_count,
        neg_t_exp_in: 0,
        history: [input.history[0], input.history[1], input.history[2]],
        normalize_table0: history_normalize_table0,
        normalize_table1: history_normalize_table1,
    };
    let history = gain_history_update(&history_in);

    let tail_in = GainTailDecayInput {
        threshold_in: input.threshold_in,
        counter_in: input.counter_in,
        update_base: phase.update_base,
        initial_gain: history_setup.initial_gain,
    };
    let tail = gain_tail_decay(&tail_in);

    let history_out = [
        history.history_out[0],
        input.history[0],
        input.history[1],
        input.history[2],
        input.history[4],
    ];

    GainOrchestrateOutput {
        pitch_gain: tail.update_base_out,
        fcb_gain: tail.initial_gain_out,
        history_out,
        threshold_out: tail.threshold_out,
        counter_out: tail.counter_out,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gain_phase_setup_zero_input() {
        let s = gain_phase_setup(0, 0, 0, 0, 0);
        assert_eq!(s.acc_main, 0);
        assert_eq!(s.update_base, 0);
        assert_eq!(s.saved_ah, 0);
    }

    #[test]
    fn gain_tail_decay_no_decay_path() {
        let r = gain_tail_decay(&GainTailDecayInput {
            threshold_in: 0,
            counter_in: 0,
            update_base: 1234,
            initial_gain: 5678,
        });
        assert!(!r.decay_applied);
        assert_eq!(r.threshold_out, 0);
        assert_eq!(r.update_base_out, 1234);
        assert_eq!(r.initial_gain_out, 5678);
    }

    #[test]
    fn gain_tail_decay_reload_at_threshold_8() {
        let r = gain_tail_decay(&GainTailDecayInput {
            threshold_in: 8,
            counter_in: 0,
            update_base: 16384,
            initial_gain: 16384,
        });
        // counter reload to 4 → acc_scratch = 4 > 0 → decay path
        // counter -= 1 → counter_out = 3
        assert!(r.decay_applied);
        assert_eq!(r.counter_out, 3);
        assert_eq!(r.threshold_out, 0);
    }

    #[test]
    fn gain_history_setup_basic() {
        let s = gain_history_setup(0, 0, 0, 7);
        // shift_count = 7 - 7 = 0
        assert_eq!(s.shift_count, 0);
        assert_eq!(s.initial_gain, 0);
    }
}
