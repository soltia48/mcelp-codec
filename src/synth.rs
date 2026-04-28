//! Excitation synthesis + LPC synthesis filter + control.

const SYNTH_CONTROL_ENABLE: i16 = 16384;
const SYNTH_CONTROL_DISABLE: i16 = 0;
const SYNTH_CONTROL_GATE: i16 = 13107;
const PREV_LAG_BLEND_GAIN: i16 = 24576;
const MIN_SYNTH_STATE: i16 = 3277;
const MAX_SYNTH_STATE: i16 = 13017;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SynthControl {
    pub base: i16,
    pub control: i16,
}

/// Main computation logic for gain.
///
/// Updates `control` via four conditional branches:
/// 1. `clamp_input` too large → DISABLE (= 0)
/// 2. `prev_control_state` smaller than gate → unchanged
/// 3. `selected_lag` ≤ `prev_lag * 24576/32768` → unchanged
/// 4. `selected_lag` < `prev_lag * 24576/32768 * 2` → ENABLE, otherwise unchanged
pub fn update_synth_control(
    selected_lag: i16,
    prev_lag: i16,
    clamp_input: i16,
    prev_control_state: i16,
) -> SynthControl {
    let base = SYNTH_CONTROL_ENABLE;

    let clamp_acc: i64 = ((clamp_input as i64) << 19) - ((SYNTH_CONTROL_ENABLE as i64) << 16);
    if clamp_acc > 0 {
        return SynthControl {
            base,
            control: SYNTH_CONTROL_DISABLE,
        };
    }

    let control_acc: i64 =
        ((prev_control_state as i64) << 17) - ((SYNTH_CONTROL_GATE as i64) << 16);
    if control_acc <= 0 {
        return SynthControl {
            base,
            control: prev_control_state,
        };
    }

    let lag_acc: i64 = (selected_lag as i64) << 16;
    let lag_delta: i64 = lag_acc - (PREV_LAG_BLEND_GAIN as i64) * (prev_lag as i64) * 2;
    if lag_delta <= 0 {
        return SynthControl {
            base,
            control: prev_control_state,
        };
    }

    let prev_lag_mix: i64 = (prev_lag as i64) * (PREV_LAG_BLEND_GAIN as i64) * 2;
    let remaining: i64 = lag_acc - (prev_lag_mix << 1);
    let control = if remaining < 0 {
        SYNTH_CONTROL_ENABLE
    } else {
        prev_control_state
    };

    SynthControl { base, control }
}

/// Combine fcb_dispatch + update_synth_control to derive the `gain` value
/// passed to pitch_enhance.
///
/// - `lag_class ∈ {0, 3, 4}` → always ENABLE (16384)
/// - `lag_class ∈ {1, 2, 5}` → `control` from update_synth_control
pub fn compute_pitch_enhance_gain(
    fcb_index: i16,
    current_lag: i16,
    clamp_input: i16,
    prev_lag: i16,
    prev_pitch_gain: i16,
) -> i16 {
    let dispatch = crate::fcb::fcb_dispatch_lag_class(fcb_index);
    let prev_control_state = prev_pitch_gain.clamp(MIN_SYNTH_STATE, MAX_SYNTH_STATE);
    let synth_control =
        update_synth_control(current_lag, prev_lag, clamp_input, prev_control_state);
    match dispatch.lag_class {
        0 | 3 | 4 => SYNTH_CONTROL_ENABLE,
        _ => synth_control.control,
    }
}

/// Excitation = 8·g_c·c[i] + 4·g_p·v[i] + round.
pub fn mix_excitation(
    adaptive_codebook: &[i16; 80],
    fixed_codebook: &[i16; 80],
    pitch_gain: i16,
    fcb_gain: i16,
) -> [i16; 80] {
    let mut output = [0i16; 80];
    for i in 0..80 {
        let mut acc: i64 = (fcb_gain as i64) * (fixed_codebook[i] as i64) * 8
            + (pitch_gain as i64) * (adaptive_codebook[i] as i64) * 4
            + 32768;
        // i32 saturation
        if acc > i32::MAX as i64 {
            acc = i32::MAX as i64;
        }
        if acc < i32::MIN as i64 {
            acc = i32::MIN as i64;
        }
        let hi = (acc >> 16) as i32;
        output[i] = hi.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
    }
    output
}

/// 1/A(z) all-pole IIR; synthesizes 80 samples and produces the history for the next sub.
///
/// `lpc_coeffs[0..10]` = $a_1..a_{10}$ (Q12); `history[0..10]` = tail of the previous sub's output.
pub fn lpc_synthesis_filter(
    excitation: &[i16; 80],
    lpc_coeffs: &[i16; 10],
    history: &[i16; 10],
) -> ([i16; 80], [i16; 10]) {
    let mut buf = [0i16; 90];
    buf[0..10].copy_from_slice(history);

    for n in 0..80 {
        let mut mac: i64 = 0;
        for k in 0..10 {
            mac += (lpc_coeffs[k] as i64) * (buf[9 + n - k] as i64); // Q12 × Q15 = Q27
        }
        // acc = excitation << 16 (Q31) - mac * 16 (Q27 → Q31) + round
        let acc: i64 = (excitation[n] as i64) * 65536 - mac * 16 + 32768;
        // i32 saturation → hi16 → i16 saturation
        let sat32 = acc.clamp(i32::MIN as i64, i32::MAX as i64);
        let hi = (sat32 >> 16) as i32;
        buf[10 + n] = hi.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
    }

    let mut output = [0i16; 80];
    output.copy_from_slice(&buf[10..90]);
    let mut new_history = [0i16; 10];
    new_history.copy_from_slice(&buf[80..90]);
    (output, new_history)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mix_excitation_zero_inputs_zero_output() {
        let zero = [0i16; 80];
        let out = mix_excitation(&zero, &zero, 8000, 8000);
        assert!(out.iter().all(|&v| v == 0));
    }

    #[test]
    fn mix_excitation_pitch_only_scales_correctly() {
        let pitch_codebook = [100i16; 80];
        let fcb_codebook = [0i16; 80];
        // g_p = Q15 1.0 = 32767; extract hi16 of 4×: 100 * 32767 * 4 / 65536 ≈ 200
        let out = mix_excitation(&pitch_codebook, &fcb_codebook, i16::MAX, 0);
        assert!(out.iter().all(|&v| v >= 199 && v <= 201));
    }

    #[test]
    fn lpc_synthesis_filter_zero_excitation_zero_lpc_zero_output() {
        let exc = [0i16; 80];
        let lpc = [0i16; 10];
        let hist = [0i16; 10];
        let (out, new_hist) = lpc_synthesis_filter(&exc, &lpc, &hist);
        assert!(out.iter().all(|&v| v == 0));
        assert!(new_hist.iter().all(|&v| v == 0));
    }

    #[test]
    fn lpc_synthesis_filter_impulse_response_with_zero_lpc() {
        // LPC all zeros (= pure passthrough); excitation = impulse → output = impulse
        let mut exc = [0i16; 80];
        exc[0] = 16384; // Q15 = 0.5
        let lpc = [0i16; 10];
        let hist = [0i16; 10];
        let (out, _) = lpc_synthesis_filter(&exc, &lpc, &hist);
        assert_eq!(out[0], 16384);
        assert!(out[1..].iter().all(|&v| v == 0));
    }
}
