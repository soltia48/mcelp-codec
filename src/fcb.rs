//! Fixed codebook (FCB) short path (κ=5) and related: dispatch + pulse synthesis + periodic enhancement.

use crate::pitch::decode_lag_fract;
use crate::tables::fcb_main::{
    BIAS, BIT_COUNT, PULSE_DATA, PULSE_POSITION_CODEBOOK, PULSE_SEEDS, T_MPY_TABLE, T_MPYA_TABLE,
};
use crate::tables::fcb_short::{CODEWORD_BIAS, DISPATCH_THRESHOLDS, TRACK_A, TRACK_B};
use crate::tables::pitch::INTERP_FILTER;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FcbDispatch {
    pub lag_class: usize,
    pub clamped_fcb_index: i16,
}

pub fn fcb_dispatch_lag_class(fcb_index: i16) -> FcbDispatch {
    let input_u = fcb_index as u16;
    let clamped_u = if input_u >= DISPATCH_THRESHOLDS[0] {
        DISPATCH_THRESHOLDS[0].wrapping_sub(1)
    } else {
        input_u
    };
    let mut lag_class: i16 = 5;
    for &t in &DISPATCH_THRESHOLDS[1..6] {
        if clamped_u < t {
            lag_class -= 1;
        } else {
            break;
        }
    }
    FcbDispatch {
        lag_class: lag_class as usize,
        clamped_fcb_index: clamped_u as i16,
    }
}

pub fn fcb_short_pulse_synthesis(track_a: usize, track_b: usize) -> [i16; 80] {
    const PULSE_MAG: i16 = 8192;
    let mut buf = [0i16; 80];
    apply_track(&mut buf, &TRACK_A[track_a], PULSE_MAG);
    apply_track(&mut buf, &TRACK_B[track_b], PULSE_MAG);
    buf
}

fn apply_track(buf: &mut [i16; 80], pulses: &[i16; 5], a: i16) {
    for &b in pulses {
        let pos = (b.unsigned_abs() as usize).saturating_sub(1);
        if pos >= 80 {
            continue;
        }
        if b > 0 {
            buf[pos] = buf[pos].saturating_add(a);
        } else if b < 0 {
            buf[pos] = buf[pos].saturating_sub(a);
        }
    }
}

/// Pitch periodicity enhancement IIR.
pub fn pitch_enhance(buf: &mut [i16; 80], lag: i16, sub_lag: i16, gain: i16) {
    let (fract, lag_adjust) = decode_lag_fract(sub_lag);
    let lag_internal = (lag as i32) - 10 + lag_adjust as i32;
    if !(0..80).contains(&lag_internal) {
        return;
    }

    let fwd = INTERP_FILTER[fract];
    let rev: [i16; 10] = match fract {
        0 => {
            let h0 = INTERP_FILTER[0];
            [
                h0[1], h0[2], h0[3], h0[4], h0[5], h0[6], h0[7], h0[8], h0[9], 0,
            ]
        }
        1 => INTERP_FILTER[2],
        2 => INTERP_FILTER[1],
        _ => unreachable!(),
    };

    let copy_count = (lag_internal as usize).min(80);
    let brc = 79 - lag_internal;
    let iter_count = (brc + 1) as usize;
    let scratch_len = 21 + copy_count + iter_count + 16;
    let mut scratch = vec![0i16; scratch_len];
    for i in 0..copy_count {
        scratch[21 + i] = buf[i];
    }

    let mut fwd_center_idx = 11usize;
    let mut buf_idx = copy_count;
    let mut scratch_write_idx = 21 + copy_count;
    let gain_64 = gain as i64;

    for _ in 0..iter_count {
        if buf_idx >= 80 {
            break;
        }
        let mut acc: i64 = 0;
        for k in 0..10 {
            let idx_fwd = (fwd_center_idx as isize) - (k as isize);
            if (0..scratch_len as isize).contains(&idx_fwd) {
                acc += (scratch[idx_fwd as usize] as i64) * (fwd[k] as i64);
            }
            let idx_rev = fwd_center_idx + 1 + k;
            if idx_rev < scratch_len {
                acc += (scratch[idx_rev] as i64) * (rev[k] as i64);
            }
        }
        let acc = acc * 2; // frct_mul
        let scaled = (acc * gain_64) >> 14;
        let combined = scaled + ((buf[buf_idx] as i64) << 16) + 0x8000;
        let hi = (combined >> 16) as i32;
        let new_val = hi.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
        buf[buf_idx] = new_val;
        if scratch_write_idx < scratch_len {
            scratch[scratch_write_idx] = new_val;
        }
        fwd_center_idx += 1;
        buf_idx += 1;
        scratch_write_idx += 1;
    }
}

/// FCB index decoding for the main path (κ < 5).
#[derive(Clone, Copy, Debug)]
pub struct FcbMainPulseDecode {
    pub bit_count: usize,
    pub bit_decomposition: [i16; 3],
    pub recursion_outputs: [i16; 3],
}

pub fn fcb_main_pulse_decode(lag_class: usize, fcb_index: i16) -> FcbMainPulseDecode {
    debug_assert!(lag_class < 5);
    let bit_count = BIT_COUNT[lag_class] as usize;
    let bias = BIAS[lag_class];
    let effective = (fcb_index.wrapping_sub(bias) as u16) as i32;

    let mut bit_decomposition = [0i16; 3];
    for i in 0..bit_count {
        bit_decomposition[i] = ((effective >> i) & 1) as i16;
    }

    let mut a_hi = (effective >> bit_count) as i16;
    let t_mpya = T_MPYA_TABLE[lag_class];
    let t_mpy = T_MPY_TABLE[lag_class];

    let mut recursion_outputs = [0i16; 3];
    for k in 0..bit_count {
        let m1 = ((i32::from(t_mpya[k]) * i32::from(a_hi)) >> 17) as i16;
        recursion_outputs[k] = a_hi.wrapping_sub(t_mpy[k].wrapping_mul(m1));
        a_hi = m1;
    }

    FcbMainPulseDecode {
        bit_count,
        bit_decomposition,
        recursion_outputs,
    }
}

/// pulse template table lookup.
pub fn fcb_main_pulse_template(lag_class: usize, bit_count: usize) -> Vec<i16> {
    debug_assert!(lag_class < 5);
    let mut pulse_idx = (PULSE_SEEDS[lag_class] as usize) * 20;

    let n_blocks = bit_count + 1;
    let mut output = Vec::with_capacity(n_blocks * 80);
    for _ in 0..n_blocks {
        for _ in 0..20 {
            output.push(PULSE_DATA.get(pulse_idx).copied().unwrap_or(0));
            pulse_idx = pulse_idx.wrapping_add(1);
        }
        output.extend(std::iter::repeat_n(0, 60));
    }
    output
}

/// Codebook offset lookup + add/subtract pulses.
pub fn fcb_main_codebook_synth(
    lag_class: usize,
    bit_count: usize,
    bit_decomposition: &[i16],
    recursion_outputs: &[i16],
    pulse_template: &[i16],
) -> [i16; 80] {
    debug_assert!(lag_class < 5);
    debug_assert_eq!(bit_decomposition.len(), bit_count);
    debug_assert_eq!(recursion_outputs.len(), bit_count);
    debug_assert_eq!(pulse_template.len(), bit_count * 80);

    let mut output = [0i16; 80];

    for k in 0..bit_count {
        // recursion_outputs is in descending order; reverse the index for ascending lookup.
        let position = recursion_outputs[bit_count - 1 - k];
        // Base: PULSE_POSITION_CODEBOOK[lag_class * 120 + k * 40 + position]
        let lookup_idx = (lag_class as i32) * 120 + (k as i32) * 40 + (position as i32);
        if lookup_idx < 0 || (lookup_idx as usize) >= PULSE_POSITION_CODEBOOK.len() {
            continue;
        }
        let start = PULSE_POSITION_CODEBOOK[lookup_idx as usize] as i32;
        if !(0..80).contains(&start) {
            continue;
        }
        let start = start as usize;
        let pulse_offset = k * 80;
        let subtract = bit_decomposition[k] == 0;

        for i in 0..(80 - start) {
            let pulse = pulse_template.get(pulse_offset + i).copied().unwrap_or(0) as i32;
            let cur = output[start + i] as i32;
            let new_val = if subtract { cur - pulse } else { cur + pulse };
            output[start + i] = new_val.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
        }
    }

    output
}

/// Main-path orchestration.
pub fn fcb_main_path(
    lag_class: usize,
    clamped_fcb_index: i16,
    lag: i16,
    sub_lag: i16,
    pitch_enhance_gain: i16,
) -> [i16; 80] {
    debug_assert!(lag_class < 5);
    let pulse_decode = fcb_main_pulse_decode(lag_class, clamped_fcb_index);
    let pulse_template = fcb_main_pulse_template(lag_class, pulse_decode.bit_count);
    let mut buf = fcb_main_codebook_synth(
        lag_class,
        pulse_decode.bit_count,
        &pulse_decode.bit_decomposition[..pulse_decode.bit_count],
        &pulse_decode.recursion_outputs[..pulse_decode.bit_count],
        &pulse_template[..pulse_decode.bit_count * 80],
    );
    if lag < 80 {
        pitch_enhance(&mut buf, lag, sub_lag, pitch_enhance_gain);
    }
    buf
}

pub fn fcb_short_path(
    clamped_fcb_index: i16,
    lag: i16,
    sub_lag: i16,
    pitch_enhance_gain: i16,
) -> [i16; 80] {
    let delta = (clamped_fcb_index as u16).wrapping_sub(CODEWORD_BIAS);
    let track_a = (delta >> 7) as usize;
    let track_b = (delta & 0x7F) as usize;
    let mut buf = fcb_short_pulse_synthesis(track_a, track_b);
    if lag < 80 {
        pitch_enhance(&mut buf, lag, sub_lag, pitch_enhance_gain);
    }
    buf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatch_extreme_indices() {
        // 0 (≪ all thresholds) → lag_class = 0
        let r = fcb_dispatch_lag_class(0);
        assert_eq!(r.lag_class, 0);
        // Large value (≫ all thresholds) → lag_class = 5, clamped = T[0] - 1
        let r = fcb_dispatch_lag_class(-1); // = 0xFFFF unsigned
        assert_eq!(r.lag_class, 5);
        assert_eq!(r.clamped_fcb_index, (DISPATCH_THRESHOLDS[0] - 1) as i16);
    }

    #[test]
    fn first_codeword_pulses() {
        // Track A row 0: 0x0001, 0x0011, 0xffd5, 0xffc3, 0x0047
        // = positions 0, 16; 0xffd5(-43)→pos42 neg; 0xffc3(-61)→pos60 neg; pos70 pos
        let buf = fcb_short_pulse_synthesis(0, 0);
        assert_eq!(buf[0], 8192); // pos 0 (= |1|-1)
        assert_eq!(buf[16], 8192); // pos 16 (= |0x11|-1)
        // Track B row 0: 0xfff4(-12)→pos11 neg, 0x0012(18)→pos17 pos, etc.
        // pos 17 becomes +8192 from Track B (Track A[0] has no pulse at pos 17).
        assert_eq!(buf[17], 8192);
    }

    #[test]
    fn short_path_skips_pitch_enhance_gain_for_large_lag() {
        // Use input that hits the short path (= small delta from bias).
        let clamped = (CODEWORD_BIAS as i16).wrapping_add(100);
        let with_skip = fcb_short_path(clamped, 100, 0, 16384);
        // Reproduce track_a/b from the same input and compare.
        let delta = (clamped as u16).wrapping_sub(CODEWORD_BIAS);
        let track_a = (delta >> 7) as usize;
        let track_b = (delta & 0x7F) as usize;
        let pure_pulse = fcb_short_pulse_synthesis(track_a, track_b);
        assert_eq!(with_skip, pure_pulse);
    }

    #[test]
    fn zero_gain_preserves_buffer() {
        let mut buf = [100i16; 80];
        let original = buf;
        pitch_enhance(&mut buf, 30, 0, 0);
        // gain=0, so scaled=0; combined = (buf << 16) + 0x8000 → hi16 = buf.
        // (The rounding bias must not add +1: hi16 of buf << 16 + 0x8000 equals buf.)
        for i in 30 - 10..80 {
            // Only this range can be modified.
            assert_eq!(buf[i], original[i], "buf[{i}]");
        }
    }
}
