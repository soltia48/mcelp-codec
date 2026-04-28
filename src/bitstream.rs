//! Bidirectional conversion between an 18-byte (= 144-bit) MCELP frame and 14 control fields.

use crate::MCELP_FRAME_BYTES;

pub const NUM_FIELDS: usize = 14;
pub const FIELD_WIDTHS: [u8; NUM_FIELDS] = [8, 12, 8, 16, 7, 5, 16, 7, 8, 16, 7, 5, 16, 7];

const SUPPRESS_FLAG_BIT: usize = 138;
const FRAME_FLAG_SRC_BIT: usize = 143;
const FRAME_PAYLOAD_BITS: usize = 139;
const FRAME_RESET_BYTE: usize = MCELP_FRAME_BYTES - 1;
const FRAME_RESET_MASK: u8 = 0x0f;
const FRAME_RESET_CODE: u8 = 0x08;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct ControlFrame {
    pub fields: [u16; NUM_FIELDS],
    pub suppress: bool,
}

impl ControlFrame {
    /// extract 14 fields + suppress from a 144-bit MSB-first stream.
    pub fn unpack(frame: &[u8; MCELP_FRAME_BYTES]) -> Self {
        let mut fields = [0u16; NUM_FIELDS];
        let mut bit_offset = 0usize;
        for (dst, &width) in fields.iter_mut().zip(FIELD_WIDTHS.iter()) {
            *dst = read_be_bits(frame, bit_offset, usize::from(width));
            bit_offset += usize::from(width);
        }
        let suppress = read_be_bits(frame, SUPPRESS_FLAG_BIT, 1) != 0;
        Self { fields, suppress }
    }
}

/// pack 14 fields + suppress into 18 bytes.
pub fn pack_control_frame(control: &ControlFrame) -> [u8; MCELP_FRAME_BYTES] {
    let mut frame = [0u8; MCELP_FRAME_BYTES];
    let mut bit_offset = 0usize;
    for (&value, &width) in control.fields.iter().zip(FIELD_WIDTHS.iter()) {
        write_be_bits(&mut frame, bit_offset, usize::from(width), value);
        bit_offset += usize::from(width);
    }
    if control.suppress {
        write_be_bits(&mut frame, SUPPRESS_FLAG_BIT, 1, 1);
        write_be_bits(&mut frame, FRAME_FLAG_SRC_BIT, 1, 1);
    }
    frame
}

/// canonicalize a received frame (bit 143 → bit 138, clear bits 139..143).
pub fn canonicalize_payload(frame: &[u8; MCELP_FRAME_BYTES]) -> [u8; MCELP_FRAME_BYTES] {
    let mut out = [0u8; MCELP_FRAME_BYTES];
    for n in 0..FRAME_PAYLOAD_BITS {
        if test_bit(frame, n) {
            set_bit(&mut out, n);
        }
    }
    if test_bit(frame, FRAME_FLAG_SRC_BIT) {
        set_bit(&mut out, SUPPRESS_FLAG_BIT);
    }
    out
}

/// `is_reset_frame(frame) = (frame[17] & 0x0F) == 0x08`.
pub fn is_reset_frame(frame: &[u8; MCELP_FRAME_BYTES]) -> bool {
    (frame[FRAME_RESET_BYTE] & FRAME_RESET_MASK) == FRAME_RESET_CODE
}

/// Convert a one-line 36-character hex string of the form "0123abcd..." into 18 bytes.
/// Whitespace is ignored. Returns `None` on failure (odd length / non-hex / wrong length).
pub fn parse_frame_hex(line: &str) -> Option<[u8; MCELP_FRAME_BYTES]> {
    let mut out = [0u8; MCELP_FRAME_BYTES];
    let mut upper: Option<u8> = None;
    let mut written = 0usize;
    for ch in line.bytes() {
        if ch.is_ascii_whitespace() {
            continue;
        }
        let nibble = match ch {
            b'0'..=b'9' => ch - b'0',
            b'a'..=b'f' => 10 + ch - b'a',
            b'A'..=b'F' => 10 + ch - b'A',
            _ => return None,
        };
        if let Some(hi) = upper {
            if written >= MCELP_FRAME_BYTES {
                return None;
            }
            out[written] = (hi << 4) | nibble;
            written += 1;
            upper = None;
        } else {
            upper = Some(nibble);
        }
    }
    if upper.is_some() || written != MCELP_FRAME_BYTES {
        return None;
    }
    Some(out)
}

#[inline]
fn test_bit(buf: &[u8], n: usize) -> bool {
    (buf[n / 8] >> (7 - n % 8)) & 1 != 0
}

#[inline]
fn set_bit(buf: &mut [u8], n: usize) {
    buf[n / 8] |= 1 << (7 - n % 8);
}

#[inline]
fn read_be_bits(buf: &[u8], start_bit: usize, width: usize) -> u16 {
    let mut value = 0u16;
    for k in 0..width {
        value = (value << 1) | u16::from(test_bit(buf, start_bit + k));
    }
    value
}

#[inline]
fn write_be_bits(buf: &mut [u8], start_bit: usize, width: usize, value: u16) {
    for k in 0..width {
        let src_shift = width - 1 - k;
        if (value >> src_shift) & 1 != 0 {
            set_bit(buf, start_bit + k);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_widths_sum_to_138() {
        let total: usize = FIELD_WIDTHS.iter().map(|&w| usize::from(w)).sum();
        assert_eq!(total, 138);
    }

    #[test]
    fn unpack_pack_canon_roundtrip() {
        let frame: [u8; 18] = [
            0xa3, 0x55, 0x12, 0x9f, 0x00, 0xff, 0x80, 0x01, 0x42, 0x7e, 0xcc, 0x33, 0x91, 0x22,
            0x55, 0xaa, 0xf0, 0x00,
        ];
        let canon = canonicalize_payload(&frame);
        let ctrl = ControlFrame::unpack(&canon);
        let repacked = pack_control_frame(&ctrl);
        assert_eq!(canon, repacked);
    }

    #[test]
    fn reset_frame_detection() {
        let mut frame = [0u8; 18];
        frame[17] = 0x08;
        assert!(is_reset_frame(&frame));
        frame[17] = 0x09;
        assert!(!is_reset_frame(&frame));
    }
}
