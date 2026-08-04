//! Transport frame handling: 18 bytes in, 14 parameter fields out.
//!
//! A frame carries 144 bits.  139 of them are payload, laid out MSB first:
//! fourteen parameter fields totalling 138 bits followed by a single
//! frame-suppression bit.  The remaining bits are transport padding, except
//! that bit 143 mirrors the suppression bit and the top nibble of the last byte
//! doubles as an in-band reset marker.

/// Bytes per transport frame.
pub const FRAME_BYTES: usize = 18;
/// Parameter fields per frame.
pub const FIELDS: usize = 14;
/// Payload bits, i.e. the fourteen fields plus the suppression bit.
pub const PAYLOAD_BITS: usize = 139;
/// Bit position of the suppression flag inside the payload.
pub const SUPPRESS_BIT: usize = 138;
/// Bit position of the transmitted copy of the suppression flag.
pub const SUPPRESS_SRC_BIT: usize = 143;

/// The fourteen parameter fields of one frame.
///
/// | field | bits | meaning                                        |
/// |-------|------|------------------------------------------------|
/// | 0     | 8    | LSF predictor mode (bit 7) + first stage index |
/// | 1     | 12   | LSF second stage, two 6-bit halves             |
/// | 2..4  | 8/16/7  | subframe 0: pitch lag, code index, gain     |
/// | 5..7  | 5/16/7  | subframe 1: lag delta, code index, gain     |
/// | 8..10 | 8/16/7  | subframe 2: pitch lag, code index, gain     |
/// | 11..13| 5/16/7  | subframe 3: lag delta, code index, gain     |
#[derive(Clone, Copy, Debug, Default)]
pub struct Params {
    /// Raw field values, in field order.
    pub field: [i16; FIELDS],
    /// Frame-suppression / erasure flag.
    pub suppress: bool,
}

/// Fields one subframe carries, and their order within its group.
pub const SUBFRAME_FIELDS: usize = 3;
pub const LAG: usize = 0;
pub const CODE: usize = 1;
pub const GAIN: usize = 2;

/// The three fields of one subframe, read out by name.
#[derive(Clone, Copy, Debug)]
pub struct Subframe {
    /// Pitch lag: absolute in the first subframe of a half-frame, relative to
    /// it in the second.
    pub lag: i16,
    /// Fixed-codebook index.
    pub code: i16,
    /// Joint gain index.
    pub gain: i16,
}

impl Params {
    /// Index of the first field belonging to half-frame `half` (0 or 1).
    pub fn half_base(half: usize) -> usize {
        2 + 6 * half
    }

    /// Index of the first of the three fields subframe `sub` of half-frame
    /// `half` carries.
    pub fn subframe_base(half: usize, sub: usize) -> usize {
        Self::half_base(half) + SUBFRAME_FIELDS * sub
    }

    /// The three fields of one subframe.
    pub fn subframe(&self, half: usize, sub: usize) -> Subframe {
        let base = Self::subframe_base(half, sub);
        Subframe {
            lag: self.field[base + LAG],
            code: self.field[base + CODE],
            gain: self.field[base + GAIN],
        }
    }
}

/// Split the nine packed words of a frame into the fourteen parameter fields.
///
/// A straight-line bit unpacker walking the payload MSB
/// first with the field widths `[8,12,8,16,7,5,16,7,8,16,7,5,16,7]`.
pub fn unpack(words: &[i16; 9]) -> Params {
    let w = |i: usize| words[i] as u16 as u32;
    let mut p = Params::default();

    p.field[0] = (w(0) >> 8) as i16;
    p.field[1] = (((w(0) & 0xff) << 4) | (w(1) >> 12)) as i16;
    p.field[2] = ((w(1) >> 4) & 0xff) as i16;
    p.field[3] = (((w(1) & 0x0f) << 12) | (w(2) >> 4)) as i16;
    p.field[4] = (((w(2) & 0x0f) << 3) | (w(3) >> 13)) as i16;
    p.field[5] = ((w(3) >> 8) & 0x1f) as i16;
    p.field[6] = (((w(3) & 0xff) << 8) | (w(4) >> 8)) as i16;
    p.field[7] = ((w(4) >> 1) & 0x7f) as i16;
    p.field[8] = (((w(4) & 0x01) << 7) | (w(5) >> 9)) as i16;
    p.field[9] = (((w(5) & 0x1ff) << 7) | (w(6) >> 9)) as i16;
    p.field[10] = ((w(6) >> 2) & 0x7f) as i16;
    p.field[11] = (((w(6) & 0x03) << 3) | (w(7) >> 13)) as i16;
    p.field[12] = (((w(7) & 0x1fff) << 3) | (w(8) >> 13)) as i16;
    p.field[13] = ((w(8) >> 6) & 0x7f) as i16;
    p.suppress = (w(8) >> 5) & 1 != 0;
    p
}

/// Copy the payload bits out of a received frame, dropping transport padding
/// and folding bit 143 down onto the suppression bit.
///
/// Mirrors the bit shuffling `main()` performs before handing a frame to the
/// decoder.
pub fn canonicalize(frame: &[u8; FRAME_BYTES]) -> [i16; 9] {
    let mut payload = [0u8; FRAME_BYTES];
    for bit in 0..PAYLOAD_BITS {
        if test_bit(frame, bit) {
            set_bit(&mut payload, bit);
        }
    }
    if test_bit(frame, SUPPRESS_SRC_BIT) {
        set_bit(&mut payload, SUPPRESS_BIT);
    }
    let mut words = [0i16; 9];
    for (i, w) in words.iter_mut().enumerate() {
        *w = (((payload[i * 2] as u16) << 8) | payload[i * 2 + 1] as u16) as i16;
    }
    words
}

/// True when a frame carries the in-band "reset the codec" marker.
///
/// The marker is looked for in the *canonicalised* payload, which is what the
/// reference decoder tests.  Because canonicalisation only keeps bits 0..=138,
/// the nibble examined here is always zero, so in this transport the marker can
/// never actually fire — a quirk of the framing layer that is reproduced here
/// so the two decoders agree bit for bit.
pub fn is_reset(words: &[i16; 9]) -> bool {
    words[8] & 0x0f == 0x08
}

#[inline]
fn test_bit(buf: &[u8], bit: usize) -> bool {
    buf[bit / 8] & (0x80 >> (bit % 8)) != 0
}

#[inline]
fn set_bit(buf: &mut [u8], bit: usize) {
    buf[bit / 8] |= 0x80 >> (bit % 8);
}

/// Parse one hex-encoded frame line, ignoring whitespace.
pub fn parse_hex_line(line: &str) -> Option<[u8; FRAME_BYTES]> {
    let mut out = [0u8; FRAME_BYTES];
    let mut n = 0;
    let mut upper: Option<u8> = None;
    for ch in line.chars() {
        if ch.is_ascii_whitespace() {
            continue;
        }
        let nib = ch.to_digit(16)? as u8;
        match upper {
            None => upper = Some(nib),
            Some(hi) => {
                if n >= FRAME_BYTES {
                    return None;
                }
                out[n] = (hi << 4) | nib;
                n += 1;
                upper = None;
            }
        }
    }
    if upper.is_some() || n != FRAME_BYTES {
        return None;
    }
    Some(out)
}

/// Pack the fourteen parameter fields back into nine words.
///
/// The exact inverse of `unpack`; the suppression flag rides in bit 5 of the
/// last word.
pub fn pack(p: &Params) -> [i16; 9] {
    let f = |i: usize| p.field[i] as u16 as u32;
    let s = p.suppress as u32;
    let w = [
        (f(0) << 8) | (f(1) >> 4),
        ((f(1) & 0xf) << 12) | (f(2) << 4) | (f(3) >> 12),
        ((f(3) & 0xfff) << 4) | (f(4) >> 3),
        ((f(4) & 0x7) << 13) | (f(5) << 8) | (f(6) >> 8),
        ((f(6) & 0xff) << 8) | (f(7) << 1) | (f(8) >> 7),
        ((f(8) & 0x7f) << 9) | (f(9) >> 7),
        ((f(9) & 0x7f) << 9) | (f(10) << 2) | (f(11) >> 3),
        ((f(11) & 0x7) << 13) | (f(12) >> 3),
        ((f(12) & 0x7) << 13) | (f(13) << 6) | (s << 5),
    ];
    let mut out = [0i16; 9];
    for (o, &v) in out.iter_mut().zip(w.iter()) {
        *o = v as u16 as i16;
    }
    out
}

/// Lay nine packed words out as the eighteen transport bytes, most significant
/// byte first — the inverse of [`canonicalize`] for a frame this crate built.
///
/// [`pack`] leaves the four padding bits and the mirrored suppression bit
/// clear, so `canonicalize(&to_bytes(&pack(p)))` returns exactly `pack(p)`.
pub fn to_bytes(words: &[i16; 9]) -> [u8; FRAME_BYTES] {
    let mut out = [0u8; FRAME_BYTES];
    for (i, &w) in words.iter().enumerate() {
        out[2 * i] = (w as u16 >> 8) as u8;
        out[2 * i + 1] = w as u8;
    }
    out
}

/// Render a frame as the 36-character hex line the `.mcelp` container uses.
///
/// This is what [`parse_hex_line`] reads back.
pub fn to_hex_line(frame: &[u8; FRAME_BYTES]) -> String {
    use std::fmt::Write as _;
    let mut line = String::with_capacity(2 * FRAME_BYTES);
    for b in frame {
        let _ = write!(line, "{b:02x}");
    }
    line
}
