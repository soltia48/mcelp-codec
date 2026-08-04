//! Linear to mu-law conversion.
//!
//! Plain G.711 mu-law with a bias of 33 and the usual complemented output, but
//! computed the way the reference does it: the biased magnitude is normalised by
//! shifting until the segment number falls out, and the segment and mantissa
//! fields are complemented individually rather than the finished byte.

/// Bias added before segmenting, in the 13-bit domain.
const BIAS: i16 = 33;
/// Largest representable biased magnitude.
const CLIP: i16 = 8191;

/// Convert one linear sample to a mu-law byte.
pub fn from_linear(sample: i16) -> u8 {
    let negative = sample < 0;
    let mut magnitude = if negative { !sample } else { sample };
    magnitude = ((magnitude >> 2) + BIAS).min(CLIP);

    let mut high = magnitude >> 6;
    magnitude >>= 1;
    let mut segment = 1;
    while high != 0 {
        high >>= 1;
        segment += 1;
        magnitude >>= 1;
    }

    let field = ((8 - segment) as i16) << 4;
    let mantissa = 15 - (magnitude & 15);
    let byte = (field | mantissa) as u8;
    if negative { byte } else { byte | 0x80 }
}

/// Convert a whole frame in place.
pub fn frame_from_linear(samples: &[i16], out: &mut [u8]) {
    for (o, &s) in out.iter_mut().zip(samples) {
        *o = from_linear(s);
    }
}
