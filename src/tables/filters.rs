//! Fixed filter coefficients: input conditioning, weighting, output.

/// Feed-forward taps of the shaping high-pass.
pub static HIGHPASS_FF: [i16; 4] = [3798, -7596, 3798, 4096];

/// Feedback taps of the same.
pub static HIGHPASS_FB: [i16; 2] = [7807, -3733];

/// Coefficients of the output filter's feed-forward part.
pub static OUT_FIR: [i16; 5] = [7472, -22381, 22381, -7472, 8192];

/// Coefficients of its feedback part.
pub static OUT_IIR: [i16; 3] = [-23033, 21666, -6815];

/// Symmetric six-tap smoothing filter.
pub static SMOOTH_FIR: [i16; 7] = [2528, 12466, 24759, 24759, 12466, 2528, 4096];

/// Feedback part of the tilt filter.
pub static SMOOTH_IIR: [i16; 5] = [16303, 26485, 21883, 9177, 1560];

/// Coefficients of the input high-pass: two recursive then three direct taps.
pub static INPUT_HIGHPASS: [i16; 5] = [7807, -3733, 3798, -7596, 3798];
