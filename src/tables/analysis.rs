//! Window and rotator tables of the spectral analysis.

/// Taper applied to both ends of the analysis window.
pub static TAPER: [i16; 10] = [0, 888, 3455, 7423, 12362, 17737, 22965, 27481, 30793, 32545];

/// Sine table the transform reads: `sin` at an index and `cos` a quarter
/// period further on.  The rotator is re-seeded from the same table every
/// 32 bins, which is what [`SINE_ANCHOR_STEP`] steps through.
pub static SINE: [i16; 98] = [
    0, 1608, 3212, 4808, 6393, 7962, 9512, 11039, 12540, 14010, 15447, 16846, 18205, 19520, 20788,
    22006, 23170, 24279, 25330, 26320, 27246, 28106, 28899, 29622, 30274, 30853, 31357, 31786,
    32138, 32413, 32610, 32729, 32767, 32729, 32610, 32413, 32138, 31786, 31357, 30853, 30274,
    29622, 28899, 28106, 27246, 26320, 25330, 24279, 23170, 22006, 20788, 19520, 18205, 16846,
    15447, 14010, 12540, 11039, 9512, 7962, 6393, 4808, 3212, 1608, 0, -1608, -3212, -4808, -6393,
    -7962, -9512, -11039, -12540, -14010, -15447, -16846, -18205, -19520, -20788, -22006, -23170,
    -24279, -25330, -26320, -27246, -28106, -28899, -29622, -30274, -30853, -31357, -31786, -32138,
    -32413, -32610, -32729, -32768, 32735,
];

/// How far along [`SINE`] each re-seed of the rotator moves.
pub const SINE_ANCHOR_STEP: usize = 16;
