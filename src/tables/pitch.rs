//! 10-tap pitch interpolation filter coefficients.

pub const INTERP_FILTER_FRACTIONS: usize = 3;
pub const INTERP_FILTER_TAPS: usize = 10;
pub static INTERP_FILTER: [[i16; 10]; 3] = [
    [ 29443,   3143,  -2783,   2259,  -1666,   1099,   -634,    308,   -120,     34],
    [ 25207,  -4402,   1211,      0,   -464,    550,   -451,    296,   -165,     91],
    [ 14701,  -5850,   3130,  -1652,    756,   -245,      0,     78,    -79,     70],
];
