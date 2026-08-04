//! Interpolation filters and lag tables of the two pitch searches.

/// Fractional-delay filter for the closed-loop lag search, read with a
/// stride of three from either end: [`LAG_INTERP_FORWARD`] for the taps
/// ahead of the sample and [`LAG_INTERP_BACKWARD`] for those behind it.
pub static LAG_INTERP: [i16; 13] = [
    29519, 24906, 13896, 2755, -3459, -3969, -1561, 534, 1023, 516, 0, -194, 0,
];

/// 1/3-resolution pitch interpolation filter, 31 taps.
pub static PITCH_INTERP: [i16; 31] = [
    29443, 25207, 14701, 3143, -4402, -5850, -2783, 1211, 3130, 2259, 0, -1652, -1666, -464, 756,
    1099, 550, -245, -634, -451, 0, 308, 296, 78, -120, -165, -79, 34, 91, 70, 0,
];

/// Short interpolation filter bank of the long-term postfilter's search.
pub static SEARCH_FILTER: [i16; 28] = [
    -189, 2873, 31650, -1598, -485, 7041, 28469, -2148, -934, 12266, 23705, -1993, -1493, 18050,
    18050, -1493, -1993, 23705, 12266, -934, -2148, 28469, 7041, -485, -1598, 31650, 2873, -189,
];

/// Long interpolation filter bank, seven phases of sixteen taps.
pub static REFINE_FILTER: [i16; 112] = [
    -41, 72, -157, 315, -580, 1023, -1875, 4439, 31915, -3391, 1595, -888, 501, -267, 130, -60,
    -78, 147, -318, 631, -1151, 2030, -3774, 9639, 29436, -5580, 2727, -1528, 859, -454, 218, -102,
    -107, 212, -456, 892, -1615, 2850, -5393, 15206, 25569, -6550, 3303, -1861, 1041, -544, 258,
    -123, -124, 253, -539, 1044, -1877, 3319, -6414, 20676, 20676, -6414, 3319, -1877, 1044, -539,
    253, -124, -123, 258, -544, 1041, -1861, 3303, -6550, 25569, 15206, -5393, 2850, -1615, 892,
    -456, 212, -107, -102, 218, -454, 859, -1528, 2727, -5580, 29436, 9639, -3774, 2030, -1151,
    631, -318, 147, -78, -60, 130, -267, 501, -888, 1595, -3391, 31915, 4439, -1875, 1023, -580,
    315, -157, 72, -41,
];

/// Lag to correlation-slot maps: [`LAG_SLOTS_LOW`] gives the first slot a
/// lag covers and [`LAG_SLOTS_HIGH`] the last.  The two overlap, so they
/// are one table with two offsets.
pub static LAG_SLOTS: [i16; 153] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
];

/// Offset of the forward half of [`LAG_INTERP`].
pub const LAG_INTERP_FORWARD: usize = 0;

/// Offset of the backward half of [`LAG_INTERP`].
pub const LAG_INTERP_BACKWARD: usize = 3;

/// Offset of the first-slot half of [`LAG_SLOTS`].
pub const LAG_SLOTS_LOW: usize = 0;

/// Offset of the last-slot half of [`LAG_SLOTS`].
pub const LAG_SLOTS_HIGH: usize = 8;
