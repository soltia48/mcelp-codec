//! Critical-band geometry and the noise-suppression curves.

/// Per-band smoothing coefficient pairs, twenty for speech frames followed
/// by twenty for silence.
pub static BAND_SMOOTH: [i16; 40] = [
    16368, 16355, 16336, 16317, 16299, 16280, 16255, 16230, 16198, 16167, 16136, 16092, 16048,
    15991, 15935, 15866, 15784, 15684, 15596, 15565, 15565, 15552, 15533, 15514, 15494, 15475,
    15450, 15424, 15392, 15360, 15328, 15283, 15238, 15181, 15123, 15053, 14970, 14867, 14778,
    14746,
];

/// Reciprocal of each band's width, Q15.
pub static BAND_INV_WIDTH: [i16; 20] = [
    32767, 10923, 10923, 10923, 10923, 8192, 8192, 8192, 6554, 6554, 5461, 4681, 4096, 3641, 3277,
    2731, 2341, 1820, 3641, 32767,
];

/// The same reciprocal halved, for the paths where the fractional multiply
/// doubles it back.
pub static BAND_INV_WIDTH_HALVED: [i16; 20] = [
    16384, 5461, 5461, 5461, 5461, 4096, 4096, 4096, 3277, 3277, 2731, 2341, 2048, 1820, 1638,
    1365, 1170, 910, 1820, 16384,
];

/// Per-band weights of the two exponentials the suppression mixes.
pub static MIX_COEFFS: [i16; 20] = [
    0, 256, 640, 1024, 1408, 1792, 2304, 2816, 3456, 4096, 4736, 5632, 6528, 7680, 8832, 10240,
    11904, 13952, 15744, 16384,
];

/// Band edges, stored as `(first, last)` bin pairs.
pub static BAND_EDGES: [i16; 40] = [
    0, 0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 16, 17, 20, 21, 24, 25, 29, 30, 34, 35, 40, 41, 47, 48, 55,
    56, 64, 65, 74, 75, 86, 87, 100, 101, 118, 119, 127, 128, 128,
];

/// Three sets of per-band suppression depths, in dB, deepest first.
pub static DEPTH_TABLES: [i16; 60] = [
    -2048, -2080, -2128, -2176, -2224, -2272, -2336, -2400, -2480, -2560, -2640, -2752, -2864,
    -3008, -3152, -3328, -3536, -3792, -4016, -4096, -1024, -1056, -1104, -1152, -1200, -1248,
    -1312, -1376, -1456, -1536, -1616, -1728, -1840, -1984, -2128, -2304, -2512, -2768, -2992,
    -3072, 0, -16, -40, -64, -88, -112, -144, -176, -216, -256, -296, -352, -408, -480, -552, -640,
    -744, -872, -984, -1024,
];

/// Per-band smoothing of the noise scale, one set for speech and one for
/// silence.  The band-0 pass reads one word before its set, so both start
/// one word in: see [`SCALE_SMOOTH_ACTIVE`] and [`SCALE_SMOOTH_IDLE`].
pub static SCALE_SMOOTH: [i16; 40] = [
    32735, 32726, 32711, 32696, 32682, 32667, 32648, 32628, 32604, 32580, 32555, 32521, 32487,
    32443, 32400, 32346, 32283, 32205, 32137, 32113, 32440, 32394, 32325, 32256, 32187, 32118,
    32026, 31933, 31818, 31703, 31588, 31427, 31265, 31058, 30851, 30597, 30298, 29929, 29606,
    29491,
];

/// Offset of the speech set within [`SCALE_SMOOTH`].
pub const SCALE_SMOOTH_ACTIVE: usize = 1;
/// Offset of the silence set within [`SCALE_SMOOTH`].
pub const SCALE_SMOOTH_IDLE: usize = 21;

/// Weight each band carries into the summary level.
pub static BAND_WEIGHT: [i16; 20] = [
    8192, 8115, 8000, 7884, 7769, 7654, 7500, 7347, 7155, 6963, 6771, 6502, 6233, 5888, 5542, 5120,
    4620, 4006, 3468, 3276,
];

/// Fractional part table of the base-2 logarithm.  Its own entries share
/// storage with the tables that follow, and a high index reads one past
/// them, so the whole run up to [`POW2_CURVE`] is kept here.
pub static LOG2_FRACTION: [i16; 65] = [
    32440, 31130, 26214, 16384, 3277, 0, 0, 0, 0, 328, 399, 507, 614, 722, 829, 973, 1116, 1295,
    1475, 1654, 1905, 2156, 2478, 2801, 3195, 3661, 4234, 4736, 4915, 5615, 21871, 32701, 0, 1455,
    2866, 4236, 5568, 6863, 8124, 9352, 10549, 11716, 12855, 13967, 15054, 16117, 17156, 18172,
    19167, 20142, 21097, 22033, 22951, 23852, 24735, 25603, 26455, 27291, 28113, 28922, 29716,
    30497, 31266, 32023, 32767,
];

/// How far the scale is allowed to move, by confidence.  A negative score
/// reads below the entry for zero, so the window starts early and
/// [`CONFIDENCE_ZERO`] is where a score of nought sits.
pub static CONFIDENCE: [i16; 13] = [
    32440, 31130, 26214, 16384, 3277, 0, 0, 0, 0, 328, 399, 507, 614,
];

/// Index of a zero score in [`CONFIDENCE`].
pub const CONFIDENCE_ZERO: usize = 4;

/// Per-band adaptation rate of the scale smoother.
pub static STATE_RATE: [i16; 20] = [
    328, 399, 507, 614, 722, 829, 973, 1116, 1295, 1475, 1654, 1905, 2156, 2478, 2801, 3195, 3661,
    4234, 4736, 4915,
];

/// Quadratic coefficients of the attenuation curve.
pub static CURVE: [i16; 36] = [
    5615, 21871, 32701, 0, 1455, 2866, 4236, 5568, 6863, 8124, 9352, 10549, 11716, 12855, 13967,
    15054, 16117, 17156, 18172, 19167, 20142, 21097, 22033, 22951, 23852, 24735, 25603, 26455,
    27291, 28113, 28922, 29716, 30497, 31266, 32023, 32767,
];

/// Interpolation table for two-to-the-power, 33 entries.
pub static POW2_CURVE: [i16; 33] = [
    16384, 16743, 17109, 17484, 17867, 18258, 18658, 19066, 19484, 19911, 20347, 20792, 21247,
    21713, 22188, 22674, 23170, 23678, 24196, 24726, 25268, 25821, 26386, 26964, 27554, 28158,
    28774, 29405, 30048, 30706, 31379, 32066, 32767,
];
