//! Long-term band-energy (noise floor) estimates.
//!
//! The encoder keeps two of these, one adapting roughly twice as fast as the
//! other.  Each holds one 32-bit mean energy per critical band and is only ever
//! updated on frames its own voice-activity decision called silent, so over
//! time it converges on the background noise rather than on the speech.  The
//! band levels the analysis works from (`bands::tracked_levels`) are the dB
//! conversion of exactly these numbers.
//!
//! Everything in here runs with saturation on, so every
//! intermediate is clipped to 32 bits rather than wrapped.

use crate::analysis::BINS;
use crate::bands::BANDS;
use crate::fixed::{acc, norm, sat, scale32, shift};
use crate::tables::BAND_INV_WIDTH;

/// The whole reciprocal of each band's width, not the halved table the dB
/// conversion uses.
/// Smallest energy an estimate is allowed to hold, and the value the reference
/// encoder's tables are documented as starting from.
const ENERGY_FLOOR: i64 = (3 << 16) | 8192;
/// Frame counts at which the adaptation weight steps up.
const RAMP: [i16; 2] = [8, 16];

/// One long-term estimate.
#[derive(Clone)]
pub struct NoiseEstimate {
    /// Mean energy per critical band.
    pub energy: [i64; BANDS],
    /// Silent frames seen so far, held at `RAMP[1] + 1`; negative until the
    /// estimate has been seeded.
    frames: i16,
    /// Weight given to the old estimate at each stage of the ramp.  A fresh
    /// estimate follows the input closely and settles down as it gains
    /// confidence.
    weight: [i16; 3],
}

impl NoiseEstimate {
    /// The faster of the encoder's two estimates.
    pub fn fast() -> Self {
        Self::with_weights([6554, 16384, 29491])
    }

    /// The slower one.
    pub fn slow() -> Self {
        Self::with_weights([3277, 13107, 26214])
    }

    fn with_weights(weight: [i16; 3]) -> Self {
        NoiseEstimate {
            energy: [0; BANDS],
            frames: -1,
            weight,
        }
    }

    /// Fold one frame's magnitude spectrum into the estimate.
    ///
    /// `headroom` is the shift the analysis normalised the frame by, and
    /// `active` is this estimate's own voice-activity decision: an active frame
    /// only ever gets to seed a fresh estimate, never to update a settled one.
    pub fn update(&mut self, spectrum: &[i16; BINS], headroom: i16, active: bool) {
        if self.frames < 0 {
            self.frames = 0;
            self.seed(spectrum, headroom);
        }
        if active {
            return;
        }

        let count = self.frames;
        self.frames = count + 1;
        let old = if count <= RAMP[0] {
            self.weight[0]
        } else if count <= RAMP[1] {
            self.weight[1]
        } else {
            self.frames = count;
            self.weight[2]
        };
        let fresh = 32767 - old;

        let mut bin = 0usize;
        for (band, (&width, energy)) in BAND_INV_WIDTH
            .iter()
            .zip(self.energy.iter_mut())
            .enumerate()
        {
            let sum = band_sum(spectrum, band, &mut bin);
            let mean = scale32(sum, width);
            let update = norm(scale32(mean, fresh), -headroom);
            let blended = sat(acc(scale32(*energy, old) + shift(update, 4)));
            *energy = blended.max(ENERGY_FLOOR);
        }
    }

    /// Start the estimate off at the current frame's band energies.
    fn seed(&mut self, spectrum: &[i16; BINS], headroom: i16) {
        let mut bin = 0usize;
        for (band, (&width, energy)) in BAND_INV_WIDTH
            .iter()
            .zip(self.energy.iter_mut())
            .enumerate()
        {
            let sum = band_sum(spectrum, band, &mut bin);
            let mean = scale32(sum, width);
            *energy = norm(mean, 4 - headroom);
        }
    }
}

/// Total of one band's magnitude bins, accumulated in the high half.
fn band_sum(spectrum: &[i16; BINS], band: usize, bin: &mut usize) -> i64 {
    let first = crate::tables::BAND_EDGES[2 * band] as usize;
    let last = crate::tables::BAND_EDGES[2 * band + 1] as usize;
    let mut total = 0i64;
    for _ in 0..=(last - first) {
        total = sat(acc(total + ((spectrum[*bin] as i64) << 16)));
        *bin += 1;
    }
    total
}
