//! The codec's constant tables: codebooks, filter coefficients and the lookup
//! curves the arithmetic reads.
//!
//! One array per table, grouped by the stage that reads it.  Constants that
//! describe the *code* rather than the codec — frame geometry, thresholds,
//! Q-format shifts — live beside the stage that uses them, not here.

pub mod analysis;
pub mod bands;
pub mod fcb;
pub mod filters;
pub mod gain;
pub mod lpc;
pub mod lsf;
pub mod pitch;
pub mod ulaw;

pub use analysis::*;
pub use bands::*;
pub use fcb::*;
pub use filters::*;
pub use gain::*;
pub use lpc::*;
pub use lsf::*;
pub use pitch::*;
pub use ulaw::*;
