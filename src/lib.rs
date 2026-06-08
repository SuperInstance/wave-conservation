//! # wave-conservation
//!
//! Spectral wave propagation on graphs. Wave speed on a graph is `√λ₂` where `λ₂`
//! is the algebraic connectivity (Fiedler value). Standing wave patterns reveal the
//! eigenvalue spectrum, and nodes where waves arrive late indicate communication
//! bottlenecks in agent networks.
//!
//! ## Core ideas
//!
//! - **Wave equation on graphs:** `u'' = -c² L u` where `L` is the graph Laplacian and `c = √λ₂`.
//! - **Coherence:** measures amplitude conservation over time.
//! - **Standing waves:** eigenmode analysis via power iteration + deflation.
//! - **Bottleneck detection:** nodes with high wave delay are communication bottlenecks.
//! - **MIDI export:** map wave patterns to MIDI (position→pitch, amplitude→velocity).

mod error;
mod graph;
mod wave;
mod coherence;
mod standing;
mod bottleneck;
mod midi;

pub use error::{WaveError, WaveResult};
pub use graph::Graph;
pub use wave::{WaveState, WaveEquation};
pub use coherence::coherence_ratio;
pub use standing::StandingWave;
pub use bottleneck::{BottleneckReport, detect_bottlenecks};
pub use midi::{MidiEvent, wave_to_midi};
