//! # wave-conservation
//!
//! **Spectral wave propagation on graphs — wave speed = √λ₂, standing waves
//! reveal the eigenvalue spectrum, and conservation ratio predicts coherence.**
//!
//! This crate implements the discrete wave equation `d²u/dt² = -L·u - γ·du/dt`
//! on graph structures, where `L` is the graph Laplacian. The wave behavior is
//! governed entirely by the Laplacian eigenvalue spectrum:
//!
//! - **Wave speed** = √λ₂ (Fiedler eigenvalue)
//! - **Resonance frequencies** = √λᵢ for each eigenvalue
//! - **Conservation ratio** = λ₂/λₙ predicts coherence halflife
//!
//! # Key Insight
//!
//! The wave equation on a graph IS a spectral probe. You don't need to compute
//! eigenvalues separately — send a wave through the graph and read the spectrum
//! from the response. Each eigenvalue appears as a resonance peak in the
//! frequency sweep, and the Fiedler eigenvalue sets the wave speed.
//!
//! # Modules
//!
//! - [`wave`] — Wave state, velocity Verlet integration, energy and coherence
//! - [`spectral`] — Eigenvalue computation, frequency sweep, Fiedler vector
//! - [`experiments`] — Graph generators and verification experiments
//!
//! # Quick Start
//!
//! ```rust
//! use wave_conservation::wave::WaveState;
//! use wave_conservation::experiments;
//!
//! let adj = experiments::path_graph(20);
//! let mut wave = WaveState::new(adj).with_damping(0.001);
//! wave.pulse(0, 1.0);
//!
//! for _ in 0..1000 {
//!     wave.step(0.01);
//! }
//! ```

pub mod spectral;
pub mod wave;
pub mod experiments;

pub use spectral::eigenvalues;
pub use wave::WaveState;
