//! Wave conservation: spectral wave propagation on graphs

pub mod spectral;
pub mod wave;

pub use spectral::eigenvalues;
pub use wave::WaveState;
