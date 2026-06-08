use std::fmt;

/// Errors that can arise during wave-graph computations.
#[derive(Debug)]
pub enum WaveError {
    /// Graph has zero nodes.
    EmptyGraph,
    /// Graph is disconnected (no wave propagation possible).
    Disconnected,
    /// Invalid parameter (e.g. negative damping).
    InvalidParameter(String),
    /// Eigenvalue computation did not converge.
    NoConvergence { iterations: usize },
    /// Index out of bounds.
    IndexOutOfBounds { index: usize, len: usize },
}

impl fmt::Display for WaveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WaveError::EmptyGraph => write!(f, "graph has zero nodes"),
            WaveError::Disconnected => write!(f, "graph is disconnected"),
            WaveError::InvalidParameter(msg) => write!(f, "invalid parameter: {msg}"),
            WaveError::NoConvergence { iterations } => {
                write!(f, "eigenvalue computation did not converge after {iterations} iterations")
            }
            WaveError::IndexOutOfBounds { index, len } => {
                write!(f, "index {index} out of bounds (len {len})")
            }
        }
    }
}

impl std::error::Error for WaveError {}

/// Convenience alias for results in this crate.
pub type WaveResult<T> = Result<T, WaveError>;
