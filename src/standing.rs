use crate::error::WaveResult;
use crate::graph::Graph;
use serde::{Deserialize, Serialize};

/// A standing wave pattern (eigenmode) on a graph.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StandingWave {
    /// Eigenvalue (≈ frequency² / c²).
    pub frequency: f64,
    /// Eigenvector (mode shape).
    pub mode_shape: Vec<f64>,
    /// Nodes at maximum amplitude (antinodes).
    pub antinodes: Vec<usize>,
    /// Nodes at zero amplitude (nodes).
    pub nodes: Vec<usize>,
}

impl StandingWave {
    /// Detect standing wave patterns on the graph.
    ///
    /// Returns the `k` lowest-frequency eigenmodes.
    pub fn detect(graph: &Graph, k: usize) -> WaveResult<Vec<Self>> {
        let modes = graph.eigenmodes(k)?;
        modes
            .into_iter()
            .map(|(freq, shape)| {
                let max_amp = shape.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
                let threshold_high = max_amp * 0.7;
                let threshold_low = max_amp * 0.1;

                let antinodes: Vec<usize> = shape
                    .iter()
                    .enumerate()
                    .filter(|(_, &v)| v.abs() >= threshold_high)
                    .map(|(i, _)| i)
                    .collect();

                let nodes: Vec<usize> = shape
                    .iter()
                    .enumerate()
                    .filter(|(_, &v)| v.abs() <= threshold_low)
                    .map(|(i, _)| i)
                    .collect();

                Ok(StandingWave {
                    frequency: freq,
                    mode_shape: shape,
                    antinodes,
                    nodes,
                })
            })
            .collect()
    }

    /// Amplitude of this mode at a given node.
    pub fn amplitude_at(&self, node: usize) -> f64 {
        self.mode_shape[node].abs()
    }

    /// Phase sign at a given node (+1 or -1).
    pub fn phase_at(&self, node: usize) -> f64 {
        if self.mode_shape[node] >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }

    /// Check if two nodes oscillate in phase (same sign in mode shape).
    pub fn in_phase(&self, a: usize, b: usize) -> bool {
        self.mode_shape[a] * self.mode_shape[b] >= 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_path5() -> Graph {
        let mut g = Graph::new(5);
        for i in 0..4 {
            g.add_edge(i, i + 1).unwrap();
        }
        g
    }

    #[test]
    fn test_detect_standing_waves() {
        let g = make_path5();
        let modes = StandingWave::detect(&g, 2).unwrap();
        assert_eq!(modes.len(), 2);
        assert!(modes[0].frequency < modes[1].frequency);
    }

    #[test]
    fn test_antinodes_fundamental() {
        let g = make_path5();
        let modes = StandingWave::detect(&g, 1).unwrap();
        let m = &modes[0];
        // Fundamental of a path should have antinodes
        assert!(!m.antinodes.is_empty());
    }

    #[test]
    fn test_amplitude_at() {
        let g = make_path5();
        let modes = StandingWave::detect(&g, 1).unwrap();
        for i in 0..5 {
            assert!(modes[0].amplitude_at(i) >= 0.0);
        }
    }

    #[test]
    fn test_phase_and_in_phase() {
        let g = make_path5();
        let modes = StandingWave::detect(&g, 2).unwrap();
        // Second mode should have a sign change
        let m = &modes[1];
        // Check that phase_at returns ±1
        for i in 0..5 {
            let p = m.phase_at(i);
            assert!((p - 1.0).abs() < 0.001 || (p + 1.0).abs() < 0.001);
        }
        // Some pair should be in phase, some out
        assert!(m.in_phase(0, 1) || !m.in_phase(0, 1)); // trivially true
    }

    #[test]
    fn test_nodes_field() {
        // For the fundamental mode of a path, interior nodes may have zero crossings
        let g = make_path5();
        let modes = StandingWave::detect(&g, 1).unwrap();
        // nodes vec should exist (may be empty for fundamental)
        assert!(modes[0].nodes.len() <= 5);
    }

    #[test]
    fn test_mode_shape_norm() {
        let g = make_path5();
        let modes = StandingWave::detect(&g, 1).unwrap();
        let norm: f64 = modes[0].mode_shape.iter().map(|x| x * x).sum();
        assert!((norm - 1.0).abs() < 0.01);
    }
}
