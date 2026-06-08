use crate::error::WaveResult;
use crate::graph::Graph;
use crate::wave::{WaveEquation, WaveState};
use serde::{Deserialize, Serialize};

/// Report for a node identified as a bottleneck.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BottleneckReport {
    pub node: usize,
    /// How slow waves reach this node (time to exceed amplitude threshold).
    pub wave_delay: f64,
    /// Normalized connectivity = degree / n.
    pub connectivity: f64,
}

/// Detect communication bottlenecks by measuring how long it takes for a wave
/// pulse to reach each node from a source.
///
/// Nodes with high wave delay are bottlenecks — they're poorly connected and
/// information arrives slowly.
pub fn detect_bottlenecks(
    graph: &Graph,
    source: usize,
    amplitude_threshold: f64,
    max_steps: usize,
) -> WaveResult<Vec<BottleneckReport>> {
    let n = graph.n();
    if n == 0 {
        return Err(crate::error::WaveError::EmptyGraph);
    }

    let we = WaveEquation::new(graph.clone(), 0.0, 0.001)?;
    let initial = WaveState {
        displacement: {
            let mut d = vec![0.0; n];
            d[source] = 1.0;
            d
        },
        velocity: vec![0.0; n],
        time: 0.0,
    };

    let mut arrival_time = vec![f64::INFINITY; n];
    arrival_time[source] = 0.0;

    let mut state = initial;
    for _step in 1..=max_steps {
        state = we.step(&state);
        for (i, arrival) in arrival_time.iter_mut().enumerate() {
            if arrival.is_infinite() && state.displacement[i].abs() > amplitude_threshold {
                *arrival = state.time;
            }
        }
    }

    let max_delay = arrival_time.iter().cloned().fold(0.0_f64, f64::max);

    let mut reports: Vec<BottleneckReport> = (0..n)
        .map(|i| {
            let delay = if arrival_time[i].is_infinite() {
                max_delay * 2.0 // never reached = maximum bottleneck
            } else {
                arrival_time[i]
            };
            BottleneckReport {
                node: i,
                wave_delay: delay,
                connectivity: graph.degree(i) as f64 / n as f64,
            }
        })
        .collect();

    // Sort by wave delay descending (worst bottlenecks first)
    reports.sort_by(|a, b| b.wave_delay.partial_cmp(&a.wave_delay).unwrap());
    Ok(reports)
}

/// Quick check: return nodes whose wave delay is more than `threshold` standard
/// deviations above the mean.
#[allow(dead_code)]
pub fn bottleneck_outliers(
    reports: &[BottleneckReport],
    stdev_threshold: f64,
) -> Vec<&BottleneckReport> {
    let n = reports.len() as f64;
    if n == 0.0 {
        return vec![];
    }
    let mean: f64 = reports.iter().map(|r| r.wave_delay).sum::<f64>() / n;
    let variance: f64 =
        reports.iter().map(|r| (r.wave_delay - mean).powi(2)).sum::<f64>() / n;
    let stdev = variance.sqrt();

    if stdev < 1e-12 {
        return vec![];
    }

    reports
        .iter()
        .filter(|r| (r.wave_delay - mean) / stdev > stdev_threshold)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_path6() -> Graph {
        let mut g = Graph::new(6);
        for i in 0..5 {
            g.add_edge(i, i + 1).unwrap();
        }
        g
    }

    #[test]
    fn test_detect_bottlenecks_path() {
        let g = make_path6();
        let reports = detect_bottlenecks(&g, 0, 0.01, 5000).unwrap();
        assert_eq!(reports.len(), 6);
        // Source should have delay 0
        let source_report = reports.iter().find(|r| r.node == 0).unwrap();
        assert_eq!(source_report.wave_delay, 0.0);
        // Farthest node should have highest delay
        assert_eq!(reports[0].node, 5);
    }

    #[test]
    fn test_connectivity_values() {
        let g = make_path6();
        let reports = detect_bottlenecks(&g, 0, 0.01, 5000).unwrap();
        // Interior nodes have connectivity 2/6
        let interior = reports.iter().find(|r| r.node == 2).unwrap();
        assert!((interior.connectivity - 2.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_bottleneck_star_graph() {
        // Star: center 0, leaves 1-4
        let mut g = Graph::new(5);
        for i in 1..5 {
            g.add_edge(0, i).unwrap();
        }
        let reports = detect_bottlenecks(&g, 0, 0.01, 5000).unwrap();
        // Source is center, all leaves should have similar delay
        let leaf_delays: Vec<f64> = reports.iter().filter(|r| r.node > 0).map(|r| r.wave_delay).collect();
        let min_d = leaf_delays.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_d = leaf_delays.iter().cloned().fold(0.0_f64, f64::max);
        assert!((max_d - min_d) / max_d < 0.3);
    }

    #[test]
    fn test_bottleneck_outliers() {
        // Barbell graph: two clusters connected by a single bridge node
        let mut g = Graph::new(7);
        // Left cluster: 0-1-2
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        g.add_edge(0, 2).unwrap();
        // Bridge: 3
        g.add_edge(2, 3).unwrap();
        g.add_edge(3, 4).unwrap();
        // Right cluster: 4-5-6
        g.add_edge(4, 5).unwrap();
        g.add_edge(5, 6).unwrap();
        g.add_edge(4, 6).unwrap();
        let reports = detect_bottlenecks(&g, 0, 0.01, 5000).unwrap();
        let outliers = bottleneck_outliers(&reports, 0.8);
        // The far right nodes (5, 6) should be outliers
        assert!(!outliers.is_empty());
    }

    #[test]
    fn test_empty_graph_bottleneck() {
        let g = Graph::new(0);
        assert!(detect_bottlenecks(&g, 0, 0.01, 100).is_err());
    }
}
