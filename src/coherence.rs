use crate::wave::WaveEquation;
use crate::wave::WaveState;

/// Compute the coherence ratio: fraction of initial amplitude conserved after `steps` steps.
///
/// Returns a value in `[0, 1]`. With no damping, the ratio should be close to 1.0
/// (energy conserved). With damping, it drops over time.
pub fn coherence_ratio(eq: &WaveEquation, initial: &WaveState, steps: usize) -> f64 {
    let e0 = eq.energy(initial);
    if e0 <= 0.0 {
        return 1.0;
    }
    let final_state = eq.simulate(initial, steps);
    let ef = eq.energy(&final_state);
    (ef / e0).clamp(0.0, 1.0)
}

/// Measure amplitude decay at a specific node over `steps` steps.
///
/// Returns the ratio of final amplitude at `node` to the initial max amplitude.
#[allow(dead_code)]
pub fn node_amplitude_ratio(eq: &WaveEquation, initial: &WaveState, node: usize, steps: usize) -> f64 {
    let max_amp = initial.displacement.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    if max_amp <= 0.0 {
        return 0.0;
    }
    let final_state = eq.simulate(initial, steps);
    final_state.displacement[node].abs() / max_amp
}

/// Compute coherence over multiple time windows, returning the ratio at each checkpoint.
#[allow(dead_code)]
pub fn coherence_timeline(eq: &WaveEquation, initial: &WaveState, checkpoints: &[usize]) -> Vec<(usize, f64)> {
    let e0 = eq.energy(initial);
    let mut results = Vec::new();
    let mut state = initial.clone();
    let mut prev_step = 0;

    for &target in checkpoints {
        let delta = target - prev_step;
        for _ in 0..delta {
            state = eq.step(&state);
        }
        let ef = eq.energy(&state);
        let ratio = if e0 > 0.0 { (ef / e0).clamp(0.0, 1.0) } else { 1.0 };
        results.push((target, ratio));
        prev_step = target;
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    fn make_path4() -> (Graph, WaveEquation) {
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        g.add_edge(2, 3).unwrap();
        let we = WaveEquation::new(g, 0.0, 0.001).unwrap();
        (we.graph.clone(), we)
    }

    #[test]
    fn test_coherence_no_damping() {
        let (_, we) = make_path4();
        let state = we.pulse(0, 1.0).unwrap();
        let ratio = coherence_ratio(&we, &state, 100);
        assert!(ratio > 0.9, "coherence with no damping should be >0.9, got {ratio}");
    }

    #[test]
    fn test_coherence_with_damping() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        g.add_edge(2, 3).unwrap();
        let we = WaveEquation::new(g, 1.0, 0.01).unwrap();
        let state = we.pulse(0, 1.0).unwrap();
        let ratio = coherence_ratio(&we, &state, 200);
        assert!(ratio < 0.5, "coherence with high damping should be <0.5, got {ratio}");
    }

    #[test]
    fn test_coherence_zero_energy() {
        let (_, we) = make_path4();
        let state = WaveState {
            displacement: vec![0.0; 4],
            velocity: vec![0.0; 4],
            time: 0.0,
        };
        let ratio = coherence_ratio(&we, &state, 100);
        assert_eq!(ratio, 1.0);
    }

    #[test]
    fn test_node_amplitude_ratio() {
        let (_, we) = make_path4();
        let state = we.pulse(0, 1.0).unwrap();
        let ratio = node_amplitude_ratio(&we, &state, 0, 200);
        assert!(ratio >= 0.0 && ratio <= 1.0);
    }

    #[test]
    fn test_coherence_timeline() {
        let (_, we) = make_path4();
        let state = we.pulse(0, 1.0).unwrap();
        let timeline = coherence_timeline(&we, &state, &[50, 100, 200]);
        assert_eq!(timeline.len(), 3);
        // Energy should be monotonically non-increasing (within numerical precision)
        for w in timeline.windows(2) {
            assert!(w[0].1 >= w[1].1 - 0.05);
        }
    }
}
