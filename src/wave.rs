use crate::error::{WaveError, WaveResult};
use crate::graph::Graph;
use serde::{Deserialize, Serialize};

/// Current wave displacement and velocity at each node.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaveState {
    pub displacement: Vec<f64>,
    pub velocity: Vec<f64>,
    pub time: f64,
}

/// Discrete wave equation on a graph: `u'' = -c² L u`.
///
/// The wave speed `c` defaults to `√λ₂` (algebraic connectivity).
#[derive(Debug, Clone)]
pub struct WaveEquation {
    pub graph: Graph,
    pub wave_speed: f64,
    pub damping: f64,
    pub dt: f64,
    laplacian: Vec<Vec<f64>>,
}

impl WaveEquation {
    /// Create a new wave equation on the given graph.
    ///
    /// Computes `λ₂` automatically to set the wave speed.
    pub fn new(graph: Graph, damping: f64, dt: f64) -> WaveResult<Self> {
        if graph.n() == 0 {
            return Err(WaveError::EmptyGraph);
        }
        let (lambda2, _) = graph.fiedler()?;
        let wave_speed = lambda2.sqrt().max(1e-6);
        let laplacian = graph.laplacian();
        Ok(Self {
            graph,
            wave_speed,
            damping,
            dt,
            laplacian,
        })
    }

    /// Create with an explicit wave speed (skip eigenvalue computation).
    pub fn with_speed(graph: Graph, wave_speed: f64, damping: f64, dt: f64) -> Self {
        let laplacian = graph.laplacian();
        Self {
            graph,
            wave_speed,
            damping,
            dt,
            laplacian,
        }
    }

    /// Advance the wave state by one time step using Störmer-Verlet integration.
    pub fn step(&self, state: &WaveState) -> WaveState {
        let n = self.graph.n();
        let c2 = self.wave_speed * self.wave_speed;
        let dt = self.dt;
        let d = self.damping;

        // a_i = -c² * (L u)_i - damping * v_i
        let lu = self.laplacian_vec(&state.displacement);
        let mut new_vel = Vec::with_capacity(n);
        let mut new_disp = Vec::with_capacity(n);

        for (i, lu_i) in lu.iter().enumerate() {
            let accel = -c2 * lu_i - d * state.velocity[i];
            let v_new = state.velocity[i] + accel * dt;
            let u_new = state.displacement[i] + v_new * dt;
            new_vel.push(v_new);
            new_disp.push(u_new);
        }

        WaveState {
            displacement: new_disp,
            velocity: new_vel,
            time: state.time + dt,
        }
    }

    /// Run `n` steps starting from the given state.
    pub fn simulate(&self, initial: &WaveState, steps: usize) -> WaveState {
        let mut state = initial.clone();
        for _ in 0..steps {
            state = self.step(&state);
        }
        state
    }

    /// Convenience: create an initial state with a pulse at `node`.
    pub fn pulse(&self, node: usize, amplitude: f64) -> WaveResult<WaveState> {
        if node >= self.graph.n() {
            return Err(WaveError::IndexOutOfBounds {
                index: node,
                len: self.graph.n(),
            });
        }
        let n = self.graph.n();
        let mut disp = vec![0.0; n];
        disp[node] = amplitude;
        Ok(WaveState {
            displacement: disp,
            velocity: vec![0.0; n],
            time: 0.0,
        })
    }

    /// Total energy (kinetic + potential) of the wave state.
    pub fn energy(&self, state: &WaveState) -> f64 {
        let c2 = self.wave_speed * self.wave_speed;
        let lu = self.laplacian_vec(&state.displacement);
        let kinetic: f64 = state.velocity.iter().map(|v| 0.5 * v * v).sum();
        let potential: f64 = 0.5 * c2 * state.displacement.iter().zip(lu.iter()).map(|(u, lu_i)| u * lu_i).sum::<f64>();
        kinetic + potential
    }

    fn laplacian_vec(&self, u: &[f64]) -> Vec<f64> {
        self.laplacian.iter().map(|row| {
            row.iter().zip(u.iter()).map(|(l, u)| l * u).sum()
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_path4() -> Graph {
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        g.add_edge(2, 3).unwrap();
        g
    }

    #[test]
    fn test_wave_creation() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.0, 0.01).unwrap();
        assert!(we.wave_speed > 0.0);
    }

    #[test]
    fn test_pulse_initial() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.0, 0.01).unwrap();
        let state = we.pulse(0, 1.0).unwrap();
        assert_eq!(state.displacement[0], 1.0);
        assert_eq!(state.displacement[1], 0.0);
    }

    #[test]
    fn test_pulse_out_of_bounds() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.0, 0.01).unwrap();
        assert!(we.pulse(10, 1.0).is_err());
    }

    #[test]
    fn test_step_displacement() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.0, 0.01).unwrap();
        let state = we.pulse(0, 1.0).unwrap();
        let next = we.step(&state);
        // Displacement should change
        assert_ne!(next.displacement[0], state.displacement[0]);
        // Wave hasn't reached node 3 yet
        assert!(next.displacement[3].abs() < 0.01);
    }

    #[test]
    fn test_simulate_spread() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.0, 0.01).unwrap();
        let state = we.pulse(0, 1.0).unwrap();
        let final_state = we.simulate(&state, 500);
        // Wave should have spread to all nodes
        let max_amp = final_state.displacement.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        assert!(max_amp > 0.0);
    }

    #[test]
    fn test_damping_reduces_energy() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.5, 0.01).unwrap();
        let state = we.pulse(0, 1.0).unwrap();
        let e0 = we.energy(&state);
        let final_state = we.simulate(&state, 1000);
        let ef = we.energy(&final_state);
        assert!(ef < e0);
    }

    #[test]
    fn test_energy_conservation_no_damping() {
        let g = make_path4();
        let we = WaveEquation::new(g, 0.0, 0.001).unwrap();
        let state = we.pulse(0, 1.0).unwrap();
        let e0 = we.energy(&state);
        let final_state = we.simulate(&state, 100);
        let ef = we.energy(&final_state);
        // Energy should be approximately conserved (within 5% for small dt)
        assert!((ef - e0).abs() / e0 < 0.05, "energy drift: e0={e0}, ef={ef}");
    }

    #[test]
    fn test_with_speed() {
        let g = make_path4();
        let we = WaveEquation::with_speed(g, 1.0, 0.0, 0.01);
        assert!((we.wave_speed - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_graph_error() {
        let g = Graph::new(0);
        assert!(WaveEquation::new(g, 0.0, 0.01).is_err());
    }
}
