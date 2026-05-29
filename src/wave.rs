/// WaveSimulation: wave equation on graphs
/// d²u/dt² = -L·u - γ·du/dt  where L is the graph Laplacian

#[derive(Clone)]
pub struct WaveState {
    pub displacement: Vec<f64>,
    pub velocity: Vec<f64>,
    adjacency: Vec<Vec<f64>>,
    pub damping: f64,
}

impl WaveState {
    pub fn new(adj: Vec<Vec<f64>>) -> WaveState {
        let n = adj.len();
        WaveState {
            displacement: vec![0.0; n],
            velocity: vec![0.0; n],
            adjacency: adj,
            damping: 0.0,
        }
    }

    pub fn with_damping(mut self, d: f64) -> Self {
        self.damping = d;
        self
    }

    pub fn n(&self) -> usize {
        self.displacement.len()
    }

    pub fn pulse(&mut self, node: usize, amplitude: f64) {
        if node < self.displacement.len() {
            self.displacement[node] += amplitude;
        }
    }

    /// Laplacian applied to displacement: L·u where L = D - A
    fn laplacian_times(&self, u: &[f64]) -> Vec<f64> {
        let n = self.adjacency.len();
        let mut result = vec![0.0; n];
        for i in 0..n {
            let mut degree = 0.0;
            for j in 0..n {
                let w = self.adjacency[i][j];
                degree += w;
                result[i] -= w * u[j];
            }
            result[i] += degree * u[i];
        }
        result
    }

    /// Step forward using velocity Verlet (symplectic) integration
    pub fn step(&mut self, dt: f64) {
        let n = self.displacement.len();

        // acceleration = -L·u - damping·v
        let compute_acc = |disp: &[f64], vel: &[f64]| -> Vec<f64> {
            let mut lu = vec![0.0; n];
            for i in 0..n {
                let mut degree = 0.0;
                for j in 0..n {
                    let w = self.adjacency[i][j];
                    degree += w;
                    lu[i] -= w * disp[j];
                }
                lu[i] += degree * disp[i];
            }
            let mut a = vec![0.0; n];
            for i in 0..n {
                a[i] = -lu[i] - self.damping * vel[i];
            }
            a
        };

        // Half-step velocity
        let a0 = compute_acc(&self.displacement, &self.velocity);
        for i in 0..n {
            self.velocity[i] += 0.5 * a0[i] * dt;
        }
        // Full-step position
        for i in 0..n {
            self.displacement[i] += self.velocity[i] * dt;
        }
        // Half-step velocity with new position
        let a1 = compute_acc(&self.displacement, &self.velocity);
        for i in 0..n {
            self.velocity[i] += 0.5 * a1[i] * dt;
        }
    }

    /// Total kinetic + potential energy
    pub fn energy(&self) -> f64 {
        let n = self.displacement.len();
        // kinetic: 0.5 * sum(v_i^2)
        let kinetic: f64 = self.velocity.iter().map(|v| 0.5 * v * v).sum();
        // potential: 0.5 * sum_{i,j} A_ij (u_i - u_j)^2 = 0.5 * u^T L u
        let mut potential = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                let d = self.displacement[i] - self.displacement[j];
                potential += self.adjacency[i][j] * d * d;
            }
        }
        kinetic + potential
    }

    /// Coherence: average correlation of displacement between adjacent nodes
    pub fn coherence(&self) -> f64 {
        let n = self.displacement.len();
        if n == 0 { return 0.0; }

        // Use normalized correlation-like measure
        let mean: f64 = self.displacement.iter().sum::<f64>() / n as f64;
        let variance: f64 = self.displacement.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        if variance < 1e-15 { return 1.0; }

        // Average correlation between neighbors
        let mut total_corr = 0.0;
        let mut count = 0;
        for i in 0..n {
            for j in 0..n {
                if self.adjacency[i][j] > 0.0 && i < j {
                    let cov = (self.displacement[i] - mean) * (self.displacement[j] - mean) / variance;
                    total_corr += cov;
                    count += 1;
                }
            }
        }
        if count == 0 { 1.0 } else { total_corr / count as f64 }
    }

    /// Conservation ratio: E_now / E_initial
    pub fn conservation_ratio(&self, initial_energy: f64) -> f64 {
        if initial_energy <= 0.0 { return 0.0; }
        self.energy() / initial_energy
    }
}

#[derive(Clone, Debug)]
pub struct WaveReport {
    pub initial_energy: f64,
    pub energy_at_step: Vec<f64>,
    pub coherence_at_step: Vec<f64>,
    pub wave_speed: f64,
    pub predicted_speed: f64,
    pub speed_error: f64,
    pub coherence_halflife: usize,
    pub cr: f64,
}
