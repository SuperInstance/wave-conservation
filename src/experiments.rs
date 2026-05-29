use crate::wave::{WaveState, WaveReport};
use crate::spectral;

/// Build a path graph adjacency matrix
pub fn path_graph(n: usize) -> Vec<Vec<f64>> {
    let mut adj = vec![vec![0.0; n]; n];
    for i in 0..n.saturating_sub(1) {
        adj[i][i + 1] = 1.0;
        adj[i + 1][i] = 1.0;
    }
    adj
}

/// Build a cycle graph
pub fn cycle_graph(n: usize) -> Vec<Vec<f64>> {
    let mut adj = path_graph(n);
    if n > 2 {
        adj[0][n - 1] = 1.0;
        adj[n - 1][0] = 1.0;
    }
    adj
}

/// Build a complete graph
pub fn complete_graph(n: usize) -> Vec<Vec<f64>> {
    let mut adj = vec![vec![1.0; n]; n];
    for i in 0..n {
        adj[i][i] = 0.0;
    }
    adj
}

/// Build a star graph (node 0 = center)
pub fn star_graph(n: usize) -> Vec<Vec<f64>> {
    let mut adj = vec![vec![0.0; n]; n];
    for i in 1..n {
        adj[0][i] = 1.0;
        adj[i][0] = 1.0;
    }
    adj
}

/// Build a barbell graph (two cliques of size k joined by a bridge)
pub fn barbell_graph(k: usize) -> Vec<Vec<f64>> {
    let n = 2 * k;
    let mut adj = vec![vec![0.0; n]; n];
    // Left clique
    for i in 0..k {
        for j in 0..k {
            if i != j { adj[i][j] = 1.0; }
        }
    }
    // Right clique
    for i in k..n {
        for j in k..n {
            if i != j { adj[i][j] = 1.0; }
        }
    }
    // Bridge
    adj[k - 1][k] = 1.0;
    adj[k][k - 1] = 1.0;
    adj
}

/// Experiment 1: Verify wave speed = √λ₂
/// Start pulse at node 0, measure when peak reaches midpoint
pub fn verify_wave_speed(adj: &[Vec<f64>]) -> WaveReport {
    let n = adj.len();
    let mid = n / 2;
    let dt = 0.01;
    let max_steps = 20000;

    let eigenvals = spectral::eigenvalues(adj);
    let lambda2 = if eigenvals.len() > 1 { eigenvals[1] } else { 0.0 };
    let predicted_speed = lambda2.sqrt();

    let mut state = WaveState::new(adj.to_vec()).with_damping(0.001);
    state.pulse(0, 1.0);
    let initial_energy = state.energy();

    let mut energy_at_step = Vec::new();
    let mut coherence_at_step = Vec::new();
    let mut wave_speed = 0.0;
    let mut speed_measured = false;

    let mut max_disp_at_mid = 0.0;
    let mut mid_peak_step = 0usize;

    for step in 0..max_steps {
        state.step(dt);

        if step % 50 == 0 {
            energy_at_step.push(state.energy());
            coherence_at_step.push(state.coherence());
        }

        // Track when peak displacement reaches midpoint
        let disp_mid = state.displacement[mid].abs();
        if disp_mid > max_disp_at_mid {
            max_disp_at_mid = disp_mid;
            mid_peak_step = step;
            if !speed_measured && disp_mid > 0.01 {
                // First time we see significant displacement at midpoint
                let time_elapsed = step as f64 * dt;
                let distance = mid as f64;
                if time_elapsed > 0.0 {
                    wave_speed = distance / time_elapsed;
                    speed_measured = true;
                }
            }
        }
    }

    // Find coherence halflife
    let initial_coherence = coherence_at_step.first().copied().unwrap_or(1.0);
    let coherence_halflife = coherence_at_step.iter().position(|&c| {
        if initial_coherence > 0.0 { c < initial_coherence * 0.5 } else { false }
    }).map(|i| i * 50).unwrap_or(max_steps);

    // Conservation ratio
    let final_energy = state.energy();
    let cr = if initial_energy > 0.0 { final_energy / initial_energy } else { 0.0 };

    let speed_error = (wave_speed - predicted_speed).abs();

    WaveReport {
        initial_energy,
        energy_at_step,
        coherence_at_step,
        wave_speed,
        predicted_speed,
        speed_error,
        coherence_halflife,
        cr,
    }
}

/// Experiment 2: CR vs coherence halflife across different graphs
pub fn cr_vs_coherence() -> Vec<(f64, usize)> {
    let mut results = Vec::new();

    let graphs: Vec<(&str, Vec<Vec<f64>>)> = vec![
        ("path10", path_graph(10)),
        ("cycle10", cycle_graph(10)),
        ("star10", star_graph(10)),
        ("complete10", complete_graph(10)),
        ("barbell5", barbell_graph(5)),
        ("path20", path_graph(20)),
        ("cycle20", cycle_graph(20)),
    ];

    for (_name, adj) in &graphs {
        let report = verify_wave_speed(adj);
        results.push((report.cr, report.coherence_halflife));
    }

    results
}

/// Experiment 3: Standing waves — drive at eigenvalue frequency
pub fn standing_waves(adj: &[Vec<f64>], freq: f64, steps: usize) -> Vec<f64> {
    let n = adj.len();
    let dt = 0.005;
    let mut state = WaveState::new(adj.to_vec()).with_damping(0.001);

    for step in 0..steps {
        // Drive node 0 at the given frequency
        let t = step as f64 * dt;
        state.displacement[0] = (freq * t).sin() * 0.5;
        state.step(dt);
    }

    state.displacement.clone()
}

/// Experiment 4: Fiedler reflection — start wave on one side of Fiedler cut
pub fn fiedler_reflection(adj: &[Vec<f64>]) -> WaveReport {
    let n = adj.len();
    let fiedler = spectral::fiedler_vector(adj);

    // Start pulse on the side where Fiedler value is negative
    let pulse_node = fiedler.iter().enumerate()
        .filter(|(_, v)| **v < 0.0)
        .map(|(i, _)| i)
        .next()
        .unwrap_or(0);

    let dt = 0.01;
    let max_steps = 15000;

    let eigenvals = spectral::eigenvalues(adj);
    let lambda2 = if eigenvals.len() > 1 { eigenvals[1] } else { 0.0 };
    let predicted_speed = lambda2.sqrt();

    let mut state = WaveState::new(adj.to_vec()).with_damping(0.001);
    state.pulse(pulse_node, 1.0);
    let initial_energy = state.energy();

    let mut energy_at_step = Vec::new();
    let mut coherence_at_step = Vec::new();

    for step in 0..max_steps {
        state.step(dt);
        if step % 50 == 0 {
            energy_at_step.push(state.energy());
            coherence_at_step.push(state.coherence());
        }
    }

    let initial_coherence = coherence_at_step.first().copied().unwrap_or(1.0);
    let coherence_halflife = coherence_at_step.iter().position(|&c| {
        if initial_coherence > 0.0 { c < initial_coherence * 0.5 } else { false }
    }).map(|i| i * 50).unwrap_or(max_steps);

    let final_energy = state.energy();
    let cr = if initial_energy > 0.0 { final_energy / initial_energy } else { 0.0 };

    WaveReport {
        initial_energy,
        energy_at_step,
        coherence_at_step,
        wave_speed: 0.0,
        predicted_speed,
        speed_error: 0.0,
        coherence_halflife,
        cr,
    }
}

/// Experiment 5: Interference — two pulses from opposite ends
pub fn interference_pattern(adj: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = adj.len();
    let dt = 0.01;
    let steps = 5000;
    let mut state = WaveState::new(adj.to_vec()).with_damping(0.0005);

    state.pulse(0, 1.0);
    if n > 1 {
        state.pulse(n - 1, 1.0);
    }

    let mut pattern = Vec::new();
    for step in 0..steps {
        state.step(dt);
        if step % 50 == 0 {
            pattern.push(state.displacement.clone());
        }
    }
    pattern
}
