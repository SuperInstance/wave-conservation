//! Integration tests for wave-conservation

use wave_conservation::*;

#[test]
fn test_wave_pulse() {
    let adj = vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]];
    let mut w = WaveState::new(adj);
    w.pulse(0, 1.0);
    assert!((w.displacement[0] - 1.0).abs() < 1e-10);
    assert!((w.displacement[1]).abs() < 1e-10);
}

#[test]
fn test_wave_step_propagates() {
    let adj = vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]];
    let mut w = WaveState::new(adj);
    w.pulse(0, 1.0);
    let d0_initial = w.displacement[0];
    w.step(0.1);
    assert!((w.displacement[0] - d0_initial).abs() > 1e-10, "wave should propagate");
}

#[test]
fn test_wave_damping_reduces_energy() {
    let adj = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
    let mut w = WaveState::new(adj).with_damping(0.5);
    w.pulse(0, 10.0);
    let e0: f64 = w.displacement.iter().map(|x| x * x).sum();
    for _ in 0..100 {
        w.step(0.01);
    }
    let ef: f64 = w.displacement.iter().map(|x| x * x).sum();
    assert!(ef < e0, "damping should reduce energy");
}

#[test]
fn test_spectral_eigenvalues_path() {
    let adj = vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]];
    let eigs = eigenvalues(&adj);
    assert_eq!(eigs.len(), 3);
    assert!(eigs[0].abs() < 1.0, "smallest eigenvalue should be near 0: {}", eigs[0]);
}

#[test]
fn test_wave_node_count() {
    let adj = vec![vec![0.0; 4]; 4];
    let w = WaveState::new(adj);
    assert_eq!(w.n(), 4);
}

#[test]
fn test_spectral_fiedler_vector() {
    let adj = vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]];
    let fiedler = spectral::fiedler_vector(&adj);
    assert_eq!(fiedler.len(), 3);
}

#[test]
fn test_spectral_conservation_ratio() {
    let adj = vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]];
    let cr = spectral::conservation_ratio(&adj);
    assert!(cr >= 0.0 && cr <= 1.0, "CR should be in [0,1]: {cr}");
}
