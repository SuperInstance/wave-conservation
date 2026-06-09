//! Tutorial: wave-conservation — spectral wave propagation on graphs.
//!
//! Run with: `cargo run --example tutorial`

use wave_conservation::wave::WaveState;
use wave_conservation::{spectral, experiments};

fn main() {
    println!("=== wave-conservation Tutorial ===\n");

    // ---- 1. Wave Propagation on a Path Graph ----
    println!("--- 1. Wave Propagation ---");
    let adj = experiments::path_graph(20);
    let mut wave = WaveState::new(adj);
    wave.pulse(0, 1.0);
    let e0 = wave.energy();
    println!("Initial energy: {:.4}", e0);

    for step in 0..5000 {
        wave.step(0.01);
        if step % 1000 == 0 {
            println!(
                "  Step {}: node[10]={:.4}, energy={:.4}",
                step, wave.displacement[10], wave.energy()
            );
        }
    }
    println!();

    // ---- 2. Wave Speed = √λ₂ ----
    println!("--- 2. Wave Speed Verification ---");
    let adj = experiments::path_graph(30);
    let eigs = spectral::eigenvalues(&adj);
    println!("λ₂ = {:.4}", eigs[1]);
    println!("Predicted wave speed (√λ₂): {:.4}", eigs[1].sqrt());

    let report = experiments::verify_wave_speed(&adj);
    println!("Measured wave speed:          {:.4}", report.wave_speed);
    println!(
        "Error: {:.4} ({:.1}%)",
        report.speed_error,
        if report.predicted_speed > 0.0 {
            100.0 * report.speed_error / report.predicted_speed
        } else {
            0.0
        }
    );
    println!();

    // ---- 3. Eigenvalues of Different Graphs ----
    println!("--- 3. Graph Zoo: Spectra ---");
    let graphs: Vec<(&str, Vec<Vec<f64>>)> = vec![
        ("Path(10)", experiments::path_graph(10)),
        ("Cycle(10)", experiments::cycle_graph(10)),
        ("Star(10)", experiments::star_graph(10)),
        ("Complete(10)", experiments::complete_graph(10)),
        ("Barbell(5)", experiments::barbell_graph(5)),
    ];

    for (name, adj) in &graphs {
        let eigs = spectral::eigenvalues(adj);
        let cr = spectral::conservation_ratio(adj);
        let res_freqs = spectral::resonance_frequencies(adj);
        println!(
            "{}: λ₂={:.4}, λₙ={:.4}, CR={:.4}, resonance_freqs=[{:.3}, {:.3}, ...]",
            name, eigs[1], eigs[eigs.len() - 1], cr,
            res_freqs.first().unwrap_or(&0.0),
            res_freqs.get(1).unwrap_or(&0.0)
        );
    }
    println!();

    // ---- 4. Standing Waves at Eigenfrequencies ----
    println!("--- 4. Standing Waves ---");
    let adj = experiments::path_graph(10);
    let eigs = spectral::eigenvalues(&adj);

    for i in 1..eigs.len().min(4) {
        let freq = eigs[i].sqrt();
        let response = experiments::standing_waves(&adj, freq, 500);
        let max_amp = response.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        println!(
            "  √λ[{}] = {:.4}: max amplitude = {:.4}",
            i, freq, max_amp
        );
    }
    println!();

    // ---- 5. Frequency Sweep — Eigenvalues from Wave Response ----
    println!("--- 5. Frequency Sweep ---");
    let adj = experiments::path_graph(10);
    let sweep = spectral::frequency_sweep(&adj, 0.1, 3.0, 100);
    let peaks = spectral::find_peaks(&sweep, 0.3);

    println!("Response peaks (frequency → eigenvalue estimate):");
    for (freq, resp) in peaks.iter().take(5) {
        println!("  f={:.4} → λ≈{:.4} (response={:.4})", freq, freq * freq, resp);
    }

    let eigs = spectral::eigenvalues(&adj);
    println!("\nActual eigenvalues:");
    for (i, &e) in eigs.iter().enumerate().take(5) {
        println!("  λ[{}] = {:.4}", i, e);
    }
    println!();

    // ---- 6. Energy Conservation (Symplectic Integrator) ----
    println!("--- 6. Energy Conservation ---");
    let adj = experiments::path_graph(20);
    let mut wave = WaveState::new(adj).with_damping(0.0);
    wave.pulse(0, 1.0);
    let e0 = wave.energy();

    for _ in 0..10000 {
        wave.step(0.002);
    }
    let ef = wave.energy();
    println!(
        "Undamped: E₀={:.6}, E_final={:.6}, drift={:.4}%",
        e0, ef, 100.0 * (ef - e0) / e0
    );

    // Now with damping
    let adj = experiments::path_graph(20);
    let mut wave = WaveState::new(adj).with_damping(0.1);
    wave.pulse(0, 1.0);
    let e0 = wave.energy();
    for _ in 0..5000 {
        wave.step(0.01);
    }
    let ef = wave.energy();
    println!(
        "Damped (γ=0.1): E₀={:.6}, E_final={:.6}, ratio={:.4}",
        e0, ef, ef / e0
    );
    println!();

    // ---- 7. Fiedler Vector and Graph Partitioning ----
    println!("--- 7. Fiedler Vector ---");
    let adj = experiments::barbell_graph(5);
    let fv = spectral::fiedler_vector(&adj);
    println!("Barbell(5) Fiedler vector:");
    for (i, v) in fv.iter().enumerate() {
        let side = if *v > 0.0 { "RIGHT" } else { "LEFT" };
        println!("  node[{}] = {:+.4} ({})", i, v, side);
    }

    // ---- 8. Conservation Ratio vs Coherence ----
    println!("\n--- 8. CR vs Coherence Halflife ---");
    let data = experiments::cr_vs_coherence();
    println!("{:<12} {:>8} {:>15}", "Graph", "CR", "Halflife");
    let names = ["path10", "cycle10", "star10", "complete10", "barbell5", "path20", "cycle20"];
    for (i, (cr, hl)) in data.iter().enumerate() {
        let name = names.get(i).unwrap_or(&"?");
        println!("{:<12} {:>8.4} {:>15}", name, cr, hl);
    }

    println!("\n=== Tutorial Complete ===");
}
