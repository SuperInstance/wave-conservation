//! Advanced: Complete spectral analysis pipeline — from graph topology to wave experiments.
//!
//! Demonstrates the full cycle: build graph → compute spectrum → verify wave
//! dynamics → sweep frequencies → analyze coherence and interference.
//!
//! Run with: `cargo run --example advanced`

use wave_conservation::{spectral, experiments};

fn main() {
    println!("=== Advanced: Complete Spectral Wave Analysis ===\n");

    // ---- Build a barbell graph (two cliques joined by a bridge) ----
    let k = 6;
    let adj = experiments::barbell_graph(k);
    let n = adj.len();
    println!("Barbell graph: {} nodes (two {}-cliques + bridge)\n", n, k);

    // ---- Full eigenvalue analysis ----
    println!("--- Eigenvalue Spectrum ---");
    let eigs = spectral::eigenvalues(&adj);
    let cr = spectral::conservation_ratio(&adj);
    println!("Eigenvalues ({}):", eigs.len());
    for (i, &e) in eigs.iter().enumerate() {
        let bar = "█".repeat((e * 3.0) as usize);
        println!("  λ[{:2}] = {:7.4} {}", i, e, bar);
    }
    println!("Conservation ratio (λ₂/λₙ): {:.4}", cr);
    println!();

    // ---- Fiedler vector partition ----
    println!("--- Fiedler Vector Partition ---");
    let fv = spectral::fiedler_vector(&adj);
    let left: Vec<usize> = fv.iter().enumerate()
        .filter(|(_, &v)| v < 0.0)
        .map(|(i, _)| i)
        .collect();
    let right: Vec<usize> = fv.iter().enumerate()
        .filter(|(_, &v)| v >= 0.0)
        .map(|(i, _)| i)
        .collect();
    println!("Left community:  nodes {:?}", left);
    println!("Right community: nodes {:?}", right);
    println!("Fiedler values: [{:.4} ... {:.4}]",
        fv.iter().cloned().fold(f64::INFINITY, f64::min),
        fv.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    println!();

    // ---- Resonance frequencies ----
    println!("--- Resonance Frequencies ---");
    let res_freqs = spectral::resonance_frequencies(&adj);
    for (i, &f) in res_freqs.iter().enumerate() {
        println!("  Mode {}: ω = {:.4} (λ = {:.4})", i + 1, f, f * f);
    }
    println!();

    // ---- Wave speed verification ----
    println!("--- Wave Speed Verification ---");
    let report = experiments::verify_wave_speed(&adj);
    println!("Predicted (√λ₂): {:.4}", report.predicted_speed);
    println!("Measured:          {:.4}", report.wave_speed);
    println!(
        "Error: {:.4} ({:.1}%)",
        report.speed_error,
        if report.predicted_speed > 0.0 {
            100.0 * report.speed_error / report.predicted_speed
        } else { 0.0 }
    );
    println!("Coherence halflife: {} steps", report.coherence_halflife);
    println!();

    // ---- Frequency sweep — find all eigenvalues from wave response ----
    println!("--- Frequency Sweep ---");
    let max_freq = eigs.last().unwrap_or(&1.0).sqrt() * 1.2;
    let sweep = spectral::frequency_sweep(&adj, 0.1, max_freq, 150);
    let peaks = spectral::find_peaks(&sweep, 0.2);

    println!("Discovered {} resonance peaks:", peaks.len());
    println!("{:<12} {:>12} {:>12}", "Frequency", "Response", "λ estimate");
    for (freq, resp) in &peaks {
        println!("{:<12.4} {:>12.4} {:>12.4}", freq, resp, freq * freq);
    }

    println!("\nActual eigenvalues for comparison:");
    for (i, &e) in eigs.iter().enumerate().take(6) {
        if e > 0.001 {
            println!("  λ[{}] = {:.4}, √λ = {:.4}", i, e, e.sqrt());
        }
    }
    println!();

    // ---- Fiedler reflection experiment ----
    println!("--- Fiedler Reflection ---");
    let f_report = experiments::fiedler_reflection(&adj);
    println!(
        "Wave launched from Fiedler-negative side: coherence halflife = {} steps",
        f_report.coherence_halflife
    );
    println!("Energy conservation: {:.4}", f_report.cr);
    println!();

    // ---- Interference pattern ----
    println!("--- Wave Interference ---");
    let pattern = experiments::interference_pattern(&adj);
    let max_amp = pattern.iter()
        .map(|row| row.iter().cloned().fold(0.0_f64, f64::max))
        .fold(0.0_f64, f64::max);
    let min_amp = pattern.iter()
        .flat_map(|row| row.iter())
        .cloned()
        .fold(0.0_f64, f64::min);

    println!("Interference from both ends: {} snapshots", pattern.len());
    println!("Max displacement: {:.4}, Min displacement: {:.4}", max_amp, min_amp);
    println!();

    // ---- Compare graphs: topology vs spectral properties ----
    println!("--- Graph Comparison ---");
    println!(
        "{:<15} {:>6} {:>8} {:>8} {:>8} {:>10}",
        "Graph", "N", "λ₂", "λₙ", "CR", "Half-life"
    );
    println!("{}", "─".repeat(60));

    let comparisons: Vec<(&str, Vec<Vec<f64>>)> = vec![
        ("Path(10)", experiments::path_graph(10)),
        ("Cycle(10)", experiments::cycle_graph(10)),
        ("Star(10)", experiments::star_graph(10)),
        ("Complete(10)", experiments::complete_graph(10)),
        ("Barbell(5)", experiments::barbell_graph(5)),
        ("Path(20)", experiments::path_graph(20)),
        ("Cycle(20)", experiments::cycle_graph(20)),
    ];

    for (name, adj) in &comparisons {
        let eigs = spectral::eigenvalues(adj);
        let cr = spectral::conservation_ratio(adj);
        let report = experiments::verify_wave_speed(adj);
        println!(
            "{:<15} {:>6} {:>8.4} {:>8.4} {:>8.4} {:>10}",
            name, adj.len(), eigs[1], eigs[eigs.len() - 1], cr,
            report.coherence_halflife
        );
    }

    println!("\n=== Analysis Complete ===");
}
