mod wave;
mod experiments;
mod spectral;

fn main() {
    println!("=== Wave Conservation: Spectral Wave Propagation ===\n");

    // Demo: path graph
    let n = 20;
    let adj = experiments::path_graph(n);
    println!("Path graph ({} nodes)", n);

    // Experiment 1: Wave speed
    let report = experiments::verify_wave_speed(&adj);
    println!("\n--- Wave Speed Verification ---");
    println!("Predicted speed (√λ₂): {:.6}", report.predicted_speed);
    println!("Measured speed:         {:.6}", report.wave_speed);
    println!("Error:                  {:.6} ({:.2}%)",
        report.speed_error,
        if report.predicted_speed > 0.0 { 100.0 * report.speed_error / report.predicted_speed } else { 0.0 });
    println!("Conservation Ratio:     {:.6}", report.cr);
    println!("Coherence halflife:     {} steps", report.coherence_halflife);

    // Experiment 2: CR vs coherence
    println!("\n--- CR vs Coherence Halflife ---");
    let cr_data = experiments::cr_vs_coherence();
    for (cr, hl) in &cr_data {
        println!("CR={:.4}  coherence_halflife={}", cr, hl);
    }

    // Experiment 3: Standing waves
    println!("\n--- Standing Wave at λ₂ frequency ---");
    let eigenvals = spectral::eigenvalues(&adj);
    if eigenvals.len() > 1 {
        let freq = eigenvals[1].sqrt();
        let sw = experiments::standing_waves(&adj, freq, 500);
        let max_amp = sw.iter().cloned().fold(0.0_f64, f64::max);
        println!("Driving at √λ₂ = {:.6}, max response amplitude: {:.6}", freq, max_amp);
    }

    // Experiment 4: Fiedler reflection
    println!("\n--- Fiedler Reflection ---");
    let f_report = experiments::fiedler_reflection(&adj);
    println!("Energy reflected back, coherence halflife: {} steps", f_report.coherence_halflife);

    // Experiment 5: Frequency sweep
    println!("\n--- Frequency Sweep (eigenvalue spectrum from wave response) ---");
    let sweep = spectral::frequency_sweep(&adj, 0.0, 4.0, 200);
    let peaks = spectral::find_peaks(&sweep, 0.5);
    println!("Response peaks at frequencies: ");
    for (freq, resp) in &peaks {
        println!("  f={:.4}  response={:.4}  (λ ≈ {:.4})", freq, resp, freq * freq);
    }
    println!("\nActual eigenvalues: ");
    for (i, &e) in eigenvals.iter().enumerate().take(6) {
        println!("  λ[{}] = {:.6}", i, e);
    }

    // Experiment 6: Interference
    println!("\n--- Wave Interference ---");
    let pattern = experiments::interference_pattern(&adj);
    let max_interference = pattern.iter().map(|row| row.iter().cloned().fold(0.0_f64, f64::max)).fold(0.0_f64, f64::max);
    println!("Max interference amplitude: {:.6}", max_interference);

    println!("\n=== All experiments complete ===");
}
