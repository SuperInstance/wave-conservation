# wave-conservation

**Spectral wave propagation on graphs — wave speed = √λ₂, standing waves reveal the eigenvalue spectrum, and the conservation ratio predicts how long coherence survives.**

> *Send a wave through a graph. Watch how fast it travels, how long it stays coherent, where it resonates. The wave doesn't know about eigenvalues — but it computes them anyway. The spectrum hides in the wave. You just have to know where to look.*

---

## The Problem

Graphs have spectral structure encoded in their Laplacian eigenvalues. But eigenvalues are abstract — you compute them with linear algebra, not physics. What if you could *measure* the spectrum by sending waves through the graph and watching what happens?

It turns out you can. The wave equation on a graph `d²u/dt² = -L·u` (where L is the graph Laplacian) produces wave behavior that is *governed entirely by the eigenvalue spectrum*. The wave speed equals √λ₂ (the Fiedler eigenvalue). The resonance frequencies are √λᵢ for each eigenvalue. And the conservation ratio λ₂/λₙ predicts how long wave coherence survives.

This crate implements the discrete wave equation on graphs and provides experiments that verify these spectral-wave correspondences. No external math dependencies — everything from eigenvalue computation to wave simulation is pure Rust.

## The Key Insight

**The wave equation on a graph IS a spectral probe. You don't need to compute eigenvalues separately — just send a wave and read off the spectrum from the response.**

Here's why. The Laplacian L has eigenvectors φ₁, φ₂, ..., φₙ with eigenvalues 0 = λ₁ < λ₂ ≤ ... ≤ λₙ. Any displacement u can be written as a superposition of eigenvectors: u = Σ cᵢφᵢ. Under the wave equation, each mode oscillates independently at frequency √λᵢ. So:

1. **Wave speed** — The lowest non-trivial mode (Fiedler) travels at speed √λ₂. Measure how fast a pulse crosses the graph, and you've measured the Fiedler eigenvalue.

2. **Standing waves** — Drive the graph at frequency ω = √λᵢ, and mode i resonates. A frequency sweep reveals all eigenvalues as peaks.

3. **Coherence** — The conservation ratio CR = λ₂/λₙ bounds how long waves stay coherent. High CR (well-connected graph) → long coherence. Low CR (path graph) → fast decoherence.

4. **Fiedler reflection** — The Fiedler vector partitions the graph into two communities. Waves launched from one side reflect at the boundary, like light at an interface.

The velocity Verlet integrator ensures symplectic energy conservation — the simulation respects the physics, not just the numerics.

## Architecture

```
 ┌──────────────────────────────────────────────────────┐
 │                  WaveState                           │
 │  d²u/dt² = -L·u - γ·du/dt                          │
 │                                                      │
 │  displacement[N] ──▶ velocity Verlet ──▶ step(dt)   │
 │  velocity[N]              integrator                │
 │  adjacency[N×N]              │                      │
 │  damping γ                   │                      │
 └──────────┬───────────────────┼──────────────────────┘
            │                   │
            ▼                   ▼
 ┌──────────────────┐  ┌─────────────────────────┐
 │    spectral      │  │     experiments          │
 │                  │  │                          │
 │ eigenvalues()    │  │ verify_wave_speed()     │
 │ fiedler_vector() │  │ cr_vs_coherence()       │
 │ conservation_    │  │ standing_waves()         │
 │   ratio()        │  │ fiedler_reflection()    │
 │ resonance_       │  │ interference_pattern()  │
 │   frequencies()  │  │                          │
 │ frequency_       │  │ Graph generators:       │
 │   sweep()        │  │  path, cycle, complete, │
 │ find_peaks()     │  │  star, barbell          │
 └──────────────────┘  └─────────────────────────┘

 Flow: adjacency → WaveState → simulate → measure speed/coherence/energy
         ↓
      spectral → eigenvalues → predict speed, find resonances
         ↓
      compare predicted vs measured → conservation verified
```

## Quick Start

```rust
use wave_conservation::wave::WaveState;
use wave_conservation::spectral;
use wave_conservation::experiments;

// Build a path graph
let adj = experiments::path_graph(20);

// Compute eigenvalues (power iteration, no external deps)
let eigs = spectral::eigenvalues(&adj);
println!("λ₂ = {:.4}", eigs[1]); // Fiedler eigenvalue
println!("Wave speed = √λ₂ = {:.4}", eigs[1].sqrt());

// Launch a wave
let mut wave = WaveState::new(adj).with_damping(0.001);
wave.pulse(0, 1.0);
let e0 = wave.energy();

for _ in 0..1000 {
    wave.step(0.01);
}

println!("Energy conservation: {:.4}", wave.conservation_ratio(e0));
```

## Tutorial

### 1. Wave Propagation on a Path Graph

Send a pulse from one end and watch it travel:

```rust
use wave_conservation::wave::WaveState;
use wave_conservation::experiments;

let adj = experiments::path_graph(20);
let mut wave = WaveState::new(adj);
wave.pulse(0, 1.0); // inject energy at node 0

// Step forward and track displacement
for step in 0..5000 {
    wave.step(0.01);
    if step % 500 == 0 {
        println!("Step {}: node[10] = {:.4}", step, wave.displacement[10]);
    }
}
```

### 2. Wave Speed = √λ₂ (Fiedler Eigenvalue)

The most fundamental correspondence — verify it experimentally:

```rust
use wave_conservation::experiments;

let adj = experiments::path_graph(30);
let report = experiments::verify_wave_speed(&adj);

println!("Predicted speed (√λ₂): {:.4}", report.predicted_speed);
println!("Measured speed:          {:.4}", report.wave_speed);
println!("Error: {:.4} ({:.1}%)",
    report.speed_error,
    if report.predicted_speed > 0.0 {
        100.0 * report.speed_error / report.predicted_speed
    } else { 0.0 }
);
```

### 3. Standing Waves Reveal the Spectrum

Drive at eigenvalue frequencies and observe resonance:

```rust
use wave_conservation::{spectral, experiments};

let adj = experiments::path_graph(10);
let eigs = spectral::eigenvalues(&adj);

// Drive at √λ₂ — should produce resonance
let freq = eigs[1].sqrt();
let response = experiments::standing_waves(&adj, freq, 10000);
let max_amp = response.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
println!("Resonance at √λ₂={:.4}: max amplitude = {:.4}", freq, max_amp);
```

### 4. Frequency Sweep — Eigenvalues from Wave Response

Sweep frequencies and find peaks that correspond to eigenvalues:

```rust
use wave_conservation::{spectral, experiments};

let adj = experiments::path_graph(10);
let eigs = spectral::eigenvalues(&adj);

// Sweep from 0.1 to 3.0 Hz
let sweep = spectral::frequency_sweep(&adj, 0.1, 3.0, 200);
let peaks = spectral::find_peaks(&sweep, 0.3);

println!("Frequency sweep peaks (should match √λᵢ):");
for (freq, resp) in &peaks {
    println!("  f={:.4}, response={:.4} (λ ≈ {:.4})", freq, resp, freq * freq);
}

println!("\nActual eigenvalues:");
for (i, &e) in eigs.iter().enumerate().take(5) {
    println!("  λ[{}] = {:.4}, √λ = {:.4}", i, e, e.sqrt());
}
```

### 5. Conservation Ratio Predicts Coherence

Higher CR → longer coherence halflife:

```rust
use wave_conservation::experiments;

let data = experiments::cr_vs_coherence();
println!("CR vs Coherence Halflife:");
for (cr, halflife) in &data {
    println!("  CR={:.4}  halflife={} steps", cr, halflife);
}
// Complete graph: CR ≈ 1.0, longest coherence
// Path graph: CR ≈ 0.1, shortest coherence
```

### 6. Energy Conservation (Symplectic Integrator)

The velocity Verlet integrator preserves energy:

```rust
use wave_conservation::wave::WaveState;
use wave_conservation::experiments;

let adj = experiments::path_graph(20);
let mut wave = WaveState::new(adj).with_damping(0.0); // no damping!
wave.pulse(0, 1.0);
let e0 = wave.energy();

for _ in 0..10000 {
    wave.step(0.002); // small dt for good conservation
}

let ef = wave.energy();
println!("Energy drift: {:.4}% ({:.6} → {:.6})",
    100.0 * (ef - e0) / e0, e0, ef);
```

### 7. Graph Zoo — Different Topologies, Different Spectra

```rust
use wave_conservation::{spectral, experiments};

let graphs = [
    ("Path(10)", experiments::path_graph(10)),
    ("Cycle(10)", experiments::cycle_graph(10)),
    ("Star(10)", experiments::star_graph(10)),
    ("Complete(10)", experiments::complete_graph(10)),
    ("Barbell(5)", experiments::barbell_graph(5)),
];

for (name, adj) in &graphs {
    let eigs = spectral::eigenvalues(adj);
    let cr = spectral::conservation_ratio(adj);
    println!("{}: λ₂={:.4}, λₙ={:.4}, CR={:.4}",
        name, eigs[1], eigs[eigs.len()-1], cr);
}
```

## API Reference

### `wave` Module

| Type/Method | Description |
|-------------|-------------|
| `WaveState::new(adj)` | Create wave state from adjacency matrix |
| `.with_damping(γ)` | Set damping coefficient (builder pattern) |
| `.pulse(node, amplitude)` | Inject energy at a node |
| `.step(dt)` | Advance one timestep (velocity Verlet) |
| `.energy()` | Total kinetic + potential energy |
| `.coherence()` | Average correlation of displacement between neighbors |
| `.conservation_ratio(e₀)` | Ratio E_now / E_initial |
| `WaveReport` | Full experiment results (speed, coherence, energy) |

### `spectral` Module

| Function | Description |
|----------|-------------|
| `eigenvalues(&adj)` | All eigenvalues via power iteration + deflation |
| `fiedler_vector(&adj)` | Eigenvector of λ₂ (graph partition) |
| `conservation_ratio(&adj)` | CR = λ₂/λₙ |
| `resonance_frequencies(&adj)` | √λᵢ for each eigenvalue |
| `frequency_response(&adj, freq)` | Steady-state amplitude at given frequency |
| `frequency_sweep(&adj, min, max, steps)` | Response curve over frequency range |
| `find_peaks(&sweep, threshold)` | Locate resonance peaks in sweep data |

### `experiments` Module

| Function | Description |
|----------|-------------|
| `path_graph(n)` | Path graph adjacency |
| `cycle_graph(n)` | Cycle graph adjacency |
| `complete_graph(n)` | Complete graph adjacency |
| `star_graph(n)` | Star graph adjacency |
| `barbell_graph(k)` | Two k-cliques joined by a bridge |
| `verify_wave_speed(&adj)` | Experiment 1: measure vs predicted speed |
| `cr_vs_coherence()` | Experiment 2: CR vs coherence across graph types |
| `standing_waves(&adj, freq, steps)` | Experiment 3: drive at eigenfrequency |
| `fiedler_reflection(&adj)` | Experiment 4: wave at Fiedler boundary |
| `interference_pattern(&adj)` | Experiment 5: two-source interference |

## Ecosystem Role

`wave-conservation` is the **wave dynamics layer** of the SuperInstance ecosystem:

- **[heat-spectral](https://github.com/SuperInstance/heat-spectral)** — Heat diffusion (parabolic PDE) on graphs
- **[analog-spectral](https://github.com/SuperInstance/analog-spectral)** — Physical dials and eigenvalue estimation
- **[dial-ecology](https://github.com/SuperInstance/dial-ecology)** — Lotka-Volterra competition for traditions
- **wave-conservation** — Wave propagation and spectral analysis on graphs (this crate)

The conservation ratio computed here connects to the spectral thermostat in `analog-spectral`, while the graph eigenvalue analysis provides structural insights for `dial-ecology`'s niche-based competition models.

## Installation

```toml
[dependencies]
wave-conservation = { git = "https://github.com/SuperInstance/wave-conservation" }
```

## License

MIT
