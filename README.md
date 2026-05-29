# wave-conservation

**Spectral wave propagation on graphs — wave speed = √λ₂, standing waves reveal the eigenvalue spectrum, conservation ratio predicts coherence.**

Pure Rust implementation of the discrete wave equation on graph structures. Propagates waves through path, cycle, complete, star, and barbell graphs, verifying that graph spectral properties (eigenvalues, Fiedler vector) directly govern wave behavior.

## What This Gives You

- **Wave speed verification** — confirms that wave speed on graphs equals √λ₂ (Fiedler eigenvalue)
- **Standing wave detection** — driving at eigenfrequencies produces resonance peaks
- **CR vs coherence** — maps conservation ratio to wave coherence halflife
- **Frequency sweep** — recovers the full eigenvalue spectrum from wave response alone
- **Fiedler reflection** — waves reflected at the algebraic connectivity boundary
- **Wave interference** — constructive/destructive patterns on graph structures
- **Zero dependencies** — pure Rust, no external math crates

## Quick Start

```rust
use wave_conservation::wave::WaveState;
use wave_conservation::spectral;

let adj = vec![vec![/* adjacency matrix */]];
let mut wave = WaveState::new(adj).with_damping(0.01);
wave.pulse(0, 1.0);  // inject at node 0

for _ in 0..1000 {
    wave.step(0.01);
}
```

```bash
cargo run  # runs all 6 experiments on a path graph
```

## Experiments

| # | Experiment | What It Shows |
|---|-----------|---------------|
| 1 | Wave Speed | Measured speed ≈ √λ₂ (predicted from spectrum) |
| 2 | CR vs Coherence | Higher conservation ratio → longer coherence halflife |
| 3 | Standing Waves | Resonance at λ₂ frequency confirms eigenvalue spectrum |
| 4 | Fiedler Reflection | Energy reflected at the algebraic connectivity boundary |
| 5 | Frequency Sweep | Response peaks recover all eigenvalues |
| 6 | Wave Interference | Multi-source interference patterns |

## API Reference

### `WaveState`

```rust
let mut ws = WaveState::new(adjacency_matrix);
ws.pulse(node_index, amplitude);  // inject energy
ws.step(dt);                       // advance one timestep (velocity Verlet)
```

### `spectral`

```rust
let eigs = spectral::eigenvalues(&adj);         // all eigenvalues (power iteration)
let response = spectral::frequency_sweep(&adj, f_min, f_max, steps);
let peaks = spectral::find_peaks(&response, threshold);
```

### Graph generators (`experiments`)

```rust
let adj = experiments::path_graph(n);
let adj = experiments::cycle_graph(n);
let adj = experiments::complete_graph(n);
let adj = experiments::star_graph(n);
let adj = experiments::barbell_graph(k);
```

## Testing

```bash
cargo test
cargo run  # interactive demo with all experiments
```

## Installation

```toml
[dependencies]
wave-conservation = { git = "https://github.com/SuperInstance/wave-conservation" }
```

## How It Fits

Part of the SuperInstance ecosystem:

- **[heat-spectral](https://github.com/SuperInstance/heat-spectral)** — Heat diffusion (parabolic PDE) on graphs
- **wave-conservation** — Wave propagation (hyperbolic PDE) on graphs (this repo)
- **[graph-neural](https://github.com/SuperInstance/graph-neural)** — Graph neural network spectral primitives

## License

MIT
