# EMMoMSuite Validation Report — `G1_sphere_CFIE_Mie`

| | |
|---|---|
| solver | EMMoMSuite v0.3.0 (Julia 1.12.3) |
| generated | 2026-10-09 23:09:15 |
| pipeline | geometry → gmsh mesh → MoM solve → RCS / far-field → this report |

## 1 · Geometry

![geometry views](geometry_views.png)

Boundary discretization: **380 triangles**, 192 nodes; geometry source `cases/geo/sphere_r0p5.geo` (gmsh/OpenCASCADE).

## 2 · Simulation Configuration

| item | value |
|---|---|
| geometry | `cases/geo/sphere_r0p5.geo` |
| mesh size | 0.15 m |
| mesh dim | 2 (surface triangulation) |
| frequency | 300.0 MHz (λ = 1.0 m) |
| formulation | CFIE (α = 0.5) |
| interfaces (n̂ = mesh triangle normal) | plus side (n̂) | minus side |
|---|---|---|
| `body` (closed) | Dielectric(εᵣ=1.0 + 0.0im, μᵣ=1.0 + 0.0im) | PEC() |
| incidence | θᵢ = 90.0°, φᵢ = 180.0°, pol = [0.0, 0.0, 1.0] |
| observation | θ ∈ [0°, 180°] (181 samples), φ cuts = 0.0°, 90.0° |
| reference | analytic Mie series, PEC sphere r = 0.5 m |

## 3 · Performance & Efficiency

| triangles/tets | nodes | unknowns | t_mesh | t_assemble | t_solve | t_RCS | t_plots |
|---|---|---|---|---|---|---|---|
| 380 | 192 | 570 | 203.8 s | 3.2 s | 1.4 s | 1.5 s | 7.1 s |

| total time | assembly rate | LU throughput |
|---|---|---|
| 217.0 s | 0.00 × 10⁹ interactions/s | 0.1 Gflop/s |
*Assembly rate counts N² impedance-matrix interactions; LU throughput uses the dense (2/3)·N³ flop model (single node, default BLAS threads).*
Machine-readable: `perf.csv`, `rcs.csv`.

## 4 · Results & Comparison

Bistatic RCS — MoM vs analytic reference:

![RCS cuts](rcs_cuts.png)

Normalized far-field |E| pattern:

![Far-field polar](farfield_polar.png)

| phi cut | RMSE vs Mie [dB] | verdict |
|---|---|---|
| 0.0° | 0.353 | pass |
| 90.0° | 0.297 | pass |

## 5 · Key Conclusions

- Electric resolution: mesh size 0.15 m = 0.15 λ at f = 300.0 MHz (90.0°, 180.0° plane-wave incidence).
- Accuracy: worst-cut RMSE vs analytic Mie reference = 0.353 dB over 2 cuts → **PASS** against the 0.5 dB acceptance line.
- Throughput: impedance assembly 0.0 G-interactions/s, dense LU 0.1 Gflop/s at N = 570 unknowns.
- Dominant cost: gmsh meshing — 203.8 s (94.0% of the 217.0 s total).
- Bistatic RCS dynamic range over the observed cuts: -6.9 … 1.4 dBsm.

## Artifacts

- geometry: `geometry_views.png` · RCS data: `rcs.csv` · performance: `perf.csv`
- plots: `rcs_cuts.png`, `farfield_polar.png`
