# EMMoMSuite Validation Report — `G1_sphere_CFIE_Mie`

| | |
|---|---|
| solver | EMMoMSuite v0.3.1 (Julia 1.12.3) |
| generated | 2026-10-10 11:51:29 |
| pipeline | geometry → gmsh mesh → MoM solve → RCS / far-field → this report |

## 1 · Geometry & Mesh

![geometry views](geometry_views.png)

![mesh views](mesh_views.png)

Boundary discretization: **790 triangles**, 397 nodes; geometry source `cases/geo/sphere_r0p5.geo` (gmsh/OpenCASCADE).

## 2 · Simulation Configuration

| item | value |
|---|---|
| geometry | `cases/geo/sphere_r0p5.geo` |
| mesh size | 0.1 m |
| mesh dim | 2 (surface triangulation) |
| frequency | 300.0 MHz (λ = 1.0 m) |
| formulation | CFIE (α = 0.5) |
| interfaces (n̂ = mesh triangle normal) | plus side (n̂) | minus side |
|---|---|---|
| `body` (closed) | air | PEC |
| incidence | θᵢ = 90.0°, φᵢ = 180.0°, pol = [0.0, 0.0, 1.0] |
| observation | θ ∈ [0°, 180°] (181 samples), φ cuts = 0.0°, 90.0° |
| reference | analytic Mie series, PEC sphere r = 0.5 m |

## 3 · Performance & Efficiency

| triangles/tets | nodes | unknowns | t_mesh | t_assemble | t_solve | t_RCS | t_plots |
|---|---|---|---|---|---|---|---|
| 790 | 397 | 1185 | 235.9 s | 3.1 s | 0.8 s | 1.2 s | 13.3 s |

| total time | assembly rate | LU throughput |
|---|---|---|
| 254.3 s | 0.000458 × 10⁹ interactions/s | 1.3 Gflop/s |

*Assembly rate counts N² impedance-matrix interactions; LU throughput uses the dense (2/3)·N³ flop model (single node, default BLAS threads).*
Machine-readable: `perf.csv`, `rcs.csv`.

## 4 · Results & Comparison

Bistatic RCS — MoM vs analytic reference:

![RCS cuts](rcs_cuts.png)

Normalized far-field |E| pattern:

![Far-field polar](farfield_polar.png)

Surface-current magnitude |J| (dB, normalized to peak):

![Current distribution](current_views.png)

| phi cut | RMSE vs Mie [dB] | verdict |
|---|---|---|
| 0.0° | 0.162 | pass |
| 90.0° | 0.137 | pass |

## 5 · Key Conclusions

- Electric resolution: mesh size 0.1 m = 0.1 λ at f = 300.0 MHz (90.0°, 180.0° plane-wave incidence).
- Accuracy: worst-cut RMSE vs analytic Mie reference = 0.162 dB over 2 cuts → **PASS** against the 0.5 dB acceptance line.
- Throughput: impedance assembly 0.000458 G-interactions/s, dense LU 1.3 Gflop/s at N = 1185 unknowns.
- Dominant cost: gmsh meshing — 235.9 s (93.0% of the 254.3 s total).
- Bistatic RCS dynamic range over the observed cuts: -6.7 … 1.5 dBsm.

## Artifacts

- geometry: `geometry_views.png` · RCS data: `rcs.csv` · performance: `perf.csv`
- plots: `rcs_cuts.png`, `farfield_polar.png`
