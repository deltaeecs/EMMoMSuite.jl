# EMMoMSuite Validation Report — `G6_ellipsoid_EFIE`

| | |
|---|---|
| solver | EMMoMSuite v0.3.1 (Julia 1.12.3) |
| generated | 2026-10-10 18:46:02 |
| pipeline | geometry → gmsh mesh → MoM solve → RCS / far-field → this report |
| parallel | MPI distributed GMRES (P = 2 ranks) |

## 1 · Geometry & Mesh

![geometry views](geometry_views.png)

![mesh views](mesh_views.png)

Boundary discretization: **1236 triangles**, 620 nodes; geometry source `cases/geo/ellipsoid.geo` (gmsh/OpenCASCADE).

## 2 · Simulation Configuration

| item | value |
|---|---|
| geometry | `cases/geo/ellipsoid.geo` |
| mesh size | 0.06 m |
| mesh dim | 2 (surface triangulation) |
| frequency | 1200.0 MHz (λ = 0.25 m) |
| formulation | EFIE |
| interfaces (n̂ = mesh triangle normal) | plus side (n̂) | minus side |
|---|---|---|
| `body` (closed) | air | PEC |
| incidence | θᵢ = 90.0°, φᵢ = 180.0°, pol = [0.0, 0.0, 1.0] |
| observation | θ ∈ [0°, 180°] (181 samples), φ cuts = 0.0°, 90.0° |

## 3 · Performance & Efficiency

| triangles/tets | nodes | unknowns | t_mesh | t_assemble | t_solve | t_RCS | t_plots |
|---|---|---|---|---|---|---|---|
| 1236 | 620 | 1854 | 162.4 s | 1.9 s | 16.0 s | 0.6 s | 3.2 s |

| total time | assembly rate | LU throughput |
|---|---|---|
| 184.1 s | 0.00179 × 10⁹ interactions/s | 0.3 Gflop/s |

*Assembly rate counts N² impedance-matrix interactions; LU throughput uses the dense (2/3)·N³ flop model (single node, default BLAS threads).*
Machine-readable: `perf.csv`, `rcs.csv`.

## 4 · Results & Comparison

Bistatic RCS — dense MoM solution:

![RCS cuts](rcs_cuts.png)

Normalized far-field |E| pattern:

![Far-field polar](farfield_polar.png)

Surface-current magnitude |J| (dB, normalized to peak):

![Current distribution](current_views.png)


Cross-method validation — MLFMA distributed GMRES (MPI distributed GMRES (P = 2 ranks); leaf = λ/2, restart = 200, tol = 10⁻⁶) vs dense MoM (LU) reference, LU total 0.6 s:

| phi cut | RMSE MoM vs MLFMA [dB] | verdict |
|---|---|---|
| 0.0° | 0.016 | pass |
| 90.0° | 0.007 | pass |

## 5 · Key Conclusions

- Electric resolution: mesh size 0.06 m = 0.24 λ at f = 1200.0 MHz (90.0°, 180.0° plane-wave incidence).
- Accuracy: no analytic reference exists for this geometry; dense MoM (LU) and fast MLFMA (GMRES) solutions agree to 0.016 dB worst-cut RMSE → **PASS** against the 1 dB cross-method acceptance line.
- Throughput: impedance assembly 0.00179 G-interactions/s, dense LU 0.3 Gflop/s at N = 1854 unknowns.
- Dominant cost: gmsh meshing — 162.4 s (88.0% of the 184.1 s total).
- Bistatic RCS dynamic range over the observed cuts: -11.3 … -5.1 dBsm.

## Artifacts

- geometry: `geometry_views.png` · RCS data: `rcs.csv` · performance: `perf.csv`
- plots: `rcs_cuts.png`, `farfield_polar.png`
