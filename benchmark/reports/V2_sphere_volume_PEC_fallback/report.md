# EMMoMSuite Validation Report — `V2_sphere_volume_PEC_fallback`

| | |
|---|---|
| solver | EMMoMSuite v0.3.1 (Julia 1.12.3) |
| generated | 2026-10-10 19:56:15 |
| pipeline | geometry → gmsh mesh → MoM solve → RCS / far-field → this report |

## 1 · Geometry & Mesh

![geometry views](geometry_views.png)

![mesh views](mesh_views.png)

Boundary discretization: **5199 boundary triangles (of 5199 tetrahedra)**, 1176 nodes; geometry source `cases/geo/sphere_r0p5.geo` (gmsh/OpenCASCADE).

## 2 · Simulation Configuration

| item | value |
|---|---|
| geometry | `cases/geo/sphere_r0p5.geo` |
| mesh size | 0.08 m |
| mesh dim | 3 (tetrahedral volume mesh) |
| frequency | 1200.0 MHz (λ = 0.25 m) |
| formulation | EFIE (all-PEC fallback, surface extracted) |
| regions (Physical Volume) | material | tag |
|---|---|---|
| `body` | PEC | 2 |
| incidence | θᵢ = 90.0°, φᵢ = 180.0°, pol = [0.0, 0.0, 1.0] |
| observation | θ ∈ [0°, 180°] (181 samples), φ cuts = 0.0°, 90.0° |
| reference | analytic Mie series, PEC sphere r = 0.5 m |

## 3 · Performance & Efficiency

| triangles/tets | nodes | unknowns | t_mesh | t_assemble | t_solve | t_RCS | t_plots |
|---|---|---|---|---|---|---|---|
| 5199 | 1176 | 1875 | 218.7 s | 3.8 s | 1.5 s | 1.9 s | 7.6 s |

| total time | assembly rate | LU throughput |
|---|---|---|
| 233.5 s | 0.000933 × 10⁹ interactions/s | 2.8 Gflop/s |

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
| 0.0° | 0.082 | pass |
| 90.0° | 0.058 | pass |

Cross-method validation — dense MoM (LU) vs MLFMA (GMRES, leaf = λ/2, block-Jacobi preconditioned, restart = 200, tol = 10⁻⁶), total 290.8 s:

| phi cut | RMSE MoM vs MLFMA [dB] | verdict |
|---|---|---|
| 0.0° | 0.010 | pass |
| 90.0° | 0.006 | pass |

## 5 · Key Conclusions

- Electric resolution: mesh size 0.08 m = 0.32 λ at f = 1200.0 MHz (90.0°, 180.0° plane-wave incidence).
- Accuracy: worst-cut RMSE vs analytic Mie reference = 0.082 dB over 2 cuts → **PASS** against the 0.5 dB acceptance line.
- Throughput: impedance assembly 0.000933 G-interactions/s, dense LU 2.8 Gflop/s at N = 1875 unknowns.
- Dominant cost: gmsh meshing — 218.7 s (94.0% of the 233.5 s total).
- Bistatic RCS dynamic range over the observed cuts: -2.2 … -0.1 dBsm.

## Artifacts

- geometry: `geometry_views.png` · RCS data: `rcs.csv` · performance: `perf.csv`
- plots: `rcs_cuts.png`, `farfield_polar.png`
