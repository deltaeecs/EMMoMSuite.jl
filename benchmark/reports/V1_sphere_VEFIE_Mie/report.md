# EMMoMSuite Validation Report — `V1_sphere_VEFIE_Mie`

| | |
|---|---|
| solver | EMMoMSuite v0.3.1 (Julia 1.12.3) |
| generated | 2026-10-10 19:30:57 |
| pipeline | geometry → gmsh mesh → MoM solve → RCS / far-field → this report |
| parallel | MPI distributed GMRES (P = 2 ranks, SAI precond, leaf = 0.1 m, rel-res = 1.5e-14) |

## 1 · Geometry & Mesh

![geometry views](geometry_views.png)

![mesh views](mesh_views.png)

Boundary discretization: **1175 boundary triangles (of 1175 tetrahedra)**, 324 nodes; geometry source `cases/geo/sphere_r0p15.geo` (gmsh/OpenCASCADE).

## 2 · Simulation Configuration

| item | value |
|---|---|
| geometry | `cases/geo/sphere_r0p15.geo` |
| mesh size | 0.04 m |
| mesh dim | 3 (tetrahedral volume mesh) |
| frequency | 600.0 MHz (λ = 0.5 m) |
| formulation | VEFIE (SWG volume discretization) |
| regions (Physical Volume) | material | tag |
|---|---|---|
| `body` | Dielectric(εᵣ = 4.0, μᵣ = 1.0) | 2 |
| incidence | θᵢ = 90.0°, φᵢ = 180.0°, pol = [0.0, 0.0, 1.0] |
| observation | θ ∈ [0°, 180°] (181 samples), φ cuts = 0.0°, 90.0° |
| reference | analytic Mie series, dielectric sphere r = 0.15 m |

## 3 · Performance & Efficiency

| triangles/tets | nodes | unknowns | t_mesh | t_assemble | t_solve | t_RCS | t_plots |
|---|---|---|---|---|---|---|---|
| 1175 | 324 | 2578 | 162.4 s | 22.6 s | 0.4 s | 1.1 s | 1.0 s |

| total time | assembly rate | LU throughput |
|---|---|---|
| 187.5 s | 0.000294 × 10⁹ interactions/s | 27.4 Gflop/s |

*Assembly rate counts N² impedance-matrix interactions; LU throughput uses the dense (2/3)·N³ flop model (single node, default BLAS threads).*
Machine-readable: `perf.csv`, `rcs.csv`.

## 4 · Results & Comparison

Bistatic RCS — MoM vs analytic reference:

![RCS cuts](rcs_cuts.png)

Normalized far-field |E| pattern:

![Far-field polar](farfield_polar.png)

Surface-current magnitude |J| (dB, normalized to peak):

![Current distribution](current_views.png)


Cross-method validation — MLFMA distributed GMRES (MPI distributed GMRES (P = 2 ranks, SAI precond, leaf = 0.1 m, rel-res = 1.5e-14); leaf = λ/2, restart = 200, tol = 10⁻⁶) vs dense MoM (LU) reference, LU total 11.2 s:

| phi cut | RMSE MoM vs MLFMA [dB] | verdict |
|---|---|---|
| 0.0° | 0.000 | pass |
| 90.0° | 0.000 | pass |

## 5 · Key Conclusions

- Electric resolution: mesh size 0.04 m = 0.08 λ at f = 600.0 MHz (90.0°, 180.0° plane-wave incidence).
- Accuracy: no analytic reference exists for this geometry; dense MoM (LU) and fast MLFMA (GMRES) solutions agree to 0.0 dB worst-cut RMSE → **PASS** against the 1 dB cross-method acceptance line.
- Throughput: impedance assembly 0.000294 G-interactions/s, dense LU 27.4 Gflop/s at N = 2578 unknowns.
- Dominant cost: gmsh meshing — 162.4 s (87.0% of the 187.5 s total).
- Bistatic RCS dynamic range over the observed cuts: -11.9 … -6.4 dBsm.

## Artifacts

- geometry: `geometry_views.png` · RCS data: `rcs.csv` · performance: `perf.csv`
- plots: `rcs_cuts.png`, `farfield_polar.png`
