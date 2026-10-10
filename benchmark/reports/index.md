# EMMoMSuite gmsh-driven RCS case family — index

![](../../docs/assets/logo.svg)

*Generated: 2026-10-10 19:57:33*

Pipeline: geometry file → gmsh mesh → MoM solve → RCS → far-field plots → publication report.

| case | geometry | IE | elements | unknowns | RMSE vs ref [dB] | MoM–MLFMA [dB] | total [s] | LU [Gflop/s] | report |
|---|---|---|---|---|---|---|---|---|---|
| [`G10_sphere_PMCHW_Mie`](G10_sphere_PMCHW_Mie/report.md) | sphere_r0p15.geo | PMCHW | 3196 | 4794 | 0.214 | 0.074 | 1164.6 | 0.1 | [report.md](G10_sphere_PMCHW_Mie/report.md) |
| [`G11_lossy_sphere_PMCHW_Mie`](G11_lossy_sphere_PMCHW_Mie/report.md) | sphere_r0p15.geo | PMCHW | 3196 | 4794 | 0.063 | 0.026 | 421.3 | 0.8 | [report.md](G11_lossy_sphere_PMCHW_Mie/report.md) |
| [`G12_dielectric_cube`](G12_dielectric_cube/report.md) | cube_r0p1.geo | VEFIE | 1602 | 3558 | — | 0.000 | 228.0 | 23.7 | [report.md](G12_dielectric_cube/report.md) |
| [`G12b_dielectric_cube_fine`](G12b_dielectric_cube_fine/report.md) | cube_r0p1.geo | VEFIE | 2725 | 5939 | — | 0.000 | 289.9 | 306.9 | [report.md](G12b_dielectric_cube_fine/report.md) |
| [`G1_sphere_CFIE_Mie`](G1_sphere_CFIE_Mie/report.md) | sphere_r0p5.geo | CFIE | 3188 | 4782 | 0.222 | 0.002 | 236.6 | 1.6 | [report.md](G1_sphere_CFIE_Mie/report.md) |
| [`G2_sphere_EFIE_Mie`](G2_sphere_EFIE_Mie/report.md) | sphere_r0p5.geo | EFIE | 3188 | 4782 | 0.018 | 0.007 | 223.8 | 1.2 | [report.md](G2_sphere_EFIE_Mie/report.md) |
| [`G3_box_EFIE`](G3_box_EFIE/report.md) | box_1m.geo | EFIE | 4136 | 6204 | — | 0.060 | 225.0 | 2.8 | [report.md](G3_box_EFIE/report.md) |
| [`G4_cylinder_EFIE`](G4_cylinder_EFIE/report.md) | cylinder.geo | EFIE | 1388 | 2082 | — | 0.033 | 173.9 | 0.4 | [report.md](G4_cylinder_EFIE/report.md) |
| [`G5_plate_EFIE`](G5_plate_EFIE/report.md) | plate.geo | EFIE | 1514 | 2271 | — | 0.357 | 177.1 | 0.5 | [report.md](G5_plate_EFIE/report.md) |
| [`G6_ellipsoid_EFIE`](G6_ellipsoid_EFIE/report.md) | ellipsoid.geo | EFIE | 1236 | 1854 | — | 0.007 | 184.1 | 0.3 | [report.md](G6_ellipsoid_EFIE/report.md) |
| [`G7_cone_EFIE`](G7_cone_EFIE/report.md) | cone.geo | EFIE | 1498 | 2247 | — | 0.029 | 196.9 | 0.3 | [report.md](G7_cone_EFIE/report.md) |
| [`G8_torus_EFIE`](G8_torus_EFIE/report.md) | torus.geo | EFIE | 1524 | 2286 | — | 0.013 | 173.6 | 1.0 | [report.md](G8_torus_EFIE/report.md) |
| [`G9_box_CFIE`](G9_box_CFIE/report.md) | box_1m.geo | CFIE | 4136 | 6204 | — | 0.041 | 380.9 | 1.1 | [report.md](G9_box_CFIE/report.md) |
| [`V1_sphere_VEFIE_Mie`](V1_sphere_VEFIE_Mie/report.md) | sphere_r0p15.geo | VEFIE | 1175 | 2578 | — | 0.000 | 187.5 | 27.4 | [report.md](V1_sphere_VEFIE_Mie/report.md) |
| [`V2_sphere_volume_PEC_fallback`](V2_sphere_volume_PEC_fallback/report.md) | sphere_r0p5.geo | EFIE | 5199 | 1875 | 0.058 | 0.006 | 233.5 | 2.8 | [report.md](V2_sphere_volume_PEC_fallback/report.md) |
