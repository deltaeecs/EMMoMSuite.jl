# EMMoMSuite gmsh-driven RCS case family — index

![](../../docs/assets/logo.svg)

*Generated: 2026-10-09 23:22:49*

Pipeline: geometry file → gmsh mesh → MoM solve → RCS → far-field plots → publication report.

| case | geometry | IE | elements | unknowns | Mie RMSE [dB] | total [s] | LU [Gflop/s] | report |
|---|---|---|---|---|---|---|---|---|
| [`G10_sphere_PMCHW_Mie`](G10_sphere_PMCHW_Mie/report.md) | sphere_r0p15.geo | PMCHW | 1784 | 2676 | 0.044 / 0.013 | 229.3 | 3.8 | [report.md](G10_sphere_PMCHW_Mie/report.md) |
| [`G12_dielectric_cube`](G12_dielectric_cube/report.md) | cube_r0p1.geo | VEFIE | 1602 | 3558 | — | 256.2 | 24.0 | [report.md](G12_dielectric_cube/report.md) |
| [`G1_sphere_CFIE_Mie`](G1_sphere_CFIE_Mie/report.md) | sphere_r0p5.geo | CFIE | 380 | 570 | 0.353 / 0.297 | 217.0 | 0.1 | [report.md](G1_sphere_CFIE_Mie/report.md) |
| [`G5_plate_EFIE`](G5_plate_EFIE/report.md) | plate.geo | EFIE | 470 | 705 | — | 200.1 | 0.6 | [report.md](G5_plate_EFIE/report.md) |
