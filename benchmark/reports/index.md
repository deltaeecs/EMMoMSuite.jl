# EMMoMSuite gmsh-driven RCS case family — index

![](../../docs/assets/logo.svg)

*Generated: 2026-10-10 17:32:53*

Pipeline: geometry file → gmsh mesh → MoM solve → RCS → far-field plots → publication report.

| case | geometry | IE | elements | unknowns | RMSE vs ref [dB] | MoM–MLFMA [dB] | total [s] | LU [Gflop/s] | report |
|---|---|---|---|---|---|---|---|---|---|
| [`G10_sphere_PMCHW_Mie`](G10_sphere_PMCHW_Mie/report.md) | sphere_r0p15.geo | PMCHW | 1784 | 2676 | 0.013 | 0.000 | 360.4 | 0.1 | [report.md](G10_sphere_PMCHW_Mie/report.md) |
| [`G12_dielectric_cube`](G12_dielectric_cube/report.md) | cube_r0p1.geo | VEFIE | 1602 | 3558 | — | 0.000 | 376.3 | 58.4 | [report.md](G12_dielectric_cube/report.md) |
| [`G1_sphere_CFIE_Mie`](G1_sphere_CFIE_Mie/report.md) | sphere_r0p5.geo | CFIE | 790 | 1185 | 0.137 | 0.000 | 287.4 | 0.2 | [report.md](G1_sphere_CFIE_Mie/report.md) |
| [`G5_plate_EFIE`](G5_plate_EFIE/report.md) | plate.geo | EFIE | 470 | 705 | — | 0.000 | 261.9 | 0.1 | [report.md](G5_plate_EFIE/report.md) |
