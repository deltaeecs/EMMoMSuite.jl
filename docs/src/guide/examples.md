# Examples

This section outlines standard benchmark problems to validate the solver.

## 1. Scattering from a PEC Sphere (Mie Series Validation)
- **Geometry**: Sphere of radius $1\lambda$.
- **Formulation**: CFIE ($\alpha=0.5$) to avoid internal resonance.
- **Validation**: Compare RCS with analytical Mie series solution.
- **Key Feature**: Demonstrates accuracy of curved surface modeling.

## 2. NASA Almond (RCS Benchmark)
- **Geometry**: "Almond" shape defined by NASA for benchmarking.
- **Formulation**: EFIE.
- **Challenge**: High dynamic range of RCS; requires precise handling of traveling waves and tip diffraction.

## 3. Parallel Plate Capacitor (Low Frequency)
- **Geometry**: Two parallel circular plates.
- **Formulation**: EFIE with Loop-Star or specialized low-frequency basis functions (if available) or standard RWG with preconditioning.
- **Challenge**: Low-frequency breakdown (conditioning issues).

## 4. Array Antenna Radiation
- **Geometry**: 4x4 Dipole array backed by a ground plane.
- **Excitation**: Voltage gap sources at dipole centers.
- **Output**: Far-field radiation pattern and input impedance.

