# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added
- **Code Style**: Applied `JuliaFormatter.jl` (v1.0.62) to all 79 `src/` files (default style, `indent=4`, `margin=100`). Configuration in `.JuliaFormatter.toml`.
- **Scripts**: Added `scripts/format_code.jl` (format `src/`) and `scripts/install_formatter.jl` (global install helper).
- **Basis Functions**: `PWC` (Piecewise Constant, 3 DOF/tet) and `RBF` (Rao-Billinghast-Farris) basis functions.
- **Geometry**: Hexahedra mesh support (`HexahedraInfo`, `Quads4Hexa`); Nastran CHEXA element reading with continuation lines.
- **Integral Equations**: VEFIE (`EFIE+SWG`, `VEFIE+PWC`, `VEFIE+RBF`, `VEFIE+PWCHex`), SCFIE mixed-media formulations (18 `assemble_impedance_matrix` methods, 17 `excitation_vector` methods).
- **Solvers**: `BlockJacobiPreconditioner` using per-cube diagonal block LU factorization (166× faster build vs sparse LU; suitable for CFIE).
- **Post-Processing**: `RadiationIntegral` and `radarCrossSection` for all basis function types (5 RCS methods).

### Performance (Phase 8.9 – 8.10) — All cases **beat Legacy** on 4 threads

| Case | Legacy | EMSuite | Speedup |
|------|--------|---------|---------|
| Jet EFIE direct fill (`N=26424`) | 20.70 s | **4.26 s** | **4.9×** |
| Jet CFIE direct fill (`N=26424`) | 168.29 s | **14.48 s** | **11.6×** |
| VEFIE SWG (`N_tet=7278`) | 46.13 s | **41.30 s** | **1.12×** |
| SCFIE total (`N=15860`) | 66.68 s | **65.67 s** | **1.02×** |

Key algorithmic improvements:
- **EFIE/MFIE/CFIE**: SIMD `@fastmath` inner loop, `@inbounds` cache-friendly access, single-pass CFIE (no extra allocation).
- **VEFIE**: Upper tet-triangle symmetry (`Z_st[j,i] = (κ_t/κ_s)·Z_ts[i,j]`); reduces kernel calls from N² → N(N+1)/2.
- **SCFIE coupling**: Reciprocity `Z_VS = Z_SV^T / κ`; cyclic row-parallel scheduling; Fss boundary parallel.

### Changed
- Refactored legacy code (`MoM_Basics`, `MoM_Kernels`, `MoM_MPI`) into a unified `EMSuite` package.
- `Driver.jl` supports multi-IE type (`EFIE`/`MFIE`/`CFIE`/`VEFIE`/`SCFIE`).
- `Configuration.jl`: added `ie_type`, `cfie_alpha`, `permittivities` fields.
- MPI parallel CFIE assembly uses independent EFIE+MFIE interactions + linear combination.
- MLFMA `Z_near` initial COO allocation capped at 10M entries/thread (prevents OOM for large problems).
- MLFMA `octree` leaf box size default aligned with legacy (0.25λ).
- Improved MLFMA accuracy factors calibrated to match legacy `MoM_Lebedev`.
- Standardized I/O and logging.

### Fixed
- Fixed singular term scaling in EFIE (1/R and ρ·ρ'/R terms).
- Fixed edge sorting logic in RWG basis functions.
- Fixed VEFIE revert: column-parallel row-scatter caused cache thrash (+11% slowdown); reverted to row-parallel before symmetry optimization.
- Fixed Julia 1.12 `threadid()` compatibility.
- Verified parallel consistency with serial execution (179/179 tests Green).

---

## [0.1.0] — Initial Release (Legacy Parity)

### Added
- **Core**: Initial project structure and core interfaces (`AbstractMesh`, `AbstractBasisFunction`, etc.).
- **Geometry**: Support for Triangle, Tetrahedra meshes. Mesh I/O for Nastran (.nas) and Gmsh (.msh).
- **Basis Functions**: RWG and SWG basis functions.
- **Integral Equations**: EFIE, MFIE, and CFIE formulations with dense direct fill.
- **Solvers**: Direct (LU) and Iterative (GMRES, BiCGSTAB) solvers.
- **Fast Algorithms**: Multilevel Fast Multipole Algorithm (MLFMA) with Lebedev integration.
- **Parallel**: MPI support for distributed matrix assembly and solving.
- **Post-Processing**: RCS, Near Field, and Far Field calculations. VTK export.
- **Documentation**: Project structure and API documentation using Documenter.jl.
