# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added
- 待定：后续版本功能预留。

## [0.1.0] - 2026-08-07

### Changed
- **Rename**: Renamed the primary package identity from `EMSuite.jl` to `EMMoMSuite.jl` to make the Method of Moments focus explicit and distinguish the refactored package from legacy naming.
- **Consolidation**: 原分散 MoM 包（`MoM_Basics` / `MoM_Kernels` / `MoM_AllinOne` / `MoM_MPI` / `MoM_Lebedev` / `MoM_Visualizing` / `MPIArray4MoMs`）功能已完整并入本包；旧仓库已归档，网格与 FEKO 基线数据收编至 `deps/fixtures/`。
- **Dependencies**: 移除对未注册旧包（`MoM_AllinOne` / `MoM_Basics` / `MoM_Kernels`）的运行时依赖，包可独立安装。
- **Metadata**: 更新 `uuid`、`authors`，包名统一为 `EMMoMSuite`。

### Fixed
- **MLFMA near-field race**: `assemble_near_field` 使用动态 `Threads.@threads` 时，迭代任务迁移线程可能导致同一 per-thread 缓冲区被并发写入，`Z_near` 间歇性损坏；改为 `:static` 调度后结果完全确定。
- **SCFIE kernel constants**: `scfie_coupling_interaction` / `scfie_sv_only_interaction` 中未定义符号 `mu0`/`eps0` 改为 `Constants.mu0`/`Constants.eps0`。
- **Tests**: 修复 `test_coverage_gaps.jl` 非法 `@test expr "msg"` 语法与 API 引用；重写 `test_integral_equations_endtoend.jl` / `test_integral_equations_frequency.jl`（原按不存在 API 编写，从未跑通）；修正复数对称矩阵的转置口径。

### Added
- **Fixtures**: 新增 `deps/fixtures/`（旧包网格 + FEKO 基线 CSV），测试与 benchmark 自包含。

### Changed
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
- Refactored legacy code (`MoM_Basics`, `MoM_Kernels`, `MoM_MPI`) into a unified `EMMoMSuite` package.
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
