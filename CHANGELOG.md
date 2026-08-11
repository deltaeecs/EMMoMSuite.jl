# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added

- **FastAlgorithms.ACA / MLACA**：部分主元 ACA（转置约定，无共轭）与 QR/SVD
  再压缩；`ACAOperator`/`MLACAOperator`（复用八叉树聚类与近场稀疏装配，GMRES
  集成）；非对称支持（CFIE，`symmetric=false` 双向压缩）；PMCHW 2N 多基函数
  支持（`PMCHWBlockEvaluator`，EJ/EM/HJ/HM 子块）；叶层分块直接 LU
  （`block_lu`/`block_lu_solve`，多 RHS，默认 `recompress=false`）。
- **MLFMA 参数配置化**：`nInterp`、`precision_digits` 构造参数（串行与 MPI），
  ILU/SPAI/BlockJacobi 便捷接线。
- **大规模本地基准**：`benchmark/run_large_fast_solvers_benchmark.jl`
  （N 至 1.1 万、多配方、稠密参照、finite/NaN 检查；本地手动运行，不进 CI）。

### Fixed

- `aca` 对整数类型矩阵的 `eps(Int)` 崩溃（改用 `eps(float(real(T)))`）。
- `block_lu` 对 PMCHW（2N 系统）的错误分区（叶层块展开为 J/M 双通道）。
- `block_lu` 再压缩在病态系统（PMCHW cond≈1e6）下因子误差复合放大，
  默认改为 `recompress=false`，再压缩仅建议良态系统显式开启。
- `get_leaf_intervals` 对 ACA/MLACA 算子缺失（原实现对 octree 调用不存在的方法）。

## [0.1.1] - 2026-08-09

### Fixed

- **Lebedev 插值权重矩阵**：定位并修复权重矩阵全部失效的根本原因——数据集生成时粗/细层
  各自重新随机生成源几何（两层数据对应不同辐射函数，插值在拟合噪声）；改为每样本生成
  一次几何并供两层共用。同时修复 `vcat(real,imag)` 只使用实部约束（改 `hcat` 全复约束）、
  源几何与层盒子失配、LVI `truncation_kernel` 参数（`ka`→`a/λ`）与 `2L+1<max` off-by-one、
  高阶回退调用不存在方法等问题。
- **论文对齐**：修复后论文 Fig.2 插值误差可复现（8 个插值点 εi≈2.5e-4~1e-3，随点数递减）。

### Added

- **FastAlgorithms.Lebedev.SHInterp**：球谐精确、局部约束、八面体群轨道压缩、笛卡尔标量 SH
  矢量插值（机器精度）、自旋加权球谐（VSH）与"数据拟合 + 精确性约束"混合权重
  （`interp_weights_*` 系列）。`LbTrainedInterp1tepInfo` 默认 `method=:sh_auto`。
- **高阶节点**：`fibonacci_grid` / `high_order_nodes`，p>131（无 Lebedev 数据集）时自动使用
  Fibonacci 准均匀格点；`LbPolesInfo` 显式携带多项式阶数，插值链路支持任意点数。
- **MLFMA 一步插值路径**：`build_octree(; interp_method=Val(:LbTrained1Step))`、
  Aggregation/Disaggregation 一步插值分支、`MLFMAOperator` 透传 `interp_method`/`near_range`。

### Performance

- 4λ PEC 球基准（N=3312，5 层）：Lebedev 路径极点数较 GL 少 31.5%、单次 MVM 快约 18~2.9×
  （随规模增长），εq 与 GL 完全一致。

## [Unreleased]

### Changed

- **License**: 项目许可证从 MIT 切换为 GPL-3.0-only，与原 MoM 系列包（`MoM_Basics` / `MoM_Kernels` / `MoM_AllinOne` / `MoM_MPI` / `MoM_Lebedev` / `MoM_Visualizing`）保持一致；MIT LICENSE 已从 git 历史中移除。
- **README**: 新增 CI 与 `master` 分支覆盖率 badge，新增英文版 [README.en.md](README.en.md)；许可证说明更新为“学术研究免费使用、不推荐商业环境”。
- **CI**: fast tests 开启 `--code-coverage=user` 并上传真实覆盖率报告，修复 Codecov 报告显示 0% 的问题。
- **CI**: coverage 任务改用 `test/runtests_light_cov.jl`（跳过 MLFMA 与大型装配，覆盖率插桩下过慢），避免全量套件插桩运行耗尽 runner 内存导致任务被终止。
- **Repo hygiene**: 从 git 历史移除非必要二进制与生成文件——`deps/InterpolationWeights/`（约 265 MB HDF5 插值权重缓存，无调用者且可按需重新生成）、历史遗留的 `sphere_mesh_data.jld2`，以及生成的测试/运行产物（`results/`、`test_results/`、`test/test_results/`、`TestPlate.*`）。

### Added

- 英文版 README（`README.en.md`）。

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
