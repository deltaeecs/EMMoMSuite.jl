# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [0.2.1] - 2026-08-27

### Added

- **覆盖率采集套件**（`test/runtests_cov.jl`）：单线程常规运行 2m10s（≤3 分钟门），
  覆盖 Core / Geometry / BasisFunctions / IntegralEquations / FastAlgorithms /
  Parallel(MPI) / IO / Materials / Utilities / Solvers / Ports / PostProcessing /
  Accuracy；最新单次采集可执行行覆盖率 **94.2%**（9726/10323）。
- **MPI 体积装配覆盖率**：VEFIE/SCFIE × SWG/PWC/PWCHex 并行装配在 P=1 下与串行
  装配逐项对照（`test_cov_volume_assembly.jl`）。
- **CI 覆盖率硬门**：coverage job 切换到 `runtests_cov.jl`，并以
  `scripts/check_coverage.jl 80` 强制 ≥80% 行覆盖，Codecov 与本地口径一致。

### Fixed

- **MLFMA 插值（issue #20）**：`interpolationCSCMatCal` 构建两层插值矩阵时，当父层方向采样数
  小于默认插值阶（`nInterp = 6`）会按方向把阶数钳制为
  `nlocalInterpTheta`/`nlocalInterpPhi`（θ 为 `L+1`、φ 为 `2(L+1)`），但函数仍用**未钳制的
  `nlocalInterp`** 去定尺寸/索引稀疏数组，导致精细八叉树（如 `leafCubeEdgel = 1e-3`）下
  `MLFMAOperator` 崩溃（θ 步 `BoundsError`；更细叶子下 φ 步稀疏长度不匹配）。修复为统一使用
  钳制值（`rawIDϕs`→`nlocalInterpPhi`、`rawIDθs`→`nlocalInterpTheta`，以及 θ 邻域循环的
  `inGlobalIDs`/`for jInter` 上界），并新增回归测试覆盖 θ/φ 两条钳制路径。
- 记录 Julia 1.12 复用未插桩编译缓存导致 `.cov` 不生成的问题及可靠做法
  （先移开 `~/.julia/compiled/v1.12/EMMoMSuite` 再采集）。

## [0.2.0] - 2026-08-13

### Added

- **MPI 混合并行算子**：`MLFMAOperatorMPI` / `AIMOperatorMPI` / `PMCHWMLFMAOperatorMPI`
  （近场按 cube/行分区、远场全层按秩分区 + 每层 Allreduce，秩内 `@threads`）；
  FFTSpectral / Lebedev / PMCHW 分支 MPI 适配；层间插值矩阵线程化。
- **分布式稠密直接求解**：ScaLAPACK（本机 MinGW/MSMPI）分布式稠密 LU
  （BLACS + `pzgesv`，库路径自动探测 + 安装指引）。
- **MPI 分布式预条件**：`DistributedBlockJacobiPreconditioner`（块按秩 1/P 分布）与
  `DistributedDiagonalPreconditioner`，接入 `distributed_gmres!` 左预条件 `Pl`。
- **PMCHW 近场原生装配**：按 octree 近邻对直接计算（不再装配全稠密 2N×2N 再提取），
  支持 MPI cube 分区。

### Fixed

- **既有 RED 门**：GD2A/GD2R（共享 EFIE 聚合/解聚 3→4 点求积统一，parity 机器精度）、
  B2（GMRES 全空间 restart，Zin 误差 0.2%）。
- **既有 broken 门**（GD2S/GD2L/GD2U/GD2V）：根因定位为对角实谱 M2L 的距离约束
  （叶层 far 偏移 <8 失效），门改用有效配置并显式验证阈值；`build_octree` 对越界
  配置给出 `@warn` 防护。
- **ScaLAPACK 可移植性**：库自动探测（环境变量 / MSYS2 路径 / PATH），缺失时清晰报错。
- **AIM 近场校正**与算子 k/eta/factor 一致；预条件块类型稳定化。

### Performance

- PMCHW matvec：N=594 MPI 3.37s → 0.107s（31×）。
- 远场翻译因子 αTrans：每秩 ~15.7GB → ~21MB（740×）；PMCHW 构造 137s → 10s。
- ScaLAPACK 分布式 LU：N≥4800 单机反超多线程 OpenBLAS。
- 端到端每秩峰值内存：17.6GB → 1.3–1.6GB（N=594 PMCHW）。

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

### Fixed

- **MLFMA 近场半径自适应（issue #22 问题 3）**：`build_octree` / `MLFMAOperator` /
  `MLFMAOperatorMPI` 的 `near_range` 默认值改为 `nothing`（由叶层推导全树统一的
  自适应半径，按 `adaptive_near_range` 保证每层最小远对 `kR_min ≥ 0.55·L`，实谱 M2L
  浮点收敛；经 GD2V 门校准），显式 `Int` 仍被尊重；叶层 kR 距离校验阈值由绝对「偏移 ≥ 8」
  改为相对「kR_min ≥ 0.55·L」，避免静默精度损失。自适应半径 > 2（电小叶层）时输出
  效率指引警告，并新增权衡基准脚本 `benchmark/mlfma_nearrange_tradeoff.jl`。
  问题 3 的性能型根治（shifted expansion 恢复标准 189 交互列表）另行立项。
- **RCS 静默 -Inf（issue #22 问题 2）**：`radarCrossSection` 全局 `k0` 未设置
  （默认 0）时不再静默返回全 `-Inf` dBsm，改为抛出可行动错误。

### Changed

- **MLFMA matvec 分配（issue #22 问题 1）**：`MLFMAOperator` 构造时预计算
  单元几何与求积数据（`element_cache`），`aggregate_leaf!` / `disaggregate_leaf!`
  不再每次 `mul!` 重建；`aggregate_upward!` / `disaggregate_downward!` 的
  per-(盒,kid) 临时数组改为按线程 scratch 复用。
- **MLFMA 预条件用法（issue #22 问题 4）**：`MLFMAOperator` 文档示例改用
  `BlockJacobiPreconditioner`（按盒分块并行 LU）替代 `lu(Z_near)`（O(N³)）。
- **统一网格对比（issue #22 问题 5）**：README（中/英）新增「统一网格：跨路径
  数值对比」指引（Gmsh 路径生成一份网格供稠密/MLFMA/PMCHW 共用），并新增
  同网格稠密 vs MLFMA parity 回归测试。

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
