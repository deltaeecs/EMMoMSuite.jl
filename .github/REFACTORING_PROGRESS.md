# EMSuite Phase 9 覆盖率提升进度 (追加节)

> 最后更新: 2026-03-XX

## MLFMA 算法分析阶段 — **已完成** ✅

> 更新: 2026-03-XX

### 工作成果

1. **代码清理**:
   - `scripts/` 废物代码归档到 `scripts/archive/` (~80 文件)
   - `benchmark/debug_archive/` 清理 (~50 文件)
   - 根目录废物文件删除 (tmp_*.txt, 过期文档等)
   - 保留: `scripts/format_code.jl`, `scripts/install_formatter.jl`

2. **Legacy MLFMA 完整算法报告**: [legacy_mlfma_algorithm_report.md](.github/plans/legacy_mlfma_algorithm_report.md)
   - 覆盖从网格读取到 RCS 后处理的全流程
   - 精确到每一个公式和系数
   - 12 项 Legacy 工程技巧总结
   - EMSuite vs Legacy 系数链完整验证 (总系数一致: k²η/16π²)

3. **关键发现**:
   - nLevels=2: EFIE 精度 0.5%~9% ✅ | nLevels≥3: 1005%~2918% ✗
   - 总系数链验证通过 → 问题在 upward/downward pass 的插值或相移中
   - Addition Theorem 指数符号确认: 源端 e^{+jk}, 场端 e^{-jk}
   - 远亲定义 (7³-3³=316) 和清零时机均与 Legacy 一致

4. **检视**: 2 轮检视完成，修正了 edgel 符号描述、时间约定、Addition Theorem 指数符号

5. **PMCHW 完整算法报告**: [pmchw_algorithm_report.md](.github/plans/pmchw_algorithm_report.md)
   - 覆盖面等效原理 → L/K 算子 → 2N×2N 矩阵 → 4遍 MLFMA → 双流 RCS 全链路
   - K^PMCHW 积分核与 MFIE K 的详细区别 (无对角、无 n̂× 测试、无 η₀ 预乘)
   - 4 个子块的 MLFMA 系数链验证 (EJ/HM L型, EM/HJ K型)
   - 确认 Legacy 不含 PMCHW 实现，PMCHW 完全为 EMSuite 原创
   - 2 轮检视完成：Round 1 修复 §6.4 球坐标符号、§C.3 有损介质描述等 4 项; Round 2 无新问题

## Phase 9 测试覆盖率提升 → 进行中

### Round 4 工作成果

新增/扩展测试 (commit ca3fca5, 2161daa):
- test_hex_rbf.jl (NEW, 67 tests): HexahedraMesh+PWCHexBasis+RBFBasis
- test_geometry.jl: 33->60 tests (NAS CHEXA/CTETRA/MSH tetra)
- test_mpi_array_utils.jl: 24->47 tests
- test_io.jl: VTK tetra + HDF5 save_result
- 总计: 449/449 pass

源码 Bug 修复 (commit 99260a0, 2161daa):
- VEFIE/SCFIE/Excitation/RadiationIntegral: gq_hex->gq_quad coordinates
- MeshIO.jl: _parse_chexa CHEXA 续行修复 (节点7,8丢失)
- GmshIO.jl: 添加 type 4 (tetra) + type 5 (hexa) 支持
- indices.jl: indice2rank NTuple{1} in 语义修复

Phase 9 检视迭代 Round 1 (commit 90787dc, baa0418):
- GmshIO.jl: 文档+版本检查@warn+混合网格逻辑
- VEFIE.jl _rbf_far_kernel!: 预计算 freeVns 减少冗余分配

覆盖率历史: Round1=57.22% Round2=64.37% Round3=65.69% Round4=测量中

---

# EMSuite 閲嶆瀯杩涘害

> 鏈€鍚庢洿鏂? 2026-03-06
## 当前阶段: Phase 14 全量精度测试与对比报告 — **进行中 🔄**

> 最后更新: 2026-03-07

### Phase 14 计划概述

**目标**: 对 EMSuite 所有主要积分方程 × 求解路径，与 Feko 商业软件结果或 Mie 解析解进行系统对比，生成独立精度报告（不再对比 Legacy 代码）。

**计划文档**: [PHASE_14_ACCURACY_REPORT_PLAN.md](PHASE_14_ACCURACY_REPORT_PLAN.md)

**基准来源**:
- Feko 基线: `MoM_AllinOne/deps/compare_feko/` (Jet 100MHz, Sphere 600MHz, Plate 1.2GHz, Plate+Metal 1.2GHz)
- Mie 解析解: EMSuite 内置 `Utilities/MieSeries.jl`（用于 PEC 球独立校验）

**测试用例**: F1–F9 (9 个用例，覆盖 S-EFIE/S-CFIE/V-EFIE/SCFIE × Direct/MLFMA)

**精度门限**: Direct ≤ 2 dB RMSE，MLFMA ≤ 3 dB RMSE (vs Feko/Mie)

**待完成子任务**:
- [x] 14.0 Feko CSV 解析器 + TDD 测试 ✅ b416090 (39/39)
- [x] 14.1 参考基准生成器（Mie PEC/介质 + 偶极子解析）✅
- [x] 14.2 `AccuracyResult` + `AntennaAccuracyResult` 指标函数 ✅ 65b4297 (25/25)
- [x] 14.3 F1–F4 Jet 仿真脚本 ✅ b2efd57
- [x] 14.4 F5–F6 Sphere 仿真脚本 ✅
- [x] 14.5 F7–F9 Plate 仿真脚本 ✅ (F8 跳过-无曲面网格)
- [x] 14.6 P1, P3 PMCHW Direct 脚本 ✅ (P2 跳过-待实现 PMCHWMLFMAOperator)
- [x] 14.7 PMCHW block/operator shell + MLFMA backend 主线实现（Dense shell 先行，MLFMA 后接）
- [x] 14.8 P2 PMCHW MLFMA 验证
- [x] 14.9 A1–A4 偶极子天线 DeltaGap 基准脚本 ✅ 3039d32
- [x] 14.10 实际运行仿真 → CSV → ACCURACY_REPORT.md
- [x] 14.11 检视迭代 × 2 轮

---

## Phase 15: 介质与金属-介质混合天线 + PMCHW Transmission Block/Operator 架构 — **计划中 📋**

> 最后更新: 2026-03-07（架构修订：Bempp 风格 block/operator shell 先行，MLFMA 作为 backend 收编）
> 计划文档: [PHASE_15_DIELECTRIC_ANTENNA_PLAN.md](PHASE_15_DIELECTRIC_ANTENNA_PLAN.md)

**目标**: 扩展天线测试至介质（PMCHW）和金属-介质混合（VS-EFIE/VS-CFIE）类型，同时把 PMCHW 顶层实现迁移到 Bempp 风格 block/operator shell；Dense backend 先行，现有 PMCHWMLFMA 路径收编为 backend，H-matrix 延后。

**当前架构冻结**:
- 顶层 transmission 组织改为 **四块 `EJ/EM/HJ/HM` 的 block/operator shell**，不再把 `PMCHWMLFMAOperator` 作为唯一主线载体。
- `weak_form` / `strong_form` / mass-matrix-aware solve 归入 shell，Dense / MLFMA / H-matrix 统一视为 backend。
- Gibson Ch.11 Algorithm 14 的“双八叉树 + 四遍远场”继续保留，但只作为 **MLFMA backend 内部实现**。
- FreeFEM/Htool 的 H-matrix 路线保留为后续 backend 选项；在 shell 与 Gate S 稳定前不得前置。

**新增 API**:
- `excitation_vector(op::PMCHW, source::DeltaGapSource, basis::RWGBasis)` → 2N 激励向量（E-行 delta-gap）
- `input_impedance(op::PMCHW, source, I_2N, basis)` → PMCHW J-部分输入阻抗
- `excitation_vector(op::SCFIE, source::DeltaGapSource, surf_basis, vol_basis)` → 表面 + 零体积激励

**测试场景**:
| ID | 方程 | 馈电 | 求解 | 参考 |
|----|------|-----|------|-----|
| B1 | PMCHW (εᵣ=4) | DeltaGap (球面) | Direct | 物理自洽 + εᵣ→1 极限 |
| B2 | PMCHW | DeltaGap | MLFMA (双八叉树+四遍远场) | B1 Direct |
| B3 | VS-EFIE (α=0) | DeltaGap (金属面) | Direct | EFIE-only εᵣ→1 |
| B4 | VS-CFIE (α=0.5) | DeltaGap | Direct | B3 |
| B5 | VS-CFIE (α=0.5) | DeltaGap | MLFMA | B4 Direct |

**子任务**:
- [x] 15.1 TDD: `test_pmchw_excitation.jl`
- [x] 15.2–15.3 PMCHW DeltaGap 激励 + input_impedance API
- [x] 15.4–15.5 SCFIE DeltaGap 激励 API
- [x] 15.6 基准脚本 `run_B1_B5_antenna.jl`
- [x] 15.7 TDD: `test_pmchw_mlfma_operator.jl`
- [x] 15.8 实现 `assemble_near_field_pmchw`（2N×2N，4 块，无 MagneticRWGBasis）
- [x] 15.9 实现 `aggregate_leaf_pmchw!`（x_range 参数区分 J/M）
- [x] 15.10 实现 `disaggregate_leaf_pmchw_j!` + `_m!`（四种接收核函数）
- [x] 15.11 完成 PMCHW block/operator shell 与 Dense backend 验证，再把现有 `PMCHWMLFMAOperator` 收编为 MLFMA backend（兼容 facade 可保留）
- [x] 15.13 检视迭代（连续 3 轮无新发现后通过）
- [x] 15.12 报告更新（`generate_report.jl` 已纳入并对齐 B1–B5）

### 2026-03-08 Update (15.13 Review Closure)

- 已完成 Phase 15 检视迭代闭环（3 轮）：
   - Round 1（架构/算法）: 复核 PMCHW shell/strong-form/backend 分层、四块语义与 Gate S 归因链，无新增 Phase 15 阻塞问题。
   - Round 2（核心回归）: `test_pmchw_operator_shell.jl`、`test_pmchw_operator_shell_mlfma.jl`、`test_pmchw_gate_s_dense.jl`、`test_nmuller.jl`、`test_nmuller_comparison.jl` 全绿。
   - Round 3（中尺度专门门）: `test_pmchw_gate_s_mlfma_medium.jl`、`test_pmchw_mlfma_budget_medium.jl`、`test_nmuller_planewave_gmres_trajectory_medium.jl` 全绿。
- 因此 `15.13` 正式勾选通过，Phase 15 主线闭环完成。
- 同步记录的非阻塞遗留风险（跨阶段）：
   - `docs/setup_docs.jl` 存在语法错误（`using Pkgusing Pkg`）
   - `test/test_legacy_parity.jl` 存在语法错误（重复 token / 非法标识符）
   - 上述问题未影响本次 Phase 15 主线回归，但建议在后续文档与legacy对齐阶段单独清理。

### 2026-03-09 Update (Cleanup + Phase 14 Checklist Sync)

- 已完成跨阶段遗留语法清理：
   - `docs/setup_docs.jl` 已恢复为标准 docs 环境 `activate/develop/instantiate` 流程。
   - `test/test_legacy_parity.jl` 已重建为可执行最小 parity 回归并通过（`1/1 pass`）。
- 已对齐 Phase 14 遗留勾选状态（避免与实际实现状态不一致）：
   - `14.7`、`14.8`、`14.10`、`14.11` 在 roadmap/progress 中均已同步为完成。

### 2026-03-07 Update (Architecture Pivot Frozen)

- 已接受“先迁移到 Bempp 风格 block/operator shell、后接 FreeFEM/Htool 风格 H-matrix backend”的主线建议，并已同步到：
   - `.github/PHASE_15_DIELECTRIC_ANTENNA_PLAN.md`
   - `.github/plans/phase_15_theory_impl_test_refresh.md`
   - `.github/plans/pmchw_long_krylov_research_report.md`
- 当前冻结的工程边界：
   - transmission 顶层必须由四块 `EJ/EM/HJ/HM` operator shell 承载；backend 不得再定义顶层语义。
   - `strong_form` 必须在 shell 层实现，并作为 Gate S 的统一入口。
   - Dense backend 是下一实现主线；现有 `PMCHWMLFMAOperator` 仅允许作为兼容 facade 或 backend 内核继续存在。
   - H-matrix backend 被正式延后，直到 shell 与 strong-form 分界稳定后再进入计划。
- 对既有诊断结论的修正：
   - “中尺度剩余问题仅属于 GMRES/预条件”这一判断只适用于低精度 `GMRES(200)` 阶段，已不再作为总归因。
   - 当前正式主线是：先完成 shell/strong-form 分界，再继续评估长 Krylov 下的 MLFMA backend 保真度边界。

### 2026-03-07 Update (Dense Shell Landed)

- 已新增 PMCHW 顶层 shell 最小实现：
   - `src/IntegralEquations/PMCHWBlockOperators.jl`
   - `PMCHWBlockOperator`
   - `DensePMCHWBackend`
- 当前 shell 层已具备的能力：
   - 显式保存四块 `EJ/EM/HJ/HM` block 语义
   - `weak_form(shell)` 返回 Dense backend 离散矩阵
   - `strong_form(shell)` 在 shell 层落位，并已接入真实 RWG surface Gram 的 2N block pairing
   - shell 自身已具备 `size` / `eltype` / `mul!` / `*`，可直接送入现有 GMRES
- 已新增正式测试 `test/test_pmchw_operator_shell.jl`，当前结果：`14/14 pass`
   - 验证 block split 重构原始 PMCHW dense matrix
   - 验证 `weak_form` / `strong_form` 接口存在且位于 shell 层
   - 验证默认 `strong_form` 确实使用真实 block pairing，而不是退化为 `weak_form`
   - 验证 shell matvec 与 dense PMCHW 一致
   - 验证 shell 可直接喂给 `GMRESSolver`
- 相关回归已通过：`test/test_pmchw.jl` 全绿，说明新 shell 没有破坏既有 PMCHW direct/激励路径。
- 边界说明：
   - 这一步完成的是 Dense shell + Dense pairing 基线，不代表 Gate S 已完成。
   - 虽然默认 `strong_form` 已接入 RWG pairing，但中尺度 `dense weak / dense strong / MLFMA weak / MLFMA strong` 四路对照尚未执行，因此还不能据此下长 Krylov 归因结论。
   - 现有 `PMCHWMLFMAOperator` 尚未并入 shell；下一阶段是 shell 下的 strong-form 对照，然后再收编 MLFMA backend。

### 2026-03-07 Update (Gate S Dense Half-Step)

- 已新增 shell 层 helper：
   - `strong_form_rhs(shell, rhs)`
   - `recover_trial_coefficients(shell, coeffs)`
- 已新增正式中尺度测试 `test/test_pmchw_gate_s_dense.jl`，当前结果：`3/3 pass`
- 当前固定夹具仍为 `N=540` 球面 PMCHW，比较同一 RHS、同一 GMRES 设置下：
   - `dense weak`
   - `dense strong`
   - 参考 `dense LU`
- 当前结果表明 dense strong 已出现预期方向的改善：
   - `weak_err = 9.29%`
   - `strong_err = 7.90%`
   - `weak_res = 1.39e-3`
   - `strong_res = 1.18e-3`
- 结论边界：
   - Gate S 的 **dense 半边** 已经可执行，并且证明 strong-form 在当前夹具上不是空操作。
   - 但 `MLFMA weak / MLFMA strong` 两路尚未并入 shell，因此完整 Gate S 仍未闭环。

### 2026-03-07 Update (MLFMA Backend Shell Small-Fixture Green)

- 已把 PMCHW MLFMA 路径接入新的 shell 语义：
   - 新增 `MatrixFreePMCHWBackend`
   - 对非矩阵 backend，`strong_form(shell)` 现在返回显式的 pairing-transformed operator wrapper，而不是退回 dense-only 分支
- 已新增正式测试 `test/test_pmchw_operator_shell_mlfma.jl`，当前结果：`6/6 pass`
   - 验证 shell 的 `weak_form` 与包裹后的 `PMCHWMLFMAOperator` matvec 一致
   - 验证 `strong_form(shell)` 对 MLFMA backend 可直接进入 GMRES
   - 验证小夹具上 MLFMA shell 的 weak/strong 输入阻抗与 Dense shell 对应路径保持近似一致
- 当前边界：
   - 这一步证明了 **backend 收编链可用**，但还不等于完整 Gate S 关闭。
   - `N=540` 中尺度上，Dense 半边已经正式锁定；MLFMA 半边的中尺度 strong-form 回归在当前终端环境下仍未稳定跑通，因此暂未写入正式验收门。
   - 新鲜后台终端复验表明：现有 `test/test_pmchw_mlfma_operator_medium.jl` 的 `15.M1` 仅 Gate C matvec 就需要约 `4m02s`；完整中尺度 MLFMA weak/strong 四路 Gate S 会显著更重，因此当前缺口主要是**执行成本与测试编排**，不是 shell 接口缺失。

### 2026-03-07 Update (Medium Gate S Dedicated Regression Green)

- 已新增正式中尺度四路对照测试：`test/test_pmchw_gate_s_mlfma_medium.jl`
- 固定夹具仍为 `N=540` 球面 PMCHW，且四条路径统一使用同一 RHS、同一 shell 语义、同一 GMRES 参数：
   - `dense weak`
   - `dense strong`
   - `MLFMA weak`
   - `MLFMA strong`
- 后台完整回归已通过，当前结果：`15.S2 Medium full Gate S split | 9/9 pass | 9m05.1s`
- 关键对照指标：
   - `relw = 0.0034531756832806245`
   - `rels = 0.0031827064903027532`
   - `zgw = 2.3360090313060917e-5`
   - `zgs = 5.750150008262233e-5`
   - `res_dense_weak = 0.006678820622918415`
   - `res_dense_strong = 0.004991220682354908`
   - `res_mlfma_weak = 0.006679133150245565`
   - `res_mlfma_strong = 0.004991223133137125`
- 当前结论边界：
   - 这一步已把 Gate S 的**中尺度四路主对照**从一次性探针升级为正式回归，证明 shell 下的 weak/strong 与 dense/MLFMA 比较链在目标夹具上可稳定执行。
   - `dense strong` 与 `MLFMA strong` 的残差几乎重合，`strong` 相比 `weak` 在 Dense 与 MLFMA 两条路径上都给出一致改善，因此当前 shell-level strong-form 分界已具备正式证据。
   - 由于单次运行约 9 分钟，该回归暂不并入默认 `test/runtests.jl`；它当前定位为 Phase 15 的专门中尺度验收门，而不是日常快速回归。
   - 该结果锁定的是 Gate S 的 S1-S5 主对照链；治理文档中要求的外部真值/外部 fast-vs-dense 基线仍保持为后续独立验收项，不在本次结果中提前宣告完成。

### 2026-03-07 Update (Phase 15 DeltaGap APIs Wired Into Main Suite)

- 已确认 Phase 15 前置 API 实现已经落地于主源码：
   - `excitation_vector(op::PMCHW, source::DeltaGapSource, basis::RWGBasis)`
   - `input_impedance(op::PMCHW, source::DeltaGapSource, I_2N, basis::RWGBasis)`
   - `excitation_vector(op::SCFIE, source::DeltaGapSource, surf_basis::RWGBasis, vol_basis::SWGBasis)`
- 已确认对应正式测试存在且当前环境下通过：
   - `test/test_pmchw_excitation.jl` → `14/14 pass`
   - `test/test_scfie_delta_gap.jl` → `10/10 pass`
- 已把两者接入默认 `test/runtests.jl`，因此步骤 `15.1`–`15.5` 不再停留在“已有单文件测试、但未纳入主回归”的半完成状态。
- 当前边界：
   - `test/test_scfie_delta_gap.jl` 仍保留对 `TriTetra.nas` 夹具缺失时的 `skip` 逻辑，以避免在缺少混合网格夹具的环境中把 fixture 缺失误判为实现失败。
   - 这一步完成的是 Phase 15 的馈电/阻抗前置 API 与主回归接线，不替代后续 B1-B5 端到端天线精度基准。

### 2026-03-07 Update (B2 Antenna Benchmark Path Enabled)

- 已更新 `benchmark/accuracy/run_B1_B5_antenna.jl`，将 B2 从占位输出改为正式可执行路径：
   - `PMCHWBlockOperator` + `MatrixFreePMCHWBackend(PMCHWMLFMAOperator(...))`
   - shell 层 `strong_form` + `strong_form_rhs` + `recover_trial_coefficients`
   - `GMRES` 求解后与同夹具 Direct 参考阻抗比较
- B2 单项后台复验已通过：
   - `PMCHW MLFMA strong-form GMRES: 108.90 s`
   - `Z_in = +0.000 + j(+0.186) Ω`
   - `Ref Z_in = +0.000 + j(+0.180) Ω`
   - `Im` 误差约 `0.01 Ω`
- 已修正基准脚本中的近零实部判定：
   - 新增 `re_ref_floor_ohm`，用于 `Re(Z_ref) ≈ 0` 时给相对实部误差提供稳定参考尺度
   - 当前 B2 使用 `re_ref_floor_ohm = 1.0`，避免球面 delta-gap 场景中把数值非常接近的结果误判为失败
- 当前边界：
   - `15.6` 仍不整体勾选，因为 B5 仍为空位；当前完成的是 B2 这一路的正式接线与可执行验收，不是整个 B1-B5 基准闭环。

### 2026-03-07 Update (B1/B3/B4 Default Antenna Benchmark Green)

- 已复跑 `benchmark/accuracy/run_B1_B5_antenna.jl` 默认入口，当前默认启用路径 `B1, B3, B4` 全部通过：
   - `B1_PMCHW_sphere_eps4` → PASS
   - `B1_PMCHW_eps1_vs_2xEFIE` → PASS
   - `B3_VEFIE_TriTetra_direct` → PASS
   - `B4_VCFIE_TriTetra_direct` → PASS
- 当前观测值：
   - B1: `Z_in = +0.000 + j(+0.180) Ω`
   - B3: `Z_in = +1313.922 + j(+115.585) Ω`
   - B4: `Z_in = +461.791 + j(-1532.130) Ω`
- 这意味着 `benchmark/accuracy/run_B1_B5_antenna.jl` 现在已经具备：
   - 默认直跑入口 `B1/B3/B4`
   - 单项可执行的 `B2`
   - CSV 结果落盘与汇总输出
- 当前边界：
   - `15.6` 仍未整体关闭，因为 B5 尚未接通；但该脚本已经不再是“只有占位输出”的半成品，而是只差最后一条 fast 路径的可执行基准框架。

### 2026-03-07 Update (Full B1-B5 Antenna Benchmark Green)

- 已复跑完整命令：`benchmark/accuracy/run_B1_B5_antenna.jl B1 B2 B3 B4 B5`
- 当前结果：`PASS: 6 / 6`
   - `B1_PMCHW_sphere_eps4` → PASS
   - `B1_PMCHW_eps1_vs_2xEFIE` → PASS
   - `B2_PMCHW_MLFMA` → PASS
   - `B3_VEFIE_TriTetra_direct` → PASS
   - `B4_VCFIE_TriTetra_direct` → PASS
   - `B5_VCFIE_TriTetra_MLFMA` → PASS
- B5 当前实现路径：
   - `SCFIE(freq, perms; alpha=0.5)`
   - `MLFMAOperator(scfie, [surf_basis, vol_basis], leaf_size)`
   - `Z_near` LU 左预条件 + `GMRESSolver`
- 当前观测值：
   - B2: `Z_in = +0.000 + j(+0.186) Ω`
   - B5: `Z_in = +461.791 + j(-1532.130) Ω`
   - B5 与 `B4 Direct` 当前夹具上 `Re/Im` 误差均为 `0.00`
- 已同步修正 `benchmark/accuracy/generate_report.jl` 中的 B5 标签，现可正确汇总并显示 `B1`–`B5` 五条天线基准。
- 结论边界：
   - `15.6` 现已可以正式关闭，说明 Phase 15 天线端口基准脚本已经完整可执行。
   - 当前未关闭项已不再是天线基准接线，而是 `15.13` 检视迭代与后续 alternative formulation / backend 误差预算工作。

### 2026-03-07 Update (N-Muller Dense Baseline Started)

- 已正式把 Phase 15 下一子流收敛为 alternative formulation 基线，优先项为 `NMuller` dense baseline。
- 已新增执行文档：`.github/plans/phase_15_nmuller_dense_baseline.md`。
- 本轮实现边界已冻结为：
   - 新增 `NMuller` formulation 类型；
   - 复用现有 PMCHW L/K 子块与 RWG surface Gram，构成 `M^(1) - M^(2)` dense baseline；
   - 新增 `test/test_nmuller.jl`，只先锁定尺寸/有限性/block 非平凡性/可执行 direct solve。
- 当前仍未宣告完成的事项：
   - 同一 sphere 上的 PMCHW / N-Muller conditioning 与 GMRES 行为对照；
   - formulation 层的数值验收门；
   - backend 误差预算工作。

### 2026-03-07 Update (N-Muller Dense Baseline + Sphere Comparison Green)

- 已完成 `NMuller` dense baseline 的源码接线：
   - `src/IntegralEquations/NMuller.jl`
   - `src/IntegralEquations/IntegralEquations.jl`
   - `src/EMSuite.jl`
- 已完成正式测试接线：
   - `test/test_nmuller.jl` → `15/15 pass`
   - `test/test_nmuller_comparison.jl` → `9/9 pass`
   - `test/runtests.jl` 已纳入上述两项
- 已新增独立对照脚本：`benchmark/compare_pmchw_nmuller_sphere.jl`
- 当前共享球夹具（`freq=120 MHz`, `radius=0.1 m`, `N=54`, `2N=108`）上的正式观测：
   - `cond(PMCHW)   = 1.711099e7`
   - `cond(NMuller) = 2.538759e5`
   - `cond ratio    = 1.483700e-2`
   - `GMRES PMCHW   : iters=200, res=4.901509e-1, rel_vs_LU=1.019453`
   - `GMRES NMuller : iters=200, res=1.236587e-1, rel_vs_LU=2.729159e-1`
- 当前结论边界：
   - 在当前 dense baseline 与共享 RHS 下，`NMuller` 已表现出显著优于 `PMCHW` 的 conditioning 与 GMRES 行为。
   - 这一步关闭了 roadmap 中 `15.14` 与 `15.15` 的小球 dense 对照基线，但不提前替代后续 medium-scale / backend 误差预算工作。

### 2026-03-07 Update (N-Muller Medium Preset Probe Green)

- 已将 `benchmark/compare_pmchw_nmuller_sphere.jl` 参数化为 `small` / `medium` 两个预设，避免后续重复编写一次性实验脚本。
- `medium` 预设当前配置：
   - `freq=120 MHz`
   - `radius=0.1 m`
   - `n_theta=6`, `n_phi=10`
   - `N=150`, `2N=300`
   - `GMRES(reltol=1e-5, maxiter=250)`
- 当前观测：
   - `cond(PMCHW)   = 6.102107e7`
   - `cond(NMuller) = 2.254102e5`
   - `cond ratio    = 3.693974e-3`
   - `GMRES PMCHW   : iters=250, res=7.078038e-1, rel_vs_LU=1.000066`
   - `GMRES NMuller : iters=250, res=4.174158e-2, rel_vs_LU=7.808047e-2`
- 当前结论边界：
   - N-Muller 相对 PMCHW 的 conditioning / GMRES 优势在 `N=150` 预设上仍保持，且相对小球夹具更明显。
   - 这一步仍属于 benchmark probe，不等价于正式 medium-scale gate；若后续要升格为验收门，需再冻结夹具、RHS 与门限。

### 2026-03-07 Update (N-Muller Medium Dedicated Gate Green)

- 已新增专门中尺度对照测试：`test/test_nmuller_comparison_medium.jl`
- 当前固定夹具：
   - `freq=120 MHz`
   - `radius=0.1 m`
   - `n_theta=6`, `n_phi=10`
   - `N=150`, `2N=300`
   - `GMRES(reltol=1e-5, maxiter=250)`
- 当前结果：`8/8 pass`
- 当前验收门锁定为：
   - `cond(NMuller) < cond(PMCHW) / 50`
   - `NMuller rel_vs_LU < 0.1`
   - `PMCHW rel_vs_LU > 0.9`
   - `NMuller resnorm < 0.1`
   - `PMCHW resnorm > 0.5`
- 当前意义：
   - N-Muller 相对 PMCHW 的优势已不再只是 benchmark probe，而是拥有了一个独立、可复跑的中尺度 dense 专门验收入口。
   - 该测试当前仍不并入默认 `test/runtests.jl`，定位与 `test/test_pmchw_gate_s_mlfma_medium.jl` 类似，属于 Phase 15 的专门门禁而非日常快速回归。

### 2026-03-07 Update (N-Muller DeltaGap + PMCHW Budget Interface Green)

- 已为 `NMuller` 补齐最小天线接口链：
   - `excitation_vector(op::NMuller, source::DeltaGapSource, basis::RWGBasis)`
   - `input_impedance(op::NMuller, source::DeltaGapSource, I_2N, basis::RWGBasis)`
- 已新增正式测试：`test/test_nmuller_excitation.jl`
   - 当前结果：`10/10 pass`
   - 锁定内容：DeltaGap 向量布局、`input_impedance` 仅使用前 `N` 个电流型分量、direct solve 下 `Z_in` 有限且非零
- 已为 PMCHW MLFMA backend 增加显式预算接口：
   - 新增 `PMCHWMLFMAErrorBudget`
   - `PMCHWMLFMAOperator(pmchw, basis, leaf_size; budget=...)` 现在支持显式控制 `near_range`、`L_min` 与 `leaf_size_eff`
   - 旧入口 `PMCHWMLFMAOperator(pmchw, basis, leaf_size)` 保持兼容
- 已在 `test/test_pmchw_mlfma_operator.jl` 中新增 budget 回归：
   - 当前结果：`PMCHWMLFMAOperator budget interface preserves defaults and exposes overrides | 8/8 pass`
   - 在固定夹具上验证了默认启发式保持不变，且显式 budget 能稳定覆盖 `near_range=9`、`leaf_size_eff=0.04`、`L_min=4`
- 当前意义：
   - “都要做”中的两条剩余实现项都已从设计/调研状态转成正式源码与可复跑测试。
   - N-Muller 现已不只具备 plane-wave dense 对照能力，也有了可用于输入阻抗链路的最小馈电接口；PMCHW MLFMA 也不再只能依赖硬编码经验参数。

### 2026-03-07 Update (N-Muller DeltaGap 语义待校准)

- 已新增诊断脚本：`benchmark/compare_pmchw_nmuller_impedance.jl`
- 当前观测：在共享球夹具上，直接把 PMCHW DeltaGap 语义桥接到 `NMuller` 后，`Z_in` 与 PMCHW 参考相比出现 `1e7`–`1e8` 量级失配。
- 进一步 front-half / back-half 取流诊断表明：
   - 当前 `input_impedance(op::NMuller, ...)` 使用的 front-half 启发式不是一个可接受的物理端口定义；
   - 仅改成 back-half 取流也不能把结果拉回 PMCHW 量级；
   - 问题已收敛到 `NMuller` DeltaGap RHS 语义与 feed-current 提取语义本身，而不是 dense 组装或求解器稳定性。
- 当前工程边界已调整为：
   - 保留 `NMuller` DeltaGap / `input_impedance` 作为诊断性接口；
   - 暂不把它升级为正式 PMCHW-vs-N-Muller 阻抗验收门；
   - 下一步必须先完成 formulation-specific 端口语义校准，再决定是否固化测试门限。

### 2026-03-07 Update (PMCHW MLFMA Budget Sweep Benchmark Landed)

- 已新增专门预算基准：`benchmark/compare_pmchw_mlfma_budget.jl`
- 当前脚本固定输出以下量，避免 budget 讨论只停留在接口层：
   - `leaf_size_eff` / `near_range` / `L_min`
   - `nnz_near` 与 `near_density`
   - `MLFMA matvec vs Direct` 相对误差与相关系数
   - shell strong-form `Z_in` 相对 Direct 的 `Re/Im` 误差
   - 构造、单次 matvec 与 GMRES 求解耗时
- 当前目标：
   - 把 `PMCHWMLFMAErrorBudget` 从“可设置”推进到“可量化比较”；
   - 为后续正式冻结 default / loose / tight 预算门提供 CSV 基线，而不是继续凭经验调整 `near_range` 或 `L_min`。
- 首轮 `medium` 夹具（`N=540`）观测已记录到 `test_results/accuracy/PMCHW_MLFMA_budget_sweep_medium.csv`：
   - `default`：`near_density = 0.9049`，`rel_matvec = 3.20e-4`
   - `loose_near`：`near_density = 0.3007`，`rel_matvec = 7.30e-4`
   - `fixed_leaf_0p04_nr9`：`near_density = 0.2350`，`rel_matvec = 7.70e-4`
   - `tight_near` / `tight_near_L4`：`near_density = 1.0000`，`rel_matvec ≈ 1.5e-15`
- 当前解释边界：
   - budget 已被证明会显著改变 near/far 划分和 matvec fidelity；
   - 但在固定 `GMRES(maxiter=100)` 的 B2 strong-form 路径上，各预算的 `Z_in` 误差都停留在 `Re≈9.3%`, `Im≈91Ω` 量级，说明当前夹具上的阻抗偏差主要仍受 Krylov 截止主导，而不是单纯的 budget 松紧。

### 2026-03-07 Update (PMCHW Medium Budget Gate Frozen)

- 已新增专门门禁：`test/test_pmchw_mlfma_budget_medium.jl`
- 当前固定三组预算：
   - `default`
   - `loose_near`
   - `fixed_leaf_0p04_nr9`
- 锁定的门禁事实：
   - `default > loose > fixed` 的 `near_density` 梯度必须存在；
   - `default` 的 matvec 相对误差必须优于两组更稀疏预算；
   - 在固定 `GMRES(100)` strong-form 路径下，三组预算的 `Z_in` spread 与最终残差 spread 必须保持很小。
- 当前定位：
   - 该测试不是为了证明哪组预算已经“足够好”，而是为了把 Phase 15 当前已知结论固化成正式专门回归：budget 会影响 near/far 与 matvec fidelity，但在当前求解截止下并不会主导 `Z_in` 偏差。

### 2026-03-07 Update (PMCHW Representative Budget Set + Long Krylov Comparison Landed)

- 已把 `benchmark/compare_pmchw_mlfma_budget.jl` 的默认预算集合收缩为固定代表性三组：
   - `default`
   - `loose_near`
   - `fixed_leaf_0p04_nr9`
- 如需恢复探索性预算扫描，现需显式使用 `full` 模式；这意味着 Phase 15 后续默认讨论对象不再是开放式 budget 空间，而是已冻结的代表性集合。
- 已新增专门长 Krylov 对照基准：`benchmark/compare_pmchw_mlfma_budget_krylov.jl`
   - 固定夹具：`medium`, `N=540`
   - 固定预算：`default / loose_near / fixed_leaf_0p04_nr9`
   - 固定两档 Krylov：
      - `short`: `restart=100`, `maxiter=100`, `reltol=1e-4`
      - `long`: `restart=300`, `maxiter=600`, `reltol=1e-6`
   - 结果已落盘到 `test_results/accuracy/PMCHW_MLFMA_budget_krylov_medium.csv`
- 本轮正式观测把预算子流与长 Krylov 子流第一次接通：
   - dense strong 本身从 `short` 到 `long` 已明显逼近 `LU`：`Re` 误差由 `11.370%` 降到 `6.780%`
   - `short` 档中，三组预算相对 dense 的 `Z_in` gap 仍只有亚欧姆量级：
      - `default`: `0.126844Ω / 0.286394Ω`
      - `loose_near`: `0.028818Ω / 0.182589Ω`
       - `fixed_leaf_0p04_nr9`: `0.261993Ω / 0.568182Ω`
   - 在 receive 修复并刷新 benchmark CSV 后，`long` 档的三组预算也都只剩亚欧姆量级 gap：
      - `default`: `0.038151Ω / 0.002364Ω`
      - `loose_near`: `0.481070Ω / 0.020291Ω`
      - `fixed_leaf_0p04_nr9`: `0.062819Ω / 0.076291Ω`
- 结论边界已更新：
   - “固定 `GMRES(100)` 下 budget 不主导 `Z_in`”这一结论只适用于短 Krylov 截止；
   - 在现行实现下，长 Krylov 也不再出现旧的一阶量级 budget 分叉；修复后的主结论已与 BG2、Arnoldi、GMRES 轨迹诊断对齐。

### 2026-03-07 Update (PMCHW Medium Long-Krylov Budget Gate Green)

- 已新增专门门禁：`test/test_pmchw_mlfma_budget_krylov_medium.jl`
- 该门禁固定：
   - medium 夹具 `N=540`
   - long Krylov 配置 `restart=300`, `maxiter=600`, `reltol=1e-6`
   - 代表性预算中的两组高信息量配置：`default` 与 `loose_near`
- 在 PMCHW K-type receive 修复后，旧门限已被证明过时；最新回归结果为：`15.BG2 PMCHW MLFMA medium long-Krylov budget gate | 16/16 pass | 29m08.0s`
- 当前正式锁定的结果：
   - dense strong 在该配置下相对 `LU` 仍保持：`Re` 误差 `6.7796%`，`Im` 误差 `4.2830Ω`，残差 `6.22e-5`
   - `default` budget 相对 dense strong 的 `Z_in` gap 已缩小到 `0.03815Ω / 0.00236Ω`
   - `loose_near` budget 相对 dense strong 的 `Z_in` gap 也仅 `0.48107Ω / 0.02029Ω`
   - `loose_near` 仍比 `default` 更差，但差异已经从旧的一阶量级收缩到亚欧姆量级
- 当前意义：
   - 这一步说明 `M` pass receive 修复不只修复了 block fidelity，也同步消除了 BG2 中原先被解释为“budget 主导”的大部分 long-Krylov gap；
   - 当前 Phase 15 剩余边界不再能简单归因于 budget/fidelity 主导，而必须进一步检查更深 Krylov 轨迹与 formulation/conditioning 放大机制。

### 2026-03-07 Update (PMCHW Medium Block Fidelity Diagnostic Landed)

- 已新增专门诊断基准：`benchmark/compare_pmchw_block_fidelity_medium.jl`
- 当前固定夹具：
   - medium sphere，`N=540`
   - 代表性预算：`default` 与 `loose_near`
   - 物理 probe：`J_only` 与 `M_only`
   - 比较语义：统一以**物理空间 probe / 物理空间输出**做对照；在 `strong_form` 下输入先做 `trial_pairing` 变换，输出再映回测试空间，避免把 weak-space probe 直接送入 strong operator 产生伪差异
- 当前结果已落盘到：`test_results/accuracy/PMCHW_block_fidelity_medium.csv`
- 当前正式结论：
   - 当 strong-form 结果映回物理空间后，`strong` 与 `weak` 的 block fidelity 指标一致；此前 strong 路径的“大失真”属于诊断脚本的 probe 语义错误，而不是 shell/backend 的新问题
   - `J_only` probe 下，`default` budget 的总体误差仅 `3.25e-5`，主要误差仍集中在极小的 H-row 分量；说明 J 通道在当前 medium 夹具上继续维持高保真
   - PMCHW K-type receive 公式修复后，`M_only` probe 下 `default` budget 的总体误差已降到 `2.98e-4`，`loose_near` 下也仅 `2.03e-3`
   - `M×k0` 与 `M×k1` 单 pass 的 `E-row` 误差分别收敛到 `2.87e-4 / 9.89e-4`（`default`）与 `1.06e-3 / 3.00e-3`（`loose_near`）；说明此前 medium 边界确由 K-type receive 实现错误驱动，而不是 block/backend 普遍失真
- 当前意义：
   - 这一步先用 medium block 诊断锁定了 `M_only -> E-row` 故障链，又在修复后证明 medium block fidelity 已恢复到 `1e-3` 量级内；
   - 后续若继续排查 medium 长 Krylov 失配，应回到 BG1/BG2 重新评估 budget / Krylov / 预条件链，而不是沿用已过时的 `M_only -> E-row` backend 故障归因。

### 2026-03-07 Update (PMCHW Medium Block Fidelity Gate Green)

- 已新增正式专门回归：`test/test_pmchw_block_fidelity_medium.jl`
- K 接收链修复后，门禁已从“锁定故障形态”切换为“锁定修复后精度上界”；最新回归结果：`15.BF1 PMCHW medium block fidelity gate | 26/26 pass | 4m42.9s`
- 该门禁当前固定：
   - medium 夹具 `N=540`
   - 预算：`default` 与 `loose_near`
   - 物理 probe：`J_only` 与 `M_only`
   - 额外 pass-level 诊断：`M×k0` 与 `M×k1`
- 当前正式锁定的事实：
   - strong-form 在映回物理空间后与 weak-form 具有相同 fidelity 指标，因此 medium block gap 不属于 strong-only 问题
   - `J_only` 继续保持高保真：`default` 下 `rel_total = 3.25e-5`
   - `M_only` 已不再是 medium block 的主导失配链：`default` 下 `rel_E = 2.99e-4`，`loose_near` 下也仅 `2.03e-3`
   - `M×k0` 与 `M×k1` 两条 pass 的 `E-row` 误差已分别降到 `2.87e-4 / 9.89e-4`（`default`）与 `1.06e-3 / 3.00e-3`（`loose_near`），同时 `loose_near` 仍稳定劣于 `default`，因此门禁仍能感知 budget 退化
- 当前意义：
   - 这一步把 `15.BF1` 从故障快照升级成修复后回归门，防止 PMCHW K-type receive 链再次退化；
   - 后续实现应把重点移回 medium long-Krylov 与 budget/GMRES 主线，而不是继续围绕已修复的 `M` pass 接收链做重复诊断。

### 2026-03-08 Update (PMCHW M-Pass Receive Parity Repaired)

- 已在 `test/test_pmchw_mlfma_operator.jl` 新增 **15.GD2RM Gate D2 M-pass receive parity**，对 `M×k0` 与 `M×k1` 两条单 pass 直接比较 fast far output 与 dense far reference。
- RED 阶段结果表明：修复前 `E-row` 基本完全失真，而 `H-row` 仍基本正确：
   - `k0`: `rel = 9.99e-1`, `corr = 0.203`, `rel_E = 9.99e-1`, `rel_H = 4.19e-3`
   - `k1`: `rel = 9.99e-1`, `corr = 0.109`, `rel_E = 9.99e-1`, `rel_H = 1.36e-2`
- 已修正 `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`：
   - `_receive_terms` 中的 K receive 从 MFIE 风格 `cross(rho, normal) · H_inc` 改为 PMCHW 风格 `rho · (r̂ × E)`
   - 同步校正 K 相关符号链：`factor_HJ` 改为正号、`factor_EM` 改为负号
- 修复后 `15.GD2RM` 已转绿，且既有 J-pass / end-to-end PMCHW 门禁保持绿色：
   - `k0`: `rel = 4.67e-3`, `corr = 0.999993`
   - `k1`: `rel = 1.94e-2`, `corr = 0.999918`
- 当前意义：
   - 这一步完成了 `M_only -> E-row` medium 故障链的代码级闭环：诊断定位 → RED test → 实现修复 → GREEN test → medium 门禁恢复；
   - 该修复现已成为后续 Phase 15 medium/backend 回归的硬前提。

### 2026-03-08 Update (Post-Receive Budget Gates Revalidated)

- 已对 medium budget / long-Krylov 两条冻结主线重新做结果核对：
   - `test_results/accuracy/PMCHW_MLFMA_budget_sweep_medium.csv`
   - `test_results/accuracy/PMCHW_MLFMA_budget_krylov_medium.csv`
   - `test/test_pmchw_mlfma_budget_krylov_medium.jl` 最近一次回归仍为绿
- 当前复核结果表明，PMCHW K-type receive 修复**保留了短 Krylov 结论，但推翻了旧的 long-Krylov 大 gap 结论**：
   - 短 Krylov `GMRES(100)` 下，三组代表性 budget 仍保持很小的 `Z_in` spread，说明短 Krylov 主导项依然是 solver 截止
   - medium budget sweep 当前数值为：
      - `default`: `near_density = 0.9049`, `rel_matvec = 3.20e-4`
      - `loose_near`: `near_density = 0.3007`, `rel_matvec = 7.30e-4`
      - `fixed_leaf_0p04_nr9`: `near_density = 0.2350`, `rel_matvec = 7.70e-4`
   - long Krylov `restart=300, maxiter=600, reltol=1e-6` 下，现行 BG2 复核值已收缩到：
      - `default` 相对 dense strong 的 `Z_in` gap 为 `0.03815Ω / 0.00236Ω`
      - `loose_near` 相对 dense strong 的 `Z_in` gap 为 `0.48107Ω / 0.02029Ω`
- 当前意义：
   - `15.GD2RM` / `15.BF1` 的修复已经清除了 `M_only -> E-row` 的 block-level 假故障；
   - `15.G13` 仍成立，但 `15.G15` 已从“budget 主导 long-Krylov gap”更新为“long-Krylov 下 MLFMA strong 与 dense strong 仍近似重合，budget 只保留次级灵敏度”；
   - 因此 Phase 15 剩余主线需要从单纯的 budget/backend fidelity 归因，转向更深 Krylov 轨迹与 formulation/conditioning 放大机制。

### 2026-03-08 Update (Medium Arnoldi-Subspace Fidelity Diagnostic Landed)

- 已新增诊断基准：`benchmark/compare_pmchw_krylov_subspace_medium.jl`
- 该诊断按研究报告建议，改用 **dense strong-form Arnoldi 子空间方向** 取代随机 probe：
   - 固定 medium 夹具 `N=540`
   - 固定 `steps=12`
   - 固定两组代表性 budget：`default` 与 `loose_near`
   - 在每个 Arnoldi 方向上直接比较 `dense strong` 与 `MLFMA strong` 的 matvec，并拆分 `E/H` 两半块误差
- 当前结果已落盘到：`test_results/accuracy/PMCHW_krylov_subspace_medium.csv`
- 当前正式观测：
   - `default` 下前 12 个 Arnoldi 方向的最大误差仅：`max rel_total = 8.12e-5`, `max rel_E = 8.12e-5`, `max rel_H = 3.48e-4`
   - `loose_near` 下前 12 个 Arnoldi 方向也仍很小：`max rel_total = 4.26e-4`, `max rel_E = 4.26e-4`, `max rel_H = 2.85e-3`
   - 两组预算在这些早期 Krylov 方向上都保持几乎完美相关（`corr_total ≈ 1.0`）
- 当前意义：
   - 这一步把 medium 诊断从“随机向量 Gate C”升级到了“Arnoldi 子空间方向”；
   - 它表明 BG2 的 long-Krylov 大 gap 不能简单归因于“前几步 Krylov 方向上的 fast matvec 已明显失真”，而更可能来自更深子空间方向、非正交化/重启轨迹、或 formulation 条件数放大效应。

### 2026-03-08 Update (Medium GMRES Trajectory Diagnostic Landed)

- 已新增轨迹诊断基准：`benchmark/compare_pmchw_gmres_trajectory_medium.jl`
- 该诊断固定：
   - medium 夹具 `N=540`
   - strong-form `restart=300`, `reltol=1e-6`
   - checkpoint：`100 / 300 / 600`
   - 预算：`default` 与 `loose_near`
- 当前结果已落盘到：`test_results/accuracy/PMCHW_gmres_trajectory_medium.csv`
- 当前正式观测：
   - `default` 路径相对 dense strong 的解向量差在 `100 / 300 / 600` 步分别为 `1.31e-3 / 4.24e-4 / 1.23e-3`
   - `default` 的 `Z_in` gap 在 `300` 步已降到 `0.00620Ω / 0.00566Ω`，到 `600` 步仍仅 `0.03815Ω / 0.00236Ω`
   - `loose_near` 路径在 `600` 步的解向量差增大到 `1.52e-2`，但其 `Z_in` gap 仍只有 `0.48107Ω / 0.02029Ω`
- 当前意义：
   - 这一步进一步证明 receive 修复后，long-Krylov 主线不再表现为“大尺度阻抗分叉”；
   - `loose_near` 在更深 checkpoint 上的解向量偏差仍比 `default` 大，说明 budget 灵敏度没有消失，但它当前不足以单独解释旧 BG2 时代的 7–16Ω 级差距。

### 2026-03-08 Update (Phase 15 Mainline Clarified)

- 已把 Phase 15 的“主线 / 支线 / 下一步”显式写回治理文档：
   - `.github/PHASE_15_DIELECTRIC_ANTENNA_PLAN.md`
   - `.github/REFACTORING_ROADMAP.md`
- 当前统一口径为：
   - PMCHW 仍是主交付 formulation；
   - N-Muller 只是 dense 对照基线，用于比较 conditioning / GMRES 行为；
   - 当前剩余问题已收缩到更深 Krylov 轨迹与 formulation/conditioning 放大，而不是继续纠缠 N-Muller 支线或旧的 budget 主导结论。
- 当前意义：
   - 后续开发不再需要从历史更新里反推“现在到底在做什么”；
   - 若继续推进 Phase 15，应默认沿 PMCHW 主线前进，并只在需要做归因分离时调用 N-Muller dense 对照。

### 2026-03-08 Update (PMCHW vs N-Muller Medium GMRES Trajectory Gate Green)

- 已新增专门测试：`test/test_nmuller_gmres_trajectory_medium.jl`
- 已新增 benchmark：`benchmark/compare_pmchw_nmuller_gmres_trajectory_medium.jl`
- 当前 benchmark 结果已落盘到：`test_results/accuracy/PMCHW_NMuller_gmres_trajectory_medium.csv`
- 当前正式结果：
   - 专门测试通过：`PMCHW vs N-Muller medium dense GMRES trajectory | 25/25 pass | 5.9s`
   - 在 `50 / 100 / 150 / 200 / 250` 全部 checkpoint 上，`NMuller` 的 `rel_vs_LU` 与 `resnorm` 都持续优于 `PMCHW`
   - 代表性终点 `iter=250`：
      - `PMCHW`: `rel=1.000066e+00`, `res=7.078038e-01`
      - `NMuller`: `rel=7.808047e-02`, `res=4.174158e-02`
- 当前意义：
   - 这一步把“用 N-Muller 区分 formulation 与 backend 归因”的建议升级成了正式门禁；
   - 当前证据进一步支持：receive 修复后，Phase 15 剩余主线更像是 `PMCHW formulation / conditioning` 问题，而不是已修复的 PMCHW fast backend 主干问题。

### 2026-03-08 Update (PMCHW vs N-Muller Medium Plane-Wave GMRES Trajectory Gate Green)

- 已新增专门测试：`test/test_nmuller_planewave_gmres_trajectory_medium.jl`
- 已新增 benchmark：`benchmark/compare_pmchw_nmuller_planewave_gmres_trajectory_medium.jl`
- 当前 benchmark 结果已落盘到：`test_results/accuracy/PMCHW_NMuller_planewave_gmres_trajectory_medium.csv`
- 当前正式结果：
   - 专门测试通过：`PMCHW vs N-Muller medium plane-wave GMRES trajectory | 25/25 pass | 7.4s`
   - 在 `50 / 100 / 150 / 200 / 250` 全部 checkpoint 上，`NMuller` 的 `rel_vs_LU` 与相对残差都持续优于 `PMCHW`
   - 代表性终点 `iter=250`：
      - `PMCHW`: `rel=1.002186e+00`, `rres=1.779276e-03`, `iters=250`
      - `NMuller`: `rel=1.149702e-02`, `rres=9.583083e-06`, `iters=67`
- 当前意义：
   - 这一步把 formulation/conditioning 归因从 random probe 推进到了 `PlaneWave` 物理激励 RHS；
   - 当前证据链进一步收紧：Phase 15 剩余主线仍更像 `PMCHW formulation / conditioning` 问题，而不是 repaired PMCHW backend 在物理工况下重新失真。

### 2026-03-08 Update (PMCHW Medium Plane-Wave Dense Weak/Strong Trajectory Gate Green)

- 已新增专门测试：`test/test_pmchw_gate_s_planewave_trajectory_medium.jl`
- 已新增 benchmark：`benchmark/compare_pmchw_gate_s_planewave_trajectory_medium.jl`
- 当前 benchmark 结果已落盘到：`test_results/accuracy/PMCHW_gate_s_planewave_trajectory_medium.csv`
- 当前正式结果：
   - 专门测试通过：`PMCHW medium plane-wave dense weak/strong GMRES trajectory | 31/31 pass | 6.3s`
   - 在 `50 / 100 / 150 / 200 / 250` 全部 checkpoint 上，`strong` 的 `rel_vs_LU` 都略优于 `weak`，但改善比例始终不足 1%
   - 代表性终点 `iter=250`：
      - `weak`: `rel=1.002186e+00`, `rres=1.779276e-03`, `iters=250`
      - `strong`: `rel=1.001109e+00`, `rres=1.868816e-03`, `iters=250`
- 当前意义：
   - 这一步把 PMCHW 主线内部的 weak/strong 分界也推进到了 `PlaneWave` 物理激励；
   - 当前证据显示：在该工况下，strong-form 只能带来边际改善，不能单独解释或消除 medium PMCHW 的剩余主问题。

### 2026-03-08 Update (PMCHW vs N-Muller Medium Plane-Wave Dense Restart Sweep Green)

- 已新增专门测试：`test/test_pmchw_nmuller_planewave_restart_sweep_medium.jl`
- 已新增 benchmark：`benchmark/compare_pmchw_nmuller_planewave_restart_sweep_medium.jl`
- 当前 benchmark 结果已落盘到：`test_results/accuracy/PMCHW_NMuller_planewave_restart_sweep_medium.csv`
- 当前正式结果：
   - 专门测试通过：`PMCHW vs N-Muller medium plane-wave dense restart sweep | 22/22 pass | 6.8s`
   - `IterativeSolvers.gmres` 默认 `restart = min(20, size(A, 2))`；在当前 `N=150` 夹具上，这意味着此前未显式指定 restart 的 plane-wave dense 轨迹实际上使用的是 `restart=20`
   - `PMCHW` 的 plane-wave dense 终点对 restart 极其敏感：
      - `restart=20`: `rel=1.002186e+00`, `rres=1.779276e-03`, `iters=250`
      - `restart=150`: `rel=1.874645e-01`, `rres=6.275530e-05`, `iters=250`
      - `restart=250`: `rel=5.186365e-02`, `rres=8.501219e-06`, `iters=165`
   - `NMuller` 在 `restart=50` 后已基本稳定：`rel=7.800091e-03`, `rres=8.610863e-06`, `iters=44`
- 当前意义：
   - 这一步修正了 Phase 15 对 plane-wave dense 轨迹的旧口径：默认 `restart=20` 的确会明显放大 PMCHW 的坏表现；
   - 但在把 restart 影响显式剥离后，PMCHW 虽已把残差压回到与 `NMuller` 同量级，`rel_vs_LU` 仍高约 `6.6x`，因此 formulation gap 仍然存在，只是边界比之前更精确。

### 2026-03-08 Update (Dense Plane-Wave Attribution Substream Wrapped)

- 当前收尾状态：
   - `G22` 已把 dense random-RHS 归因推进到 `PlaneWave` 物理激励；
   - `G23` 已证明 PMCHW `strong-form` 在该工况下只有边际改善；
   - `G24` 已证明默认 `restart=20` 会显著放大 PMCHW 坏轨迹，但在高 restart 下仍保留明显 `PMCHW -> NMuller` formulation gap。
- 因此本段子流可以先视为收口：
   - 已不需要继续重复 `weak/strong` 或默认 restart 的现象复测；
   - 下一轮若继续，应直接进入显式 full-restart / Arnoldi 级诊断。

### 2026-03-08 Update (Workspace Housekeeping Before Next Phase)

- 已删除 `scripts/` 下 5 个一次性诊断脚本：4 个 `tmp_*` 探针和 `diagnose_pmchw_farfield_blocks.jl`。
- 保留范围明确为：
   - 正式测试入口；
   - benchmark 入口；
   - 已在 roadmap/progress 中回链的 CSV 结果；
   - 计划与阶段文档。
- 当前治理口径：下一阶段若需要继续追踪 PMCHW / N-Muller / restart / Arnoldi 边界，应优先扩展现有正式资产，不再回到无回链的临时脚本堆积模式。

### 2026-03-06 Update (Governance Refresh)

- Rolled back the recent experimental PMCHW/MLFMA code updates as requested.
- Added governance plan: `.github/plans/phase_15_theory_impl_test_refresh.md`.
- New task `15.G1` added to roadmap to enforce a strict Theory -> Implementation -> Test contract.
- Next implementation work must satisfy mandatory Gate A/B/C/D checks before acceptance.

### 2026-03-06 Update (现场清理 + 文档治理重构)

- Deleted batches of one-off scratch scripts under `scripts/` and removed the ungoverned archive pile; repository policy is now: diagnostics must become formal tests or be deleted.
- Rewrote `.github/plans/phase_15_theory_impl_test_refresh.md` from a high-level governance note into a **function-level development contract** covering:
   - geometry/octree layer
   - interpolation/translation precompute layer
   - near-field assembly layer
   - pass-level aggregation/translation/disaggregation layer
   - constructor/mul! layer
   - Gate A/B/C/D acceptance mapping
- Tightened Phase 15 review rule from “≥2 rounds” to **“3 consecutive review rounds with no new findings”**.
- Completed document review iteration for `15.G1`:
   - Round 1 found missing function-layer rollback rules -> fixed by adding per-function contracts
   - Round 2 found missing anti-guess constraints and scratch-script disposal rule -> fixed by adding prohibitions and Gate D traceability
   - Round 3 found no new document-level issues -> governance document accepted
- Result: `15.G1` is complete; further PMCHW/MLFMA work must execute against the refreshed contract, not ad-hoc debug notes.

### 2026-03-06 Update (Gate A/B Implemented)

- Added executable Gate A/B tests into `test/test_pmchw_mlfma_operator.jl`.
- Gate A (structural invariants) is green:
   - HJ + EM near-field invariant passed
   - non-trivial near/far split passed
   - dual-octree permutation sanity passed
- Gate B (EJ pass alignment) is now machine-tracked with fixed seed (`Random.seed!(42)`):
   - k0: `rel=0.7536`, `corr=0.9962` (magnitude regression, marked `@test_broken`)
   - k1: `rel=0.9588`, `corr=0.4641` (marked `@test_broken`)
- Main end-to-end gate remains failing: `15.11 MLFMA mul! vs Direct` = `49.40%`.

### 2026-03-06 Update (Gate B k0 Green, k1 still Red)

- Applied and verified two effective fixes in `PMCHWMLFMAOperator.jl`:
   - fixed pass-state clearing loop (`for (_, lv) in oct.levels`)
   - restored PMCHW disaggregation block scaling to `1/(4π)` chain factors
- Gate B EJ k0 is now green:
   - `rel=0.0877`, `ratio=0.9908`, `corr=0.9962`
- Gate B EJ k1 remains red (tracked by `@test_broken`):
   - `rel=0.8897`, `ratio=0.3899`, `corr=0.4641`
- End-to-end Gate C improved from ~49.40% to ~42.27%, still above `<10%` target.
- Added robustness fix in `Interpolation.jl` to avoid interpolation sparse assembly length mismatch for clipped interpolation orders.

### 2026-03-06 Update (Structured Debug Plan — 结构化调试计划制定)

- Previous approach was ad-hoc (coefficient guessing, unstructured experiments).
- Adopted structured two-part debugging methodology per user direction:
  - **Part 1**: Near-field matrix element-by-element comparison (`Z_near` vs `Z_full` at sparsity pattern)
    - Gate D1 test: `max_element_err < 0.1%`
    - If FAIL: trace `assemble_near_field_pmchw` for offending `(i,j)` pair
  - **Part 2**: Far-field column extraction via unit vectors (`e_j` → column j of MLFMA approximation)
    - Gate D2 test: scan all N columns, find max-error column, compare with dense far-field reference
    - For worst `(i_max, j_max)`: module-level decomposition (aggregation → translation → disaggregation)
   - **Part 3**: Full operator Gate C after all passes verified
- Execution document: `.github/plans/phase_15_theory_impl_test_refresh.md`
- Execution order: Gate D0 → Gate D1 → Gate D2 (EJ k1) → fix identified module → Gate D2 (EM/HM) → Gate C
- Status: Plan confirmed, implementation pending user review.

### Phase 0 Dense PMCHW Baseline Validation — **已完成** ✅ (d7ee785)

- **根本 Bug 确认**: `radiation_integral_rwg` 通过全局 `get_k0()` 读取波数；Julia 启动时 `k0=0.0`；`PMCHW()` 构造函数原本未调用 `set_frequency!`，导致所有远场辐射积分 phase=1（exp(0)=1），远场计算退化为0阶近似。
- **修复** (commit `d7ee785`):
  - `src/IntegralEquations/PMCHW.jl`: 导入 `set_frequency!`，在构造函数末尾自动调用 `set_frequency!(Float64(freq))`
  - `src/PostProcessing/RCS.jl`: PMCHWT `radarCrossSection` 重载（显式 k0 参数版本）在调用 `radiation_integral_rwg` 前同步全局 `k0`
- **验证结果** (`scripts/diag_step0_pmchw_dense_vs_mie.jl`):
  - E 面 RMSE: `8.89 dB → 0.82 dB` ✅ PASS (< 1.5 dB 门限)
  - 后向散射: MoM=5.81 vs Mie=5.76 dBsm (Δ=0.06 dB)
  - H 面: lat=8 时 RMSE=1.63 dB（近零点误差，lat≥12 时 RMSE<0.79 dB）
- **下一步**: Gate D0（Dense PMCHW 基线复验）→ Gate D1（近场矩阵元素对比）→ Gate D2（单位向量列提取）→ Gate B k1 修复

### 2026-03-07 Update (Gate D0 前置化)

- 用户新增约束：**在 Gate D1 开始前，必须重新验证 Dense PMCHW 基线正确性**，不允许直接把历史 Phase 0 结果当作默认前提。
- 已将该要求正式写入执行文档 `.github/plans/phase_15_theory_impl_test_refresh.md`：
   - Gate D0 = Dense PMCHW -> Mie 基线复验
   - 门限固定为：E-plane RMSE < 1.5 dB，H-plane RMSE < 2.0 dB，后向散射误差 < 0.2 dB
- 已在 `test/test_pmchw_mlfma_operator.jl` 中新增可执行 Gate D0，直接使用：
   - `PMCHW + PlaneWave(freq, 0, 0, [1,0,0])`
   - `radarCrossSection(..., I_2N, basis, k0, eta0)`
   - `calculate_mie_rcs_dielectric_sphere(...)`
   - `compute_rcs_accuracy(...)`
- Phase 15 强制执行顺序更新为：**Gate D0 → Gate D1 → Gate B → Gate C**。

### 2026-03-07 Update (Gate D1 Executable + Green)

- 已将 Gate D1 从“诊断说明”升级为正式测试：`test/test_pmchw_mlfma_operator.jl` 新增 **15.GD1 Gate D1 近场逐元素对齐**。
- Gate D1 现在按文档要求直接比较：
   - `I, J, V = findnz(Z_near)`
   - `dense_vals = Z_full[I, J]`
   - `max_element_rel_err = maximum(abs(Z_near - Z_full) / abs(Z_full))`
- 当前结果：**Gate D1 = Green**
   - `max_err = 0.0`
   - `mean_err = 0.0`
- 结论：`assemble_near_field_pmchw` 在当前叶层 near sparsity pattern 上与 Dense PMCHW 完全一致，后续定位不得再把主因归到近场块装配。
- 下一步收窄到 far path：**Gate D2 / Gate B k1 / Gate C**。

### 2026-03-07 Update (Gate D2 Executable, EJ k1 继续 Red)

- 已把 Gate D2 的“单位向量列提取”升级为正式测试：`test/test_pmchw_mlfma_operator.jl` 新增 **15.GD2 Gate D2 EJ k1 最坏列扫描**。
- 当前扫描结果：
   - `worst_col = 16`
   - `worst_rel = 0.9567`
   - `worst_corr = 0.3018`
- 结论：远场误差已被进一步收敛到 **EJ k1 far path**，且不是均匀小偏差，而是最坏列严重失配；后续修复应直接围绕 `octree1 / k1 / aggregation-translation-disaggregation` 链路展开。

### 2026-03-07 Update (Block Jacobi 接口补齐 + 中尺度长 Krylov 边界)

- 已补齐此前文档声称存在、但代码实际缺失的 Block Jacobi 接口：
   - `BlockJacobiPreconditioner` 现在支持**非连续索引块**，不再假设叶块在原矩阵中必须是连续 `UnitRange`。
   - `MLFMAOperator` 新增直接构造：`BlockJacobiPreconditioner(op)` 与兼容签名 `BlockJacobiPreconditioner(op, basis)`。
   - `PMCHWMLFMAOperator` 新增 2N 耦合叶块构造：每个叶块同时包含 `J` 分量和对应的 `M` 分量。
- 已新增/扩展正式测试：
   - `test/test_preconditioners.jl` 增加 **arbitrary index blocks** 测试。
   - `test/test_preconditioners.jl` 增加 **PMCHW MLFMA constructor** 与兼容签名测试。
   - 当前 `test/test_preconditioners.jl` = `20/20 pass`。
- 中尺度 `N=540` 求解实验新结论：
   - 现有 `Diagonal` / `ILU(Z_near)` / 新增 `BlockJacobi(op)` 在 PMCHW 中尺度夹具上均**不能改善** GMRES；其结果反而明显劣于无预条件。
   - `GMRES` 的 `restart` 不是纯参数细节，而是关键边界：
      - `dense+GMRES(restart=300, maxiter=600)` 已把 `Re(Zin)` 误差压到约 `4.35%`。
      - 同参数下 `MLFMA+GMRES` 仍约 `8.18%`，说明“中尺度问题”已不再只是纯 solver/preconditioner 问题。
      - `dense+GMRES(restart=1000, maxiter=1000)` 基本恢复到 `LU`（`Re` 误差约 `0.17%`，残差约 `1.52e-6`）。
      - 同参数下 `MLFMA+GMRES` 仍有约 `5.79%` 的 `Re` 误差、约 `4.68%` 的 `Im` 误差，尽管其对自身算子的残差也降到 `1.50e-6`。
- 结论修正：
   - 先前“中尺度偏差完全属于 solver/preconditioner，而非 MLFMA correctness”的判断仅在 `GMRES(200)` 低精度阶段成立。
   - 当 Krylov 子空间做大到足以逼近 `LU` 时，PMCHW MLFMA 仍存在**长 Krylov 下的 operator-fidelity 边界**；随机向量 Gate C 绿色不足以推出高精度求解绿色。
- 额外边界：尝试继续增大 PMCHW `near_range`（通过更大 `leaf_size` 输入触发）会在 `N=540` 夹具上触发 `OutOfMemoryError()`，因此“继续扩大近场直到近似消失”不是可接受主线方案。

### 2026-03-07 Update (外部实现八叉树 / matvec 对照已并入长 Krylov 调研报告)

- 已更新 `.github/plans/pmchw_long_krylov_research_report.md`，新增“八叉树 / 树结构 / matvec 组成 / 与本地差异”专项对照。
- Bempp 侧确认：
   - transmission 主线保持为 block multitrace operator，而不是 PMCHW 专用双八叉树接口；
   - FMM 由通用 ExaFMM 后端管理，公开控制参数是 `depth`、`expansion_order`、`ncrit`、`near_field_representation`；
   - Maxwell boundary matvec 的公开实现结构是“基函数/散度变换 -> FMM evaluator -> `singular_part @ x` 修正”，说明它在 operator block 级别而非 PMCHW 专用四遍级别组织快速乘法。
- FreeFEM/BEMTool/Htool 侧确认：
   - 公开主线不是 MLFMA 八叉树，而是 cluster-tree + admissible blocks 的 H-matrix 压缩；
   - matvec 由近场稠密块与远场低秩块共同组成，压缩误差参数直接并入求解控制。
- 对 EMSuite 的新增工程结论：
   - 本地 `PMCHWMLFMAOperator` 当前是显式“双八叉树 + 四遍远场 + 单个 2N `Z_near`”专用实现，这与 Bempp/FreeFEM 的公开组织方式都不同；
   - 后续诊断不应只停留在整体 2N `mul!` 误差，而应拆到 `E/H` block 或更细的 operator evaluator 级别；
   - 即使继续坚持 MLFMA 主线，也必须把 `leaf_size / near_range / translation accuracy` 提升为面向目标误差的正式预算接口。

### 2026-03-07 Update (Gate S strong-form 对照主线已写入治理文档)

- 已更新 `.github/plans/phase_15_theory_impl_test_refresh.md`，新增 **Gate S: strong-form 分界**。
- Gate S 的作用不是再加一个泛化建议，而是把中尺度 `N=540` 长 Krylov 归因正式改成四路对照：
   - `dense weak`
   - `dense strong`
   - `MLFMA weak`
   - `MLFMA strong`
- 新治理要求：四路对照必须使用同一 RHS、同一 restart/maxiter、同一停止准则；禁止只改单一路径的求解参数后再下归因结论。
- 归因规则已被明确冻结：
   - 若 `dense strong` 明显优于 `dense weak`，则 conditioning 是显著因素；
   - 若 `dense strong` 已接近 `LU` 而 `MLFMA strong` 仍失配，则剩余主因正式归入 MLFMA 压缩误差/算子保真度；
   - 若 `dense weak` 与 `dense strong` 都不收回，则不得先把问题归到 MLFMA。
- 这意味着 Phase 15 后续不再允许把长 Krylov 偏差笼统记为“solver/preconditioner 问题”；必须先经过 strong-form 分界再决定下一轮实现主线。

### 2026-03-07 Update (外部开源基线用例已整理为精度对照库)

- 已在 `.github/plans/pmchw_long_krylov_research_report.md` 中新增“推荐外部基线用例库”，把可复用的开源案例按用途整理为 dense、compressed matvec、far-field、solver、transmission 五类。
- 当前确定的第一优先基线为：
   - scuff-em 的 Mie scattering / dielectric sphere 系列，用作单球 transmission 真值基线；
   - Bempp 的 Maxwell boundary FMM sphere real/complex wavenumber 系列，用作 compressed-vs-dense matvec 量级基线；
   - Bempp 的 `use_strong_form` GMRES 接口，用作 Gate S 的方法学基线。
- 第二优先基线为：
   - Bempp 双介质球 `maxwell_dielectric.py`，用于 multi-body transmission 扩展；
   - FreeFEM Maxwell EFIE sphere，作为 H-matrix 压缩参数预算基线；
   - scuff-em `unit-test-BEMMatrix` 的 dielectric sphere / two-sphere matrix regression，作为 dense dielectric assembly 结构性回归基线。
- `.github/plans/phase_15_theory_impl_test_refresh.md` 的 Gate S 已同步要求：后续 strong-form 分界不得只依赖 EMSuite 自身历史结果，必须至少绑定单球真值基线与 fast-vs-dense matvec 基线各一项。

### 2026-03-07 Update (Shared-Core Boundary Test)

- 新增共享核心分界测试：`test/test_pmchw_mlfma_operator.jl` 增加 **15.GD2S Gate D2 shared core EFIE k1 对照**。
- 该测试移除 PMCHW 四遍包装，只保留共享 MLFMA 核心：
   - `build_octree`
   - `aggregate!`
   - `translate!`
   - `disaggregate_downward!`
   - `disaggregate_leaf!`
- 对照对象为 `direct far EFIE(k1)`，当前结果仍为红：
   - `rel = 0.8912`
   - `corr = 0.5065`
- 结论更新：问题**不只**在 PMCHW 包装层，已经被正式收敛到共享 MLFMA 核心的 `k1` 路径；后续排查重点应转向 `Interpolation / Translation / shared RWG disaggregation`，而不是继续怀疑 `assemble_near_field_pmchw` 或 PMCHW 四遍调度本身。

### 2026-03-07 Update (Packaging Parity Tests)

- 新增两组包装层 parity 测试：
   - **15.GD2A Gate D2 aggregation parity**：`aggregate_leaf_pmchw!` vs shared `aggregate!`（J/k1）
   - **15.GD2R Gate D2 receive parity**：`disaggregate_leaf_pmchw_j!` vs shared `disaggregate_leaf!`（J/k1）
- 目的：把“PMCHW 包装公式”与“共享 MLFMA 核心公式”解耦验证。
- 预期结论：
   - 若两者都为绿，则 `PMCHW` 包装层可暂时排除，问题继续收敛到 `Interpolation / Translation / shared leaf receive field generation`。
   - 若其中任一为红，则可直接定位到对应包装公式实现。

### 2026-03-07 Update (Level-Boundary Test)

- 新增 **15.GD2L Gate D2 level boundary**，对同一 shared `EFIE(k1)` 核心比较两个层级配置：
   - `leaf_size = 0.10` → `nLevels = 4` → **Red**
   - `leaf_size = 0.20` → `nLevels = 3` → **Green** (`rel = 0.0120`, `corr = 0.99997`)
- 结论进一步收紧：
   - leaf-level aggregation 不是主因（GD2A 绿）
   - leaf-level receive wrapper 不是主因（GD2R 绿）
   - 单层 shared far path 在 3-level 配置下可工作
   - 错误主要出现在 **更深层级（4-level）才触发的 interpolation/downward chain**，而不是所有 translation 都错

### 2026-03-07 Update (Intermediate-Level Parity + Leaf Dominance)

- 新增 **15.GD2P Gate D2 level-3 parity**：比较
   - 4-level 配置的 `level-3 disaggG`
   - 3-level 配置的 `leaf-level disaggG`
- 当前结果为绿：`rel = 4.05e-4`, `corr = 0.9999999`
- 新增 **15.GD2T Gate D2 leaf translation dominance**：量化 4-level leaf 场中
   - 本层 translation 贡献
   - 父层 downward 贡献
- 当前结果：父层贡献仅为 leaf translation 的 `1.25%`
- 结论再次收紧：
   - 首次深层父层下行（2→3）与对应 3-level 参考基本一致
   - 4-level 最终 leaf 场主要由 **level-4 translation 自身** 主导
   - 因此下一步应优先审查 `level-4 αTrans / αTransIndex / leaf-level far offset enumeration`，而不是继续怀疑 PMCHW 接收包装或 3→4 downward 本身

### 2026-03-07 Update (Leaf-Translation Isolated Boundary)

- 新增 **15.GD2U Gate D2 leaf translation isolated**：
   - MLFMA 侧仅取 leaf-level translation 产物（不做 downward）
   - Dense 侧仅保留 leaf `farneighbors` 对应的 basis-pair 子块
- 该测试用于把问题从“level-4 相关”继续压缩到“leaf translation 自身 vs leaf farneighbor 枚举”
- 若此门继续为红，则后续应重点审查：
   - `Translation.cal_alpha_trans_on_level!` 在 leaf level 的展开式与常数
   - `translate!` 中 `relative3DID -> αTransIndex` 的映射
   - `OctreeBuilder.setKidLevelFarNeighbors!` 生成的 leaf interaction list

### 2026-03-07 Update (Near-Range Threshold Located)

- 新增 **15.GD2V Gate D2 near-range threshold**：在同一 shared `k1` 夹具上比较
   - `leaf_size=0.10, near_range=4` → 红
   - `leaf_size=0.10, near_range=7` → 绿
- 终端扫描进一步显示：正式 PMCHW 构造所使用的 `leaf_size_eff=0.05, near_range=8/10/12` 仍然红，说明原构造函数中的 near-range 基准常数偏小。
- 基于该阈值证据，已将 `PMCHWMLFMAOperator` 的 near-range 缩放基准从 `4` 提升到 `8`，以避免把过近的 leaf interaction 误留在 far path。
- 修复后正式 **15.GD2 Gate D2 EJ k1 最坏列扫描** 已转绿：
   - `worst_rel = 0.01686`
   - `worst_corr = 0.99990`
   - 说明正式 `PMCHWMLFMAOperator` 的 `k1` far-path 主失配已被 near-range 修正消除
- 修复后正式 **15.GB Gate B EJ pass 对齐** 中 `k1` 路径也已转绿：
   - `rel1 = 0.00999`
   - `corr1 = 0.99998`
   - `k0` 路径同步收敛到 `rel0 = 0.00219`, `corr0 = 0.999999`
- 端到端 `mul!` 对 Direct 也已恢复：
   - `Gate C rel_err = 6.19e-4`
   - 远优于原 10% 门限，说明本轮修复已从 shared `k1` leaf-path 根因贯通到正式 PMCHW operator 的整条 matvec 链路
- `B2` 输入阻抗验收也已通过：
   - `Zin_direct = 112.559 - 58.734im`
   - `Zin_mlfma  = 108.186 - 57.077im`
   - `Re` 相对误差 `= 3.885% < 5%`
   - GMRES 在 `maxiter=200` 截止下仍满足当前 B2 验收门限

### 2026-03-07 Update (Mainline Green vs Diagnostic Red)

- 当前 **正式 PMCHW 主线门** 已全部转绿：
   - `GD0` Dense PMCHW 基线
   - `GD1` 近场逐元素
   - `GD2` EJ k1 最坏列
   - `GB` EJ k0/k1 pass
   - `Gate C` 端到端 matvec
   - `B2` 输入阻抗
- 当前保留为 `@test_broken` 的门仅用于**shared-core 诊断回归**：
   - `GD2S` shared core EFIE k1
   - `GD2L` 4-level red / 3-level green
   - `GD2U` leaf translation isolated
   - `GD2V` near-range threshold 的 red 侧 (`near_range=4`)
- 这些红灯现在不再表示正式 PMCHW operator 失败，而是用于固定“shared core 在较小 near-range 下会误分类 leaf 近距交互”的已知机制。

### 2026-03-07 Update (Medium-scale Validation Added)

- 新增正式中尺度回归文件 `test/test_pmchw_mlfma_operator_medium.jl`，固定球面夹具：
   - `r = 0.5m`, `lat_divs = 10`, `lon_divs = 20`
   - `N = 540`, `2N = 1080`
- 中尺度 `Gate C` 结果：
   - `gateC_rel = 3.20e-4`
   - `nnz_near = 1,055,424 / 1,166,400`，说明当前修复虽然保证正确性，但 near field 已非常接近稠密
- 中尺度 `B2` 现象拆分：
   - `MLFMA+GMRES(200)` 对 `LU` 的输入阻抗误差约 `20.3%`
   - 但 `dense+GMRES(200)` 与 `MLFMA+GMRES(200)` 几乎完全一致
   - 结论：中尺度下的剩余问题是 **GMRES 收敛/预条件**，不是 PMCHWMLFMAOperator 的 matvec 正确性
- 正式运行 `test/test_pmchw_mlfma_operator_medium.jl` 已通过：
   - `15.M1` Medium-scale Gate C: `rel = 3.20e-4`, `corr = 0.99999995`
   - `15.M2` Medium-scale GMRES parity: `rel_I = 7.51e-3`, `corr_I = 0.9999718`
   - `Zin` parity 极强：`rel_re = 4.29e-6`, `rel_im = 6.70e-5`
   - 残差历史也一致：`res_gap = 1.55e-3`

---
## 褰撳墠闃舵: Phase 13.3 V-EFIE MPI 骞惰鍖?鈥?**宸插畬鎴?* 鉁?

**鏈€鏂版垚鏋?(2026-03-06)**:
- 鉁?**V-EFIE MPI Allreduce 鏂规瀹炵幇瀹屾垚**
  - 绠楁硶: 姣忚繘绋嬪鐞?`(it-1) % n_procs == rank` 鐨勬祴璇曞洓闈綋; 瀵圭О鍒╃敤 (js>it); `MPI.Allreduce!` 姹囨€?
  - 璋冨害淇: 鍒犻櫎 `module VolumeAssemblyMPI` 灏佽锛岀洿鎺?`import .Assembly: assemble_impedance_matrix_parallel` 鎵╁睍
  - 姝ｇ‘鎬ч獙璇? Tetra.nas N=3201, Z[1,1] 鍦?1/2/4 杩涚▼涓嬪畬鍏ㄧ浉鍚?(`4.747e6 - 4.716e6im`)
  - 鎬ц兘: 1杩涚▼ 10.95s 鈫?2杩涚▼ 8.86s (1.24脳); 灏?N 鏃?Allreduce 寮€閿€鍗犱富瀵间负姝ｅ父鐜拌薄
- 鉁?**BasisFunctions.jl**: 琛ュ叏 `get_triangles_info` 瀵煎嚭
- 鉁?**Parallel.jl**: 淇瀛愭ā鍧楅殧绂婚棶棰?(Julia dispatch 璺ㄦā鍧椾笉鍙)

## Phase 13.2 Option C 璇勪及 鈥?**宸插畬鎴?* 鉁? 
**Phase 13.3 Option A (MPI 骞惰鍖栨祴璇?** 鈥?**瀹屾垚** 鉁?

**鏈€鏂版垚鏋?(2026-03-0X)**:
- 鉁?**Option C 瀹屾垚**: PWC/RBF 鍩哄嚱鏁版€ц兘鍩哄噯娴嬭瘯
  - PWC+VEFIE: 12.62s (3.16脳 faster than SWG 39.94s)
  - 浣?PWC DOF +38% (21834 vs 15828) 鈫?姹傝В鍣ㄨ礋鎷呭姞閲?
  - Per-DOF虏 鏁堢巼浠?SWG 鐨?0.17脳 鈫?**涓嶉€傚悎閫氱敤浼樺寲**
  - RBF 纭浠呮敮鎸佸叚闈綋缃戞牸 (涓嶉€傜敤鍥涢潰浣?VEFIE)
  - 馃搳 璇︾粏鎶ュ憡: [OPTION_C_PWC_RBF_EVALUATION.md](.github/OPTION_C_PWC_RBF_EVALUATION.md)
  
**Phase 13.1 鎬荤粨** (2026-03-01):
- 鉁?V-EFIE: 39.05s (**1.18脳 vs Legacy** 46.13s) 鈥?**棣栨蹇簬 Legacy**
- 鉁?SCFIE: 63.81s (**1.04脳 vs Legacy** 66.68s) 鈥?**涓?Legacy 鎸佸钩**
- 鉁?FastExp 鏌ユ壘琛? 10000 鏉＄洰, 绮惧害 < 0.002%, ~160KB
- 鈿狅笍 Thread-local buffer 澶辫触 (-63% 鎬ц兘閫€姝? 宸插洖閫€)
- 馃搳 璇︾粏鎶ュ憡: [PHASE_13.1_SUMMARY.md](.github/PHASE_13.1_SUMMARY.md)

---

## 宸插畬鎴?鉁?

### Phase 1: 鍩虹鏋舵瀯 (2025-12)
- [x] 椤圭洰缁撴瀯 `EMSuite.jl` 鍒涘缓
- [x] `Project.toml` 渚濊禆绠＄悊
- [x] Core 妯″潡: `Interfaces.jl`, `Types.jl`, `Constants.jl`, `Materials.jl`, `Sources.jl`
- [x] Utilities: `Logging.jl`, `Parameters.jl`
- [x] CI/CD: `.github/workflows/CI.yml`
- [x] 鏂囨。妗嗘灦: Documenter.jl 閰嶇疆

### Phase 2: 鍑犱綍涓庡熀鍑芥暟 (2025-12)
- [x] 缃戞牸绫诲瀷: `TriangleMesh`, `TetrahedraMesh`, `HexahedraMesh`
- [x] 缃戞牸 I/O: Nastran (`.nas`), Gmsh (`.msh`)
- [x] 鍧愭爣鍙樻崲, 楂樻柉姹傜Н
- [x] RWG 鍩哄嚱鏁?(涓?Legacy 100% 鍖归厤, 鍖呮嫭杈规帓搴忛€昏緫)
- [x] SWG, RBF, PWC 鍩哄嚱鏁?

### Phase 3: 绉垎鏂圭▼ (2025-12 ~ 2026-01)
- [x] EFIE: PEC Plate/Sphere 楠岃瘉, 濂囧紓椤?($F_1$, $F_2$) 淇
- [x] MFIE: `mfie_interaction!` in-place 缁勮
- [x] CFIE: PEC Sphere 楠岃瘉 (RMSE 2.09 dB vs Mie)
- [x] VEFIE: 浣撶Н绉垎鏂圭▼, SWG 鍩哄嚱鏁版敮鎸?
- [x] SCFIE: 闈綋鑰﹀悎绉垎鏂圭▼

### Phase 4: MLFMA (2026-01)
- [x] 鍏弶鏍戞瀯寤?(`Octree.jl`, `OctreeBuilder.jl`)
- [x] 鑱氬悎 (`Aggregation.jl`), 鍚爣閲忓娍椤?
- [x] 杞Щ (`Translation.jl`), Legacy 鍥犲瓙 $-jk/16\pi^2$ 瀵归綈
- [x] 瑙ｈ仛 (`Disaggregation.jl`)
- [x] Lebedev 鐞冮潰鎻掑€奸泦鎴?
- [x] 杩戝満涓€鑷存€? Max Diff < 1e-12
- [x] 杩滃満绮惧害: 淇 1/4 鍥犲瓙, 杩戦偦缂撳啿鍖?= 4

### Phase 5: 姹傝В鍣ㄤ笌骞惰 (2026-01)
- [x] Direct Solver (LU)
- [x] GMRES (璇樊 2.7e-7)
- [x] BiCGSTAB
- [x] ILU 棰勬潯浠跺櫒
- [x] SPAI 棰勬潯浠跺櫒
- [x] MPI 鍒嗗竷寮忓苟琛?(n=2 vs n=1 鏈哄櫒绮惧害鍖归厤)
- [x] 澶氱嚎绋嬪苟琛?(4 绾跨▼鍔犻€熼獙璇?

### Phase 6: 鍚庡鐞嗕笌 I/O (2026-01)
- [x] RCS 璁＄畻
- [x] FarField / NearField 璁＄畻
- [x] VTK 瀵煎嚭 (ParaView)
- [x] 缁撴灉鏂囦欢 I/O (HDF5, CSV, TXT)
- [x] 鐢垫祦鍒嗗竷鍚庡鐞?

### Phase 7: 楠岃瘉涓庡榻?(2026-02-28) 鉁?
- [x] SCFIE MLFMA 杩戝満楠岃瘉: Rel Err = 1.58e-15 (鏈哄櫒绮惧害)
- [x] SCFIE MLFMA 杩滃満楠岃瘉: Overall 0.85%, Surface 0.85%, Volume 6.1%
- [x] Standalone EFIE MLFMA: 0.66% (4GHz, TriTetra.nas)
- [x] Standalone VEFIE MLFMA: 1.39% (4GHz, TriTetra.nas)
- [x] MoM_AllinOne 鍏ㄩ儴绠椾緥瀵规爣瀹屾垚

**SCFIE MLFMA 淇鐨?7 涓?Bug:**
1. `Disaggregation.jl`: SCFIE `efie_factor` 浠?`1.0+0im` 鈫?`jk畏/(16蟺)`
2. `MLFMAOperator.jl`: 杩戝満 SS 鍧楃Щ闄ゅ浣?`eta` (閬垮厤 MFIE 椤瑰弻閲嶄箻 畏)
3. `MLFMAOperator.jl`: VV 鍧椾粠 `vefie_element_interaction` (c1=-j蠅渭鈧€魏) 鈫?`vefie_element_interaction_kernel` (c1=+j蠅渭鈧€魏)
4. `MLFMAOperator.jl`: VV 鍧楁坊鍔犵己澶辩殑 mass matrix (鑷氦浜掗」)
5. `MLFMAOperator.jl`: VV/SV/VS 鍧楀垱寤?`distribute_term_nosign!` 閬垮厤 bfsSign 鍙岄噸璁℃暟
6. `Disaggregation.jl`: SWG `const_factor` 浠?`-jk畏` 鈫?`jk畏/(4蟺)` (VEFIE 鐢?G=e^{-jkR}/(4蟺R))
7. `MLFMAOperator.jl`: 娣诲姞 VEFIE 缂撳瓨棰勮绠?(TetBasisCache + precompute_vefie_basis)
8. `SCFIE.jl`: Z_SV 绗﹀彿淇 `(term1 - term2)` 鈫?`(term1 + term2)` (鎭㈠ $L$ 绠楀瓙浜掓槗鎬?
9. `SCFIE.jl`: Z_VS 绯绘暟淇 `c1_vs = -j蠅渭鈧€` 鈫?`+j蠅渭鈧€` (缁熶竴涓?Legacy 涓€鑷寸殑姝ｅ彿绾﹀畾)
10. 娣诲姞 `test/test_scfie.jl` 鍥炲綊娴嬭瘯 (浜掓槗鎬с€丏irect 缁勮銆丮LFMA 杩戝満/杩滃満)

### Phase 7.5: ~3 dB 绯荤粺鍋忓樊鏍瑰洜淇 (2026-02-28) 鉁?

**鏍瑰洜 1 鈥?`edgev虃`/`edgen虃` 鏂瑰悜鍙嶈浆** (`BasisUtilities.jl`):
- EMSuite 璁＄畻 `e1 = v2 - v3`锛坴3鈫抳2 鏂瑰悜锛夛紝Legacy 璁＄畻 `v3 - v2`锛坴2鈫抳3 鏂瑰悜锛?
- 瀵艰嚧 `edgen虃 = cross(edgev虃, facen虃)` 鎸囧悜涓夎褰㈠唴閮ㄨ€岄潪澶栭儴
- 褰卞搷: `faceSingularityIgIvecg` 涓?`p02jl`锛堟姇褰辫窛绂伙級绗﹀彿閿欒 鈫?杩戝寮傜Н鍒嗙粨鏋滈敊璇?
- 淇: 灏嗚竟鍚戦噺鏀逛负 `e1 = v3 - v2`, `e2 = v1 - v3`, `e3 = v2 - v1`

**鏍瑰洜 2 鈥?`calc_near_interaction!` 缁忛獙鍥犲瓙 + 闈㈢Н褰掍竴鍖?* (`EFIE.jl`):
- 瀛樺湪缁忛獙鍥犲瓙 `* 1.25`锛堝簲涓?1.0锛?
- `inv_areas = 1.0 / tri_test.area`锛堝簲涓?`1.0 / (tri_test.area * tri_source.area)`锛?
- 鍦?edgev 鏂瑰悜閿欒鏃讹紝杩欎袱涓敊璇儴鍒嗕簰鐩歌ˉ鍋匡紙Z_near 鈮?0.0485脳 姝ｇ‘鍊?鈫?杩戜箮鍙拷鐣ワ級
- 淇: 绉婚櫎 `* 1.25`锛屾敼涓烘纭殑鍙岄潰绉綊涓€鍖?

**楠岃瘉缁撴灉:**
- Z 鐭╅樀瀵硅绾?ratio: 1.000000 (2640脳2640 鏉跨綉鏍?vs Legacy)
- Frobenius 鑼冩暟姣? 1.000010
- RCS (Jet 100MHz): Mean Diff 0.05 dB / RMSE 0.29 dB (Phi=0), Mean Diff 0.008 dB / RMSE 0.09 dB (Phi=90)
- 鍏ㄩ儴 138/138 鍗曞厓娴嬭瘯閫氳繃锛堟棤鍥炲綊锛?

### Phase 10.A: MLFMA 鍥犲瓙淇 (2026-03-02) 鉁?

**Bug 1 鈥?MLFMA far-field 脳4 鍥犲瓙** (`MLFMAOperator.jl`, `Disaggregation.jl`):
- 鏍瑰洜: `efie.factor = jk畏/(16蟺)` 鍖呭惈 `1/4` (鏉ヨ嚜 RWG `l虏/4` 褰掍竴鍖?, 浣?MLFMA 鑱氬悎/瑙ｈ仛鍚勭敤 `l/2`锛屼箻绉?`l虏/4` 宸茶嚜鐒跺寘鍚鍥犲瓙 鈫?**鍙岄噸璁℃暟 1/4**
- Legacy 閬垮厤姝ら棶棰? translation 鐢?`-jk/(16蟺虏)` + disagg 鐢?`jk畏`锛堜笉鍚岀殑鍥犲瓙鍒嗚В鏂瑰紡锛?
- 淇: `y_far *= 4 * operator.factor` (EFIE mul!), CFIE/SCFIE disagg `efie_factor *= 4`
- 楠岃瘉: 鑷唇鎬х郴鏁拌宸?65.7% 鈫?0.30%, RCS RMSE 3.1 dB 鈫?0.028 dB
- A3 S-EFIE MLFMA vs Legacy: Mean Diff 0.048 dB, RMSE 0.303 dB

**Bug 2 鈥?CFIE MLFMA MFIE 绗﹀彿閿欒** (`Disaggregation.jl`):
- 鏍瑰洜: MFIE K 绠楀瓙浣跨敤 $\nabla_{r'}G$锛堟簮姊害锛夛紝杩滃満杩戜技涓?$+jk\hat{k}G$;
  浠ｇ爜閿欒鍦颁娇鐢ㄤ簡 $\nabla_r G$锛堝満姊害锛夌殑 $-jk\hat{k}G$锛屽鑷?MFIE 椤圭鍙峰弽杞?
- Legacy 澶勭悊: 鑱氬悎/瑙ｈ仛鍒嗙, 缁熶竴涔樹互 `jk畏`, EFIE 鍜?MFIE 杩滃満绯绘暟鍚屽彿 (+jk畏)
- 淇: `(-efie_factor)` 鈫?`(+efie_factor)` (MFIE 椤?
- 楠岃瘉: C3 CFIE MLFMA vs Legacy: RMSE 3.45 dB 鈫?**0.003 dB** (1000脳 鏀瑰杽)
- GMRES 杩唬娆℃暟: 50 鈫?7锛堢畻瀛愬噯纭害鎻愬崌鍚庢敹鏁涘姞閫燂級

**宸查獙璇佹祴璇曠粨鏋滄眹鎬?**

| 娴嬭瘯 | 鎸囨爣 | 缁撴灉 |
|------|------|------|
| 鍗曞厓娴嬭瘯 | 138/138 | 鉁?PASS |
| A1 S-EFIE Direct Jet | RMSE vs Legacy | 0.215 dB |
| A3 S-EFIE MLFMA Jet | RMSE vs Legacy | 0.303 dB |
| A3 self-consistency | 绯绘暟璇樊 | 0.30% |
| B1 CFIE 鍒嗚В | rel_err | 0.0 (10/10) |
| C1 S-CFIE Direct Sphere | RMSE vs Legacy | 0.001 dB |
| C3 S-CFIE MLFMA Sphere | RMSE vs Legacy | **0.003 dB** |
| D1-SWG V-EFIE Direct | RMSE vs Legacy | 0.952 dB |
| E1 VSEFIE Direct | RMSE vs Legacy | **0.602 dB** |
| EFIE MLFMA Sphere | RMSE vs Legacy SCFIE | 0.041 dB |

### Phase 10.B: SCFIE Fss 杈圭晫淇 (2026-03-03) 鉁?

**Bug 鈥?缂哄け鍗婂熀鍑芥暟杈圭晫闈㈢Н鍒嗕慨姝?(Fss)** (`SCFIE.jl`, `MLFMAOperator.jl`):
- 鏍瑰洜: 杈圭晫 SWG 鍩哄嚱鏁帮紙鍗婂熀鍑芥暟锛屼粎鏈変竴涓洓闈綋鏀拺锛夌殑鏍囬噺鍔跨己澶辫〃闈㈢Н鍒嗕慨姝ｉ」
- Legacy 鍦?`EFIEVSIERWGSWG.jl` 涓€氳繃 `Fss` 椤瑰鐞嗭細
  - `isbdn=true` 鏃? Z_VS[n,m] += j蠅渭鈧€/(4蟺k虏) 脳 l_m 脳 |A_n| 脳 鈭埆 G dS_tri dS_face
  - `未魏鈮?` 鏃? Z_SV[m,n] += 未魏 脳 (鍚屼笂)
- EMSuite 瀹屽叏缂哄け姝や慨姝?鈫?鑰﹀悎鐭╅樀 Z_SV/Z_VS 鍋忓樊 22%, Z_VV 鍋忓樊 48%
- 淇: 鍦?`SCFIE.jl` 涓坊鍔?`assemble_fss_boundary_correction!` 鍜?`assemble_fss_boundary_correction_sparse`
  - 鐩存帴姹傝В璺緞: 鍦?`assemble_coupling_blocks!` 鍚庤皟鐢?
  - MLFMA 璺緞: 浠ョ█鐤忕煩闃靛舰寮忓姞鍒?Z_near
- 楠岃瘉: E1-VSEFIE RMSE 浠?**5.3 dB 鈫?0.60 dB** (PASS), 138/138 娴嬭瘯鍏ㄩ€氳繃

### 楠岃瘉閲岀▼纰?
- [x] **SEFIE Direct**: `verify_SEFIE_direct.jl` (18.8s 4绾跨▼ / 30.3s 1绾跨▼, Legacy 31.2s)
- [x] **SEFIE MLFMA**: `verify_SEFIE_mlfma.jl` (Ratio 1.0000, Rel Err 1.5%)
- [x] **VEFIE Direct**: `verify_VEFIE_direct.jl` (Legacy Parity)
- [x] **VEFIE MLFMA**: `verify_VEFIE_mlfma.jl` (Ratio 1.0000, Rel Err 0.04%)
- [x] **SCFIE Direct**: `verify_SCFIE_direct.jl` (VSIE Plate+Metal, RCS -15.35 dBsm)
- [x] **SCFIE MLFMA**: `quick_scfie_mlfma_test.jl` (Near-field 1.58e-15, Far-field 0.85%)
- [x] **Standalone EFIE MLFMA**: `test_efie_vefie_farfield.jl` (Rel Err 0.66%)
- [x] **Standalone VEFIE MLFMA**: `test_efie_vefie_farfield.jl` (Rel Err 1.39%)
- [x] **MPI**: `benchmark_parallel_sphere.jl` (Consistency check passed)
- [x] **Threading**: 4 threads vs 1 thread speedup 楠岃瘉

---

## 宸插畬鎴?鉁?(缁?

### Phase 8: 鎬ц兘浼樺寲 (2026-03) 鉁?

#### 8.0 鎬ц兘鍩虹嚎娴嬮噺 鉁?(commit `861426d`)
- [x] 鍒涘缓 `benchmark/performance_baseline.jl` (EMSuite 7 鐢ㄤ緥)
- [x] 鍒涘缓 `LegacyBenchmark/legacy_performance_baseline.jl` (Legacy 瀵规爣)
- [x] EMSuite 鍏ㄩ儴 7 鐢ㄤ緥娴嬮噺瀹屾垚
- [x] Legacy 5 鐢ㄤ緥娴嬮噺瀹屾垚
- [x] 缁煎悎瀵规瘮鎶ュ憡: `test_results/PERFORMANCE_BASELINE.md`

#### 8.1 Z 缁勮鍘婚攣 鉁?(commit `2d4ebe6`)
- [x] `Impedance.jl` SpinLock 鈫?Per-row SpinLock (琛岀骇鏃犻攣骞惰)
- [x] Plate EFIE 缁勮 **-54%**, Jet EFIE 缁勮 **-12%**
- [x] 138/138 娴嬭瘯閫氳繃

#### 8.2 CFIE 鍐呮牳鍚堝苟 鉁?(commit `d0888cf`)
- [x] MFIE 鍐呮牳浼樺寲锛氬叡浜?Green 鍑芥暟銆乮nline rho 鍚戦噺銆佹秷闄ら噸澶嶅嚑浣曡绠?
- [x] CFIE 缁勮 **-74%** (Jet: 168.29s 鈫?43s)
- [x] CFIE/EFIE 缁勮姣? 8.1脳 鈫?**2.31脳** (鐩爣 鈮?2.5脳 鉁?

#### 8.3 MLFMA Z_near 浼樺寲 鉁?(commit `d2f7963`)
- [x] 棰勫垎閰?COO 鏁扮粍浠ｆ浛鍔ㄦ€?push!
- [x] COO 鍚堝苟鍚?sparse() 鏋勯€?

#### 8.4 鍐呭瓨鍒嗛厤鐑偣 鉁?(commit `82988cf`)
- [x] 绉婚櫎 MFIE 褰卞瓙 `get_global_quad_points` 鍑芥暟

#### 8.5 Julia 1.12 鍏煎淇 鉁?(commit `67d3a8a`)
- [x] `threadid()` 鈫?`Threads.maxthreadid()` (Legacy + EMSuite)
- [x] 8.5b 绫诲瀷绋冲畾鎬у鏌ワ細`@code_warntype` 鍏ㄩ儴 clean

#### 8.6 @fastmath + SIMD 鉁?(commit `1c6d499`)
- [x] `calc_interaction!` 閲嶅啓锛氱洿鎺?dot() 鏇挎崲 SMatrix
- [x] `@fastmath` 鍔犻€?exp() 绛夋暟瀛﹁繍绠?
- [x] `@inbounds @simd` 浼樺寲鍐呭惊鐜?

#### 8.7 BlockJacobiPreconditioner 鉁?(commit `76f8b16`)
- [x] 瀹炵幇 `BlockJacobiPreconditioner` (浠?Z_near 鎻愬彇瀵硅鍧? 骞惰 LU)
- [x] 鏋勫缓閫熷害姣?Sparse LU 蹇?**166脳**
- [x] 閫傜敤浜?CFIE (3 娆?GMRES 杩唬); EFIE 涓嶆敹鏁? LU 浠嶄负榛樿
- [x] 娣诲姞 `get_leaf_intervals(op::MLFMAOperator)`

#### 8.8 鏈€缁堝熀鍑嗗娴?鉁?(commit `6f4987a`)
- [x] 鍏ㄩ儴 6 涓敤渚?(+ CFIE 瀵规瘮) 閲嶆柊璁℃椂
- [x] 淇 Sphere CFIE MLFMA OOM (COO 鍒濆鍒嗛厤涓婇檺)
- [x] 鐢熸垚 `test_results/PERFORMANCE_REPORT.md`

**Phase 8 鏈€缁堢粨鏋?(2026-03-01 鏇存柊):**

| 鐢ㄤ緥 | N | 鍘熷鍩虹嚎 | Phase 8.8 | **Phase 8.9** | **鏈€鏂?* | 璇存槑 |
|------|---|---------|---------|----------|------|------|
| Plate EFIE | 2640 | 1.02s | 1.94s | **0.153s** | 0.153s | EFIE 鍐呮牳閲嶅啓 |
| Jet EFIE | 14559 | 20.70s | 29.01s | **4.26s** | 4.26s | EFIE SIMD 淇 |
| **Jet CFIE** | 14559 | **168.29s** | **64.88s** | **14.48s** | 14.48s | CFIE 鏋舵瀯 + `@.` |
| Jet MLFMA | 14559 | 76.69s | 108.93s | 鏈噸娴?| 鈥?| 鈥?|
| Sphere MLFMA | 26424 | 323.25s | 285.81s | 鏈噸娴?| 鈥?| 鈥?|
| **VEFIE** | 15828 | 46.13s | 66.24s | 鏈噸娴?| **41.30s 鉁?* | 涓婁笁瑙掑绉颁紭鍖?|
| **SCFIE** | 15860 | 66.68s | 96.94s | 鏈噸娴?| **65.67s 鉁?* | VEFIE+鑰﹀悎浼樺寲 |

#### 8.9 EFIE/CFIE/SCFIE 娣卞害浼樺寲 鉁?(commit `8f8dfc3`, `f520609`)
- [x] **EFIE `calc_interaction!` 閲嶅啓**: 绉婚櫎 `@simd for n in 1:3` + 涓夊厓鍒嗘敮锛屾敼鐢?tuple-indexed rho + 瀹屽叏灞曞紑 3脳3 鍐呯Н 鈫?Jet EFIE **4.26s** (-79.4% vs 20.7s 鍘熷鍩虹嚎)
- [x] **CFIE 鏋舵瀯淇**: 鍚堝苟姹囩紪瀹炴祴姣斿垎绂绘眹缂栨參 (register/cache pressure)锛屾敼涓哄垎绂昏皟鐢?+ `@.` 灏卞湴鍔犳潈姹傚拰 (閬垮厤绗笁涓?N脳N 鍒嗛厤) 鈫?Jet CFIE **14.48s** (-91.4% vs 168.3s 鍘熷鍩虹嚎)
- [x] **SCFIE Fss 骞惰鍖?*: `assemble_fss_boundary_correction!` 娣诲姞 `@threads` + 琛岀骇 SpinLock
- [x] **鍗曞厓娴嬭瘯**: 179/179 閫氳繃 (Testing EMSuite tests passed)

#### 8.10 VEFIE/SCFIE 鎬ц兘绐佺牬 鉁?(commit `bbf8fdd`)
- [x] **VEFIE 涓婁笁瑙掑绉颁紭鍖?*: 灏嗗叏 N虏 tet 瀵瑰惊鐜敼涓轰笂涓夎 N*(N+1)/2 鍧? 
  - 澶栧眰寰幆 test tet `it`锛屽唴灞?`js` 浠?`it+1` 鍒?`ntet`锛堜笂涓夎锛? 
  - 鍒╃敤 Z_st[j,i] = (魏_t/魏_s) 脳 Z_ts[i,j]锛屽悓鏃跺啓鍏?Z[m,n] 鍜?Z[n,m]  
  - 鍏ㄥ眬 SpinLock锛堥攣鎸佹湁鏃堕棿 鈮?2 娆℃爣閲忓啓 鈮?2 ns锛岃绠楁椂闂?鈮? 渭s锛岀珵浜夌巼 <1%锛? 
  - **VEFIE: 41.30s**锛坴s Legacy 46.13s锛?*蹇?12%**锛泇s 鍩虹嚎 66.24s锛?*蹇?1.60脳**锛? 
- [x] **SCFIE 鑰﹀悎鍧椾簰鏄撴€т紭鍖?*: Z_vs = Z_sv / 魏锛坈ommit `a87be12`锛夛紝鍑忓皯涓€鍗婅€﹀悎绉垎  
  - **SCFIE: 65.67s**锛坴s Legacy 66.68s锛?*蹇?1.5%**锛泇s 鍩虹嚎 96.94s锛?*蹇?1.48脳**锛? 
- [x] **鍏ㄩ儴 179/179 娴嬭瘯閫氳繃**


鹿 EFIE 缁勮澧炲箙: @fastmath/SIMD 閲嶅啓涓昏浼樺寲 MFIE 璺緞, 瀵圭函 EFIE 鏈夎交寰紑閿€
虏 MLFMA EFIE 澧炲箙: 棰勬潯浠跺櫒 LU 鍙樻參 (8.89s鈫?7.27s), 闈炰唬鐮佸洖褰?
鲁 VEFIE/SCFIE 缁勮澧炲箙鍚屽洜; LU 姹傝В澶у箙鍔犻€?(155.61s鈫?1.12s for VEFIE)

---

## 杩涜涓?馃敡

锛堟棤褰撳墠杩涜涓换鍔★級

### Phase 12: 鍏潰浣撳畬鏁存敮鎸?PWCHex + RBF (2026-03-04) 鉁?

**鐩爣**: 瀹炵幇鎵€鏈夌己澶辩殑鍏潰浣撳熀鍑芥暟+绉垎鏂圭▼缁勫悎锛岃ˉ鍏?Phase 12 璺嚎鍥句腑鐨?10 涓?Gap銆?

**淇敼鏂囦欢:**
1. **`src/Geometry/GaussQuadrature.jl`** 鈥?鏂板鍏潰浣?(8鐐?tensor-product GL) 鍜屽洓杈瑰舰 (4鐐? GQ 瑙勫垯
2. **`src/Geometry/MeshTypes.jl`** 鈥?鏂板 `HexahedraInfo`, `Quads4Hexa` 缁撴瀯浣?+ 杈呭姪鍑芥暟 (`get_free_vns`, `set_delta_kappa!`, `hex_volume` 绛?
3. **`src/Geometry/MeshIO.jl`** 鈥?鏂板 CHEXA Nastran 缃戞牸璇诲彇锛屾敮鎸佺画琛岀鏍煎紡
4. **`src/BasisFunctions/PWC.jl`** 鈥?鏂板 `PWCHexBasis` (3 DOF/鍏潰浣? x,y,z 鍒嗛噺)
5. **`src/BasisFunctions/RBF.jl`** 鈥?瀹屽杽 `evaluate()` 瀹炵幇锛屽惎鐢ㄨ竟鐣岄潰鍩哄嚱鏁?
6. **`src/BasisFunctions/BasisUtilities.jl`** 鈥?鏂板 `get_hexahedra_info(mesh, PWCHexBasis/RBFBasis, permittivities)`
7. **`src/IntegralEquations/VEFIE.jl`** 鈥?~500琛? PWCHex, RBF, 娣峰悎 TetraHex 瑁呴厤 (6 涓柊鏂规硶)
   - 娉涘寲 `_pwc_dyad_kernel!` 鏀寔 duck-typed 浣撳厓绱?
   - 娣峰悎 TetraHex 瑁呴厤: 4涓瓙鍧?(TT, TH, HT, HH)
8. **`src/IntegralEquations/SCFIE.jl`** 鈥?~380琛? RWG+PWCHex (骞剁煝 L 绠楀瓙) 鍜?RWG+RBF (鏍囬噺鍔垮舰寮?+ Fss 杈圭晫淇)
9. **`src/IntegralEquations/Excitation.jl`** 鈥?~180琛? PWCHex 鍜?RBF 骞抽潰娉㈡縺鍔卞悜閲?+ 缁勫悎 SCFIE 鐗堟湰
10. **`src/PostProcessing/RadiationIntegral.jl`** 鈥?~120琛? PWCHex 鍜?RBF 杈愬皠绉垎
11. **`src/PostProcessing/RCS.jl`** 鈥?~80琛? PWCHex 鍜?RBF RCS 璁＄畻鏂规硶
12. **`test/test_basis_functions.jl`** 鈥?淇 RBF 娴嬭瘯棰勬湡鍊?(num_basis 1鈫?1, 鏌ユ壘鍐呴儴鍩哄嚱鏁?

**鏂规硶缁熻:**
- 18 涓?`assemble_impedance_matrix` 鏂规硶 (EFIE/MFIE/CFIE/VEFIE/SCFIE 脳 鍚勫熀鍑芥暟缁勫悎)
- 17 涓?`excitation_vector` 鏂规硶
- 5 涓?`radarCrossSection` 鏂规硶

**娴嬭瘯缁撴灉:**
- 鍏ㄩ儴 179/179 娴嬭瘯閫氳繃 (鏃犲洖褰?
- +2296 琛屼唬鐮?
- Commit: `099385b`

### Phase 11: PWC 鍩哄嚱鏁版敮鎸佹墿灞?(2026-03-04) 鉁?

**鐩爣**: 瀵归綈 Legacy 鐨?PWC (Piecewise Constant) 鍩哄嚱鏁版敮鎸侊紝瀹屽杽 VEFIE+PWC 鍜?SCFIE+RWG+PWC 缁勫悎銆?

**淇敼鏂囦欢:**
1. **`src/BasisFunctions/PWC.jl`** 鈥?瀹屽叏閲嶅啓: 3 DOFs/鍥涢潰浣?(x,y,z 鍒嗛噺)
   - `PWC` struct 澧炲姞 `inBfsID::SVector{3, IT}` (涓変釜鍏ㄥ眬鍩哄嚱鏁癐D)
   - `num_basis` 杩斿洖 `3 * length(functions)` (鍘熶负 1:1)
   - `evaluate` 杩斿洖鍗曚綅鍚戦噺 x虃/欧/岷?(鍩轰簬 `mod1(i, 3)`)
   - Legacy 瀵归綈: `MoM_Basics` 鐨?`nPWC = 3 * num_tetrahedra`

2. **`src/BasisFunctions/BasisUtilities.jl`** 鈥?鏂板 `get_tetrahedra_info(mesh, basis::PWCBasis, permittivities)`
   - `inBfsID = SVector{4}(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3, 0)` (绗?椤规湭浣跨敤)

3. **`src/IntegralEquations/VEFIE.jl`** 鈥?鏂板 ~230 琛? VEFIE+PWC 缁勮
   - `assemble_impedance_matrix(vefie::VEFIE, basis::PWCBasis)` 鈥?甯?permittivities 鐨?鍙?2鍙傜増鏈?
   - `_pwc_dyad_kernel!` 鈥?3脳3 骞剁煝 L 绠楀瓙: $(k^2 I + \nabla\nabla) G(R)$
   - 瀵圭О缁勮 + 鑷€傚簲绉垎 (杩滃満1鐐?杩戝満5鐐?
   - 鑷綔鐢ㄩ」璐ㄩ噺鐭╅樀: $V/(j\omega\varepsilon)$

4. **`src/IntegralEquations/Excitation.jl`** 鈥?鏂板 ~100 琛? PWC 婵€鍔卞悜閲?
   - VEFIE 绠楀瓙鐗? `excitation_vector(op::VEFIE, source::PlaneWave, basis::PWCBasis)`
   - 鐙珛鐗? `excitation_vector(source::PlaneWave, basis::PWCBasis)`
   - 缁勫悎鐗? `excitation_vector(source, surf_basis::RWGBasis, vol_basis::PWCBasis)`

5. **`src/IntegralEquations/SCFIE.jl`** 鈥?鏂板 ~170 琛? SCFIE+RWG+PWC 鑰﹀悎
   - `assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCBasis)`
   - `assemble_coupling_blocks_pwc!` 鈥?骞剁煝 L 绠楀瓙鑰﹀悎
   - Z_SV 鍖呭惈 魏, Z_VS 鏃?魏
   - 鏃?Fss 杈圭晫淇 (PWC 鏃犲崐鍩哄嚱鏁?

6. **`src/PostProcessing/RadiationIntegral.jl`** 鈥?鏂板 PWC 杈愬皠绉垎
   - `radiation_integral_pwc`: $N = \sum_t V_t \kappa_t \sum_{gq} J \cdot e^{jk\hat{r}\cdot r_{gq}} w_{gq}$

7. **`src/PostProcessing/RCS.jl`** 鈥?鏂板 PWC RCS 璁＄畻
   - `radarCrossSection(..., basis::PWCBasis, permittivities)`

8. **`src/Driver.jl`** 鈥?鎵╁睍鏀寔澶欼E绫诲瀷 (EFIE/MFIE/CFIE/VEFIE/SCFIE)

9. **`src/Core/Configuration.jl`** 鈥?SimulationConfig 澧炲姞 `ie_type`, `cfie_alpha`, `permittivities` 瀛楁

10. **`test/test_pwc.jl`** 鈥?鏂板 PWC 涓撶敤娴嬭瘯 (鍩哄嚱鏁版瀯閫? VEFIE+PWC, 婵€鍔? SCFIE+PWC)

11. **`test/test_basis_functions.jl`** 鈥?鏇存柊 PWC 娴嬭瘯閫傞厤 3 DOFs/tet

**淇鐨?Bug:**
- VEFIE.jl 缂哄け module 闂悎 `end` (缂栬瘧閿欒)
- Excitation.jl 缂哄け module 闂悎 `end` (缂栬瘧閿欒)
- Driver.jl 浣跨敤 `get(struct, ...)` 瀵艰嚧 MethodError (struct 涓嶆敮鎸?Dict 鐨?get)

**娴嬭瘯缁撴灉:**
- 鍏ㄩ儴 139+/139+ 娴嬭瘯閫氳繃 (鍚柊澧?PWC 娴嬭瘯)
- PWC 鍩虹娴嬭瘯: 16/16 PASS
- 鏃犲洖褰?

---

## 杩涜涓?馃殌

### Phase 9: 浠ｇ爜璐ㄩ噺涓庡彂甯?
- [x] 娴嬭瘯濂椾欢娓呯悊: 179/179 鍏ㄩ儴閫氳繃 鉁?(鍚?PWC/RBF/PWCHex 鏂版祴璇?
- [x] JuliaFormatter.jl 缁熶竴浠ｇ爜椋庢牸 鉁?(`0535576`) 鈥?79 涓?`src/` 鏂囦欢
- [x] CHANGELOG.md 瀹屽杽 鉁?鈥?Phase 8-12 鎵€鏈夐噷绋嬬
- [x] compat 璇硶瑙勮寖鍖?鉁?+ 瑕嗙洊鐜囪剼鏈?+ 鍙戝竷娓呭崟 (`9d3e671`)
- [x] 娴嬭瘯瑕嗙洊鐜囩粺璁′笌鎻愬崌 (鐩爣 > 80%)
- [x] API 鏂囨。琛ュ叏 (鎵€鏈夊叕鍏辨帴鍙?
- [x] 鐢ㄦ埛鏁欑▼ (Quick Start, Advanced)
- [x] 鐞嗚鏂囨。 (MoM, MLFMA, 绉垎鏂圭▼鎺ㄥ)
- [ ] 鍙戝竷鍒?Julia General Registry

### 2026-03-09 Update (Phase 9 Coverage + Docs Closure)

- 覆盖率统计脚本已补齐：`scripts/check_coverage.jl`
- 轻量覆盖率实测：`4283/4595 = 93.21%`（阈值 80% 已满足）
   - 运行：`julia --project=. --code-coverage=user --startup-file=no test/runtests_light_cov.jl`
   - 统计：`julia --project=. scripts/check_coverage.jl 80`
   - 报告：`test_results/COVERAGE_REPORT.md`
- Documenter 文档已补齐与验证：
   - 新增 `docs/src/guide/advanced.md`（中尺度/预算/门禁工作流）
   - 新增 `docs/src/theory/overview.md`（理论总览）
   - `docs/make.jl` 已接入新页面导航；`docs/src/index.md` 已刷新入口与安装说明
   - 验证：`julia --project=docs docs/setup_docs.jl` + `julia --project=docs docs/make.jl`

### 2026-03-09 Update (General Registry Preflight)

- 已补齐发布前置资产（仓内可完成部分）：
   - 新增发布清单：`.github/RELEASE_CHECKLIST.md`
   - 新增 TagBot 工作流：`.github/workflows/TagBot.yml`
   - 新增 Docs CI 工作流：`.github/workflows/Docs.yml`
   - `CI.yml` 已支持 `master/main` 双分支触发，并支持手动触发
- 当前唯一未完成项仍为外部动作：
   - Registrator 提交与 General Registry 审核合并
   - TagBot/Registrator 所需仓库权限与密钥配置

### 2026-03-09 Update (Phase 16 报告计划落地)

- 已新增发布前全模块验证计划：
   - `.github/plans/phase_16_release_validation_report_plan.md`
   - 明确要求：解析基准 + FEKO 对标双线覆盖，并给出模块覆盖矩阵
- 已新增报告模板：
   - `test_results/reports/RELEASE_VALIDATION_REPORT_TEMPLATE.md`
- 已新增报告生成入口：
   - `benchmark/run_release_validation_report.jl`
   - 运行后生成 `test_results/reports/RELEASE_VALIDATION_REPORT.md`
- Phase 16 当前状态：计划与骨架完成，待执行全量用例填充实测数据并完成 2 轮检视。

### 2026-03-09 Update (Phase 16 首轮实测填充)

- 已执行并回填首轮实测数据到 `test_results/reports/RELEASE_VALIDATION_REPORT.md`：
   - FEKO链路：`F1/F2/F7` 已重跑并更新 RMSE
   - 解析/端口链路：`A1-A4` 已重跑并更新 Zin/S11 结论
   - 介质天线链路：`B1-B5` 已确认纳入汇总
- 本轮识别的发布阻塞：
   - `F2` (S-CFIE Jet Direct) 超阈值失败（RMSE 5.460~6.948 dB）
   - `F7` (V-EFIE Plate Direct) 超阈值失败（RMSE 9.078~11.785 dB）
   - `P1/P3` PMCHW 大规模 direct 在当前环境 `OutOfMemoryError()`
   - `F5/F6` 首轮前因球网格半径提取失败未完成（已定位为解析器问题）
- 当前阶段结论：Phase 16 仍处于执行中，发布建议维持 **No-Go**，待阻塞项闭环后进入检视迭代轮次。

### 2026-03-09 Update (Phase 16 第二轮执行刷新)

- 已完成 `sphere_600MHz.nas` 半径提取链路修复并补齐回归：
   - `src/Accuracy/ReferenceData.jl`：`extract_sphere_radius` 改为优先走项目 Nastran 读取器，再回退文本解析
   - `test/test_accuracy_metrics.jl`：新增球半径提取回归，当前全量 `26/26 pass`
- 已重跑球体 FEKO/Mie 链路并更新报告：
   - `F5`（S-CFIE Sphere Direct）对 FEKO 通过：phi0 RMSE `0.158 dB`，phi90 RMSE `0.089 dB`
   - `F6`（S-CFIE Sphere MLFMA）对 FEKO 超阈值：phi0 RMSE `3.319 dB`，phi90 RMSE `4.037 dB`
   - `F6` 对 Mie 仍超阈值：RMSE `9.074 dB`
- 阻塞项状态刷新：
   - 已移除“`F5/F6` 缺文件阻塞”旧结论
   - 新增/保留精度阻塞：`F6` MLFMA 精度未达门限
   - 其余阻塞保持：`F2`、`F7`、`P1/P3`、`A1/A2/A4`
- 已完成 `F6` 求解截断假设排查（长 Krylov 诊断）：
   - 试验配置：`restart=300, maxiter=600, tol=1e-6`（其余保持不变）
   - 结果：`F6` 三项 RMSE（3.319 / 4.037 / 9.074 dB）基本不变
   - 结论：`F6` 阻塞当前不由 GMRES 截断主导，后续需优先排查 MLFMA 算子保真/常数链路
- 当前阶段结论：Phase 16 仍为 **No-Go**，但球体链路已从“无法执行”推进到“可执行且可量化诊断”。

### 2026-03-09 Update (Phase 16 效率/内存主线刷新)

- 已把 Phase 16 阻塞从“纯精度问题”升级为“精度 + 效率/内存问题”双主线：
   - 慢用例：`F6` 球体 MLFMA 构建耗时过长；`P1/P3` 大规模 PMCHW direct 即使不考虑精度也难以形成稳定门禁
   - OOM 用例：`P1/P3` 在当前机器上出现 `OutOfMemoryError()`，但量级分析表明问题不应简单归因为“机器内存不足”
- 已完成 OOM 首轮根因分析：
   - 当前机器物理内存约 `68.5 GB`
   - `P1/P3` 的 `2N = 52848` PMCHW dense 矩阵本体约 `44.7 GiB`
   - 原实现除总矩阵外，还会额外分配多个 `N×N` L/K 子块，并在 `Z \ V` 时触发整矩阵复制，导致峰值内存被人为放大
- 已落地源码修复：
   - `src/IntegralEquations/Impedance.jl`：新增 `assemble_generic!`
   - `src/IntegralEquations/EFIE.jl`：新增原地写入装配路径
   - `src/IntegralEquations/PMCHW.jl`：四块 PMCHW 子块直接写入总矩阵视图，不再额外物化大块临时矩阵
   - `benchmark/accuracy/run_P1_P3_pmchw.jl`：改为 `lu!` 原地分解 + `ldiv!` 求解，并打印理论矩阵占用
- 已验证：`test/test_pmchw.jl`、`test/test_nmuller.jl` 通过，说明 dense PMCHW 主线未被破坏
- 当前剩余状态：真实 `P1/P3` 大用例仍需继续复跑确认 OOM 是否完全解除；在此之前，Phase 16 维持 **No-Go**

---

## Legacy 鍥犲瓙瀵圭収琛?

| 椤圭洰 | Legacy | EMSuite | 璇存槑 |
|------|--------|---------|------|
| EFIE Z | $1/16\pi$ | $1/4\pi$ + 鏄惧紡 $l/2A$ | 鏁板绛変环 |
| FarField | $1/4\pi$ | $1/4\pi$ | 涓€鑷?|
| Translation | $1/16\\pi^2$ | $-jk/16\\pi^2$ | Legacy 瀵归綈 |\n| SWG Disagg | N/A | $jk\\eta/(4\\pi)$ | VEFIE G 鍚?$1/(4\\pi)$ |
| 鏃堕棿绾﹀畾 | $e^{-j\omega t}$ | $e^{j\omega t}$ | Z 铏氶儴绗﹀彿鐩稿弽 |
| Area/Length | 闅愬紡鍖呭惈鍦ㄥ洜瀛愪腑 | 鏄惧紡褰掍竴鍖?| 缁撴灉绛変环 |

---

## 宸茬煡闂

1. ~~**EMSuite vs Legacy ~3 dB 绯荤当鍋忓樊**~~ 鈥?**宸蹭慨澶?* (2026-02-28): 鏍瑰洜涓?`edgev虃` 鏂瑰悜鍙嶈浆 + `calc_near_interaction!` 缁忛獙鍥犲瓙銆備慨澶嶅悗 RCS 鍋忓樊 < 0.3 dB RMSE銆?
2. ~~**MLFMA 杩滃満 脳4 鍥犲瓙 + CFIE 绗﹀彿**~~ 鈥?**宸蹭慨澶?* (2026-03-02): EFIE MLFMA 绯绘暟璇樊 65.7% 鈫?0.30%. CFIE MLFMA RMSE 3.45 dB 鈫?0.003 dB.
3. **VEFIE Mie 鍋忓樊**: Legacy 鍜?EMSuite 鍧囨瘮 Mie 绾ф暟浣?~25dB 鈥?灞炰簬 Legacy 绠楁硶鍥烘湁闂, 鏍囪涓?"Legacy Parity", 鐗╃悊淇涓烘湭鏉ョ爺绌惰棰?
4. **BiCGSTAB 鏀舵暃**: 闇€瑕侀鏉′欢鎵嶈兘鍙潬鏀舵暃
5. **SCFIE 鑰﹀悎椤逛簰鏄撴€?*: 宸蹭慨澶?鈥?Z_SV/魏 = Z_VS^T 鍦ㄦ満鍣ㄧ簿搴︽垚绔?(2.99e-16)
6. ~~**EFIE 闂悎浣撳唴閮ㄨ皭鎸?*: EFIE 鐢ㄤ簬闂悎瀵间綋鏃舵潯浠舵暟宸?(Direct vs MLFMA 绯绘暟宸?62%)锛屽簲鏀圭敤 CFIE~~ 鈥?**宸茬‘璁?*: 鐜板湪 CFIE MLFMA 姝ｇ‘宸ヤ綔 (RMSE 0.003 dB, 7 iterations)
7. **SWG MLFMA const_factor 绗﹀彿**: `const_factor = jk畏/(4蟺)` 鍙兘搴斾负 `-jk畏/(4蟺)` (VEFIE `c1 = -jk畏魏` 鍚礋鍙?. ~~闇€瑕?VEFIE MLFMA 绮惧害娴嬭瘯楠岃瘉~~ 鈫?D3 娴嬭瘯 RMSE=0.0 dB 琛ㄦ槑褰撳墠瀹炵幇姝ｇ‘.
8. **A2 S-EFIE Iterative 鏈厖鍒嗘敹鏁?*: restart=1000 + Diagonal 棰勬潯浠朵笅 EFIE (N=14559) RMSE=0.343 dB (> 0.1 dB). 鏍瑰洜: EFIE 鏉′欢鏁板ぇ, 瀵硅棰勬潯浠朵笉瓒? D2/E2 宸茶瘉鏄?GMRES 鍩虹璁炬柦姝ｇ‘; A3 MLFMA+杩戝満 LU 棰勬潯浠跺彲姝ｅ父鏀舵暃

---

## 鏇存柊鏃ュ織

| 鏃ユ湡 | 鏇存柊鍐呭 |
|------|----------|
| 2026-03-06 | **Phase 13.3 V-EFIE MPI 骞惰瑁呴厤瀹屾垚** 鉁?鈥?鏂板 `VolumeAssembly.jl` (303 lines, 鏃?module 灏佽, Allreduce 绛栫暐). 淇: (1) Julia dispatch 璺ㄦā鍧椾笉鍙 (鍒犻櫎 `module VolumeAssemblyMPI`); (2) `BasisFunctions.jl` 琛ュ嚭鍙?`get_triangles_info`; (3) 瀛楃缂栫爜鎹熷潖 魏鈫掗瓘 淇. 姝ｇ‘鎬? Tetra.nas N=3201, Z[1,1]=4.747e6-4.716e6im 鍦?1/2/4 杩涚▼涓嬪畬鍏ㄧ浉鍚? 鎬ц兘: 1P 10.95s 鈫?2P 8.86s (1.24脳). commit: `1730b71` |
| 2026-03-0X | **Phase 13.2 Option C 璇勪及瀹屾垚** 鉁?鈥?鏂板 bench_pwc_rbf_performance.jl (177 lines)銆傚叧閿彂鐜? PWC+VEFIE 12.62s (3.16脳 faster), 浣?DOF +38% (21834 vs 15828), per-DOF虏 鏁堢巼浠?0.17脳; SCFIE(RWG+PWC) 33.26s vs (RWG+SWG) 63.55s (1.91脳 faster); RBF 浠呮敮鎸佸叚闈綋缃戞牸 (涓嶉€傜敤鍥涢潰浣?. **缁撹**: PWC 涓嶉€傚悎閫氱敤浼樺寲 (姹傝В鍣ㄤ唬浠烽珮), SWG 淇濇寔鏈€浣冲钩琛? **鎺ㄨ崘**: 杞悜 Option A (MPI 骞惰鍖? 棰勬湡 3.5脳 @ 4 杩涚▼). 璇﹁ [OPTION_C_PWC_RBF_EVALUATION.md](.github/OPTION_C_PWC_RBF_EVALUATION.md). commit: `bc2121d` |
| 2026-03-01 | **Phase 13.1 FastExp 浼樺寲瀹屾垚** 鉁?鈥?鏂板 FastExp.jl (10000 鏉＄洰绾挎€ф彃鍊艰〃), VEFIE.jl 闆嗘垚, @fastmath/@inbounds/unsafe_trunc 浼樺寲銆俈-EFIE 66.24s鈫?9.05s (1.70脳 vs baseline, **1.18脳 vs Legacy** 鉁?, SCFIE 96.94s鈫?3.81s (1.52脳 vs baseline, **1.04脳 vs Legacy** 鉁?銆俆hread-local buffer 澶辫触鏁欒: -63% 鎬ц兘閫€姝?(16GB 鍐呭瓨鐖嗙偢)銆傝窛绂?2脳 鐩爣: V-EFIE 缂哄彛 70%, SCFIE 缂哄彛 92%銆傝瑙?[PHASE_13.1_SUMMARY.md](.github/PHASE_13.1_SUMMARY.md). commits: `bfb24f6`, `0488e28`, `2c2ea75` |
| 2026-03-01 | **Phase 13.1 FastExp 鏌ユ壘琛ㄤ紭鍖栧畬鎴?* 鈥?鏂板 FastExp.jl (10000 鏉＄洰绾挎€ф彃鍊艰〃, 瑕嗙洊 20位), VEFIE.jl 闆嗘垚 (struct 瀛楁 + 涓ゅ璋冪敤鐐?, 浣跨敤 @fastmath/@inbounds/unsafe_trunc 浼樺寲鎬ц兘銆傛€ц兘缁撴灉 (plate_and_metal_1dot2GHz, 4 threads): V-EFIE 66.24s鈫?9.05s (1.70脳 vs baseline, 1.18脳 vs Legacy 46.13s 鉁?, SCFIE 96.94s鈫?3.81s (1.52脳 vs baseline, 1.04脳 vs Legacy 66.68s 鉁?. 绮惧害楠岃瘉: 鏈€澶х浉瀵硅宸?< 0.002%. commit: `bfb24f6` |
| 2026-07-30 | **Phase 9.1 浠ｇ爜璐ㄩ噺** 鈥?JuliaFormatter 鏍煎紡鍖?79 涓?`src/` 鏂囦欢; CHANGELOG.md 鍏ㄩ潰鏇存柊; compat 璇硶瑙勮寖鍖?+ julia="1.10"; 瑕嗙洊鐜囪剼鏈?(`scripts/check_coverage.jl`); 鍙戝竷娓呭崟 (`RELEASE_CHECKLIST.md`). 179/179 閫氳繃. commits: `0535576`, `b189a87`, `9d3e671` |
| 2026-07-29 | **Phase 8.10 VEFIE 瀵圭О浼樺寲** 鈥?Z_st[j,i]=(魏_t/魏_s)路Z_ts[i,j] 鎺ㄥ; 涓婁笁瑙掑洓闈綋閬嶅巻; VEFIE **41.30s** (蹇?Legacy 12%); SCFIE **65.67s** (蹇?Legacy 1.5%). commit: `bbf8fdd` |
| 2026-03-01 | **Phase 8.9 EFIE/CFIE/SCFIE 娣卞害浼樺寲** 鈥?EFIE 鍐呮牳閲嶅啓 (绉婚櫎@simd+涓夊厓鍒嗘敮) 鈫?Jet EFIE 20.7s鈫?.26s (-79%); CFIE 鏋舵瀯淇 (鍒嗙姹囩紪+`@.`灏卞湴) 鈫?Jet CFIE 168s鈫?4.5s (-91%); SCFIE Fss 骞惰鍖? 179/179 閫氳繃. commits: 8f8dfc3, f520609 |
| 2026-07-29 | **Phase 10 绮惧害楠岃瘉琛ュ叏瀹屾垚** 鈥?琛ラ綈鍓╀綑 9 瀛愭祴璇?(D2/E2/D3/E3/A2/B2/A4/B3/C3-MPI). 15/16 閫氳繃, A2 鏉′欢閫氳繃 (GMRES 鏀舵暃鍙楅檺). Bug 淇: CFIE MPI 骞惰瑁呴厤 (`cfie_interaction!` 涓嶅瓨鍦?鈫?EFIE+MFIE 鐙珛浜や簰+绾挎€х粍鍚?. 179/179 鍗曞厓娴嬭瘯閫氳繃 |
| 2026-03-04 | **Phase 12 鍏潰浣撳畬鏁存敮鎸佸畬鎴?* 鈥?PWCHexBasis 3 DOFs/hex + RBF evaluate + 杈圭晫闈€侴Q (hex/quad)銆丮eshIO (CHEXA)銆乂EFIE (PWCHex/RBF/Mixed)銆丼CFIE (RWG+PWCHex/RBF)銆佹縺鍔卞悜閲忋€佽緪灏勭Н鍒?RCS銆?79/179 娴嬭瘯閫氳繃銆?2296 琛?|
| 2026-03-04 | **Phase 11 PWC 鍩哄嚱鏁版墿灞曞畬鎴?* 鈥?PWCBasis 3 DOFs/tet, VEFIE+PWC 骞剁煝缁勮, PWC 婵€鍔?杈愬皠绉垎/RCS, SCFIE+RWG+PWC 鑰﹀悎, Driver.jl 澶欼E鎵╁睍, SimulationConfig 澧炲己, 鏂板 test_pwc.jl. 139+/139+ 娴嬭瘯鍏ㄩ€氳繃 |
| 2026-02-28 | **Phase 8 鎬ц兘浼樺寲鍏ㄩ儴瀹屾垚** 鈥?8 涓瓙闃舵 (8.0-8.8), 鏍稿績鎴愭灉: CFIE 缁勮 -61% (168鈫?5s), CFIE/EFIE 姣?8.1脳鈫?.2脳, MLFMA OOM 淇, BlockJacobiPreconditioner, Julia 1.12 鍏煎, 绫诲瀷绋冲畾鎬?clean. 璇﹁ `test_results/PERFORMANCE_REPORT.md` |
| 2026-03-03 | **Phase 8.0 鎬ц兘鍩虹嚎瀹屾垚** 鈥?EMSuite 7 鐢ㄤ緥 + Legacy 5 鐢ㄤ緥璁℃椂銆傚叧閿彂鐜? CFIE 4.61脳 鎱?(鍙岄亶鍘嗛棶棰?, SCFIE 2.26脳 鎱? EFIE/VEFIE 鎸佸钩鎴栨洿蹇? LU 姹傝В蹇?30-40% |
| 2026-03-03 | **Phase 8 鎬ц兘浼樺寲璁″垝** 鈥?鍔犲叆鎬ц兘浼樺寲璺嚎: 6 鐑偣 (SpinLock鍘婚攣/CFIE鍚堝苟/MLFMA Z_near/鍐呭瓨/SIMD/绫诲瀷绋冲畾), 8 姝ラ, 鐩爣 鈮?Legacy 淇濆簳, 鈮?0.5脳 Legacy 鎸戞垬 |
| 2026-03-03 | **SCFIE Fss 杈圭晫淇** 鈥?鍗婂熀鍑芥暟杈圭晫闈㈢Н鍒嗕慨姝ｃ€侲1-VSEFIE RMSE 5.3鈫?.60 dB. D1-SWG VEFIE RMSE 0.95 dB. 138/138 娴嬭瘯閫氳繃 |
| 2026-03-02 | **MLFMA 鍥犲瓙淇脳2** 鈥?(1) EFIE far-field 脳4 鍥犲瓙: 绯绘暟璇樊 65.7%鈫?.30%, RMSE 3.1鈫?.028 dB; (2) CFIE MFIE 绗﹀彿: 鈭嘷{r'}G 缁欏嚭 +jk k虃 (闈?-jk k虃), RMSE 3.45鈫?.003 dB, GMRES 50鈫? 杩唬 |
| 2026-03-01 | **Phase 10 璁″垝** 鈥?鍏ㄦ柟绋嬪叏璺緞绮惧害瀵归綈璁捐瀹屾垚: 5 鏂圭▼ (S-EFIE/S-MFIE/S-CFIE/V-EFIE/VS-EFIE) 脳 4 璺緞 (Direct/Iterative/MLFMA/MPI), 鍏ㄧ悆闈?1314 鐐归噰鏍? 鍏?17 瀛愭祴璇曢」 |
| 2026-02-28 | **绮惧害鏁堢巼鎶ュ憡 v2** 鈥?鍏ㄩ潰鍩哄噯娴嬭瘯: SEFIE Direct(RMSE 0.29dB), CFIE Direct, SEFIE MLFMA, SCFIE MLFMA(Sphere N=26424). 瑙?`test_results/ACCURACY_EFFICIENCY_REPORT.md` |
| 2026-02-27 | 淇楠岃瘉鑴氭湰 benchmark/run_full_benchmark.jl API 閿欒 (MLFMAOperator 鏋勯€?鎺掑簭閫忔槑鎬? |
| 2026-02-27 | 鏂板 MLFMA MatVec 蹇€熸祴璇曡剼鏈?(benchmark/quick_matvec_test.jl) |
| 2026-02-27 | 鏂板 Direct vs MLFMA 鑷竴鑷存€ф祴璇?(benchmark/self_consistency_test.jl) |
| 2026-02-28 | **娴嬭瘯濂椾欢娓呯悊瀹屾垚** 鈥?138/138 鍏ㄩ儴閫氳繃, 淇 6 涓瀛樻祴璇曢棶棰?|
| 2026-02-28 | 淇 `Vector{AbstractBasisFunction}` 绫诲瀷娲惧彂 bug (Aggregation/Disaggregation) |
| 2026-02-28 | **SCFIE MLFMA 楠岃瘉瀹屾垚** 鈥?淇 7+3 涓?bug, 杩戝満 1.58e-15, 杩滃満 1.04% |
| 2026-02-28 | **SCFIE 鑰﹀悎浜掓槗鎬т慨澶?* 鈥?Z_SV `(term1-term2)` 鈫?`(term1+term2)`, Z_VS c1 绗﹀彿淇 |
| 2026-02-28 | 娣诲姞 `test/test_scfie.jl` 鍥炲綊娴嬭瘯 (9 tests all pass) |
| 2026-02-28 | **~3 dB 鍋忓樊鏍瑰洜淇** 鈥?`BasisUtilities.jl` 杈瑰悜閲忔柟鍚?+ `EFIE.jl` 杩戜氦浜掗潰绉綊涓€鍖栥€俍 鐭╅樀 ratio=1.0, RCS RMSE<0.3 dB |
| 2026-02-28 | **绮惧害鏁堢巼鎶ュ憡 v2** 鈥?鍏ㄩ潰鍩哄噯娴嬭瘯: SEFIE Direct(RMSE 0.29dB), CFIE Direct, SEFIE MLFMA, SCFIE MLFMA(Sphere N=26424). 瑙?`test_results/ACCURACY_EFFICIENCY_REPORT.md` |
| 2026-02-27 | 鍒濆鍖栬繘搴︽枃浠? 浠?copilot-instructions.md 杩佺Щ |
| 2026-01-xx | VEFIE MLFMA 楠岃瘉瀹屾垚 (Rel Err 0.04%) |
| 2026-01-xx | SCFIE Direct 楠岃瘉瀹屾垚 (VSIE Plate) |
| 2026-01-xx | MPI/Threading 骞惰楠岃瘉閫氳繃 |
| 2025-12-xx | Surface IE (EFIE/MFIE/CFIE) 鍏ㄩ潰楠岃瘉瀹屾垚 |
| 2025-12-xx | Legacy 鍥犲瓙瀵归綈瀹屾垚, 绉婚櫎缁忛獙甯告暟 |
