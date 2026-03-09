# Phase 16 发布前全模块效率与精度验证计划

> 创建日期: 2026-03-09  
> 状态: 计划已冻结，执行中  
> 关联目标: 发布前形成一份覆盖所有核心模块的统一验证报告（效率 + 精度）

---

## 0. 开发原则遵循声明

本 Phase 适用并遵循以下原则（来自 `.github/copilot-instructions.md`）：

- 原则 1（TDD 工作流）：新增报告聚合逻辑先定义输出结构，再补实现与回归。
- 原则 2（严格 Legacy 对齐）：Legacy 相关对标仅以 `MoM_Kernels` / `MoM_Basics` / `MoM_AllinOne` 为真值源。
- 原则 3（差异排查流程）：若与 Legacy/FEKO/解析解不一致，按几何-常数-积分-组装顺序定位。
- 原则 4（进度同步）：本计划执行后必须同步 `REFACTORING_ROADMAP.md` 与 `REFACTORING_PROGRESS.md`。
- 原则 5（提交规范）：每个逻辑闭环单独提交（脚本、报告、文档分开）。
- 原则 7（检视迭代）：报告完成后做至少 2 轮 clean review。
- 原则 9（计划文档规范）：本计划明确基准、DoD、检视与回链。

---

## 1. Legacy / 参考基准锚点

### 1.1 Legacy 基准锚点

- `MoM_AllinOne/meshfiles/`：Jet、Sphere、TriTetra 等历史主夹具。
- `test_results/legacy_baseline/`：已沉淀的 Legacy 对照数据（CSV）。
- `benchmark/verify_*` 与 `test/test_legacy_parity.jl`：关键常数与链路对齐锚点。

### 1.2 FEKO 基准锚点

- `C:/Users/12253/OneDrive/MoM/MoM_AllinOne/deps/compare_feko/`
- 最低必测：
  - `jet_100MHzRCS.csv`
  - `sphere_600MHzRCS.csv`
  - `plate_1dot2GHzRCS.csv`
  - `plate_metal_1dot2GHzRCS.csv`

### 1.3 解析基准锚点

- Mie：`Utilities/MieSeries.jl`（PEC 球/介质球）。
- 可解析天线门：半波偶极子输入阻抗与方向图参考（已有 A1-A4 框架）。

---

## 2. 报告目标与覆盖边界

输出一份统一发布验证报告，包含：

1. 模块效率：组装、求解、总耗时、内存行为（可得时）。
2. 模块精度：对 Legacy / FEKO / 解析基准的 RMSE、关键角误差、阻抗误差。
3. 功能覆盖矩阵：模块是否被用例命中与结果状态。

### 2.1 模块覆盖要求

- Geometry / Mesh I/O
- BasisFunctions（RWG, SWG, PWC, PWCHex, RBF）
- IntegralEquations（EFIE, MFIE, CFIE, VEFIE, SCFIE, PMCHW, NMuller）
- FastAlgorithms（MLFMA, Lebedev）
- Solvers（LU, GMRES, BiCGSTAB, Preconditioners）
- Ports / Excitations（PlaneWave, DeltaGap）
- PostProcessing（RCS, Far/Near field, input impedance）
- Parallel（MPI/Threading 路径，至少专门门禁验证）
- IO（CSV/HDF5/VTK 输出链）

---

## 3. 用例矩阵（含解析与 FEKO）

| 组别 | 用例ID | 方程/路径 | 对标源 | 指标 |
|------|--------|-----------|--------|------|
| 解析 | M1 | PEC Sphere EFIE Direct | Mie | RMSE(dB), backscatter |
| 解析 | M2 | Dielectric Sphere PMCHW Direct | Mie | RMSE(dB), 关键角差 |
| FEKO | F1 | Jet 100MHz SEFIE Direct | FEKO CSV | RMSE(dB), Max|Diff| |
| FEKO | F2 | Sphere 600MHz SCFIE/CFIE | FEKO CSV | RMSE(dB) |
| FEKO | F3 | Plate 1.2GHz VEFIE | FEKO CSV | RMSE(dB) |
| FEKO | F4 | PlateMetal 1.2GHz SCFIE | FEKO CSV | RMSE(dB) |
| 端口 | A1-A4 | Dipole 系列 | 解析/参考阻抗 | Re/Im 误差, Dmax |
| 介质端口 | B1-B5 | PMCHW/VS-* 天线链 | Direct/脚本参考 | Zin gap, pass/fail |
| 形成对照 | N1 | PMCHW vs NMuller dense | 内部真值（LU） | cond 比、GMRES轨迹 |
| 后端保真 | G1 | PMCHW MLFMA budget/krylov | Dense shell | rel_I, Zin gap |

---

## 4. 执行入口与产物

### 4.1 执行入口（现有）

- `benchmark/run_full_benchmark.jl`
- `benchmark/run_full_accuracy_benchmark.jl`
- `benchmark/accuracy/run_B1_B5_antenna.jl`
- `benchmark/compare_pmchw_mlfma_budget.jl`
- `benchmark/compare_pmchw_mlfma_budget_krylov.jl`
- `benchmark/compare_pmchw_nmuller_sphere.jl`

### 4.2 专门回归入口（现有）

- `test/test_pmchw_gate_s_mlfma_medium.jl`
- `test/test_pmchw_mlfma_budget_medium.jl`
- `test/test_pmchw_mlfma_budget_krylov_medium.jl`
- `test/test_pmchw_block_fidelity_medium.jl`
- `test/test_nmuller_comparison_medium.jl`

### 4.3 报告产物（本 Phase）

- 主报告：`test_results/reports/RELEASE_VALIDATION_REPORT.md`
- 附录（CSV）：继续复用 `test_results/accuracy/*.csv`
- 覆盖率报告：`test_results/COVERAGE_REPORT.md`

---

## 5. 报告结构规范

报告统一分为：

1. 环境信息（Julia、线程、日期、git commit）
2. 执行摘要（通过率、关键风险）
3. 效率结果表（按模块/用例）
4. 精度结果表（解析、FEKO、Legacy）
5. 功能覆盖矩阵
6. 风险清单与发布建议

---

## 6. DoD（完成定义）

满足以下条件才可勾选完成：

- 报告文件已生成：`test_results/reports/RELEASE_VALIDATION_REPORT.md`
- 报告中同时包含：
  - 至少 2 组解析对标结果（Mie/解析阻抗）
  - 至少 4 组 FEKO 对标结果
  - 全模块覆盖矩阵（命中状态）
- 关键门限（建议值，可按历史阈值微调）：
  - Direct vs FEKO/Mie：RMSE <= 2.5 dB
  - MLFMA vs Dense/参考：RMSE <= 3.0 dB 或已解释的 budget 边界
  - 端口误差：Re < 5%，Im < 20 ohm（按用例定义）
- 检视迭代 >= 2 轮且无新增阻塞问题

---

## 7. 检视迭代计划

### Round 1（完整性）

- 检查所有用例是否执行并进入报告。
- 检查模块覆盖矩阵是否存在漏项。

### Round 2（可信性）

- 抽样复核 FEKO/解析对标原始 CSV 与报告指标计算一致性。
- 复核失败项是否有明确归因（算法/预算/收敛）。

连续两轮无新增阻塞问题后，进入发布动作。

---

## 8. 回链

- Roadmap 勾选位置：`Phase 9` 最后一项（发布到 General Registry）前置质量门
- Progress 更新点：新增“Phase 16 发布验证计划与执行日志”
- 模板文件：`test_results/reports/RELEASE_VALIDATION_REPORT_TEMPLATE.md`

---

## 9. 当前执行快照（2026-03-09 第二轮）

- 已完成：
  - `F5/F6` 球体链路可执行化（半径提取修复后已可稳定运行）
  - `F5` Direct 对 FEKO 通过（phi0/phi90 RMSE 均低于 2 dB）
  - `test_results/reports/RELEASE_VALIDATION_REPORT.md` 与 `test_results/accuracy/ACCURACY_REPORT.md` 已回填第二轮数据
- 当前阻塞（No-Go 保持）：
  - `F2`（S-CFIE Jet Direct）FEKO 超阈值
  - `F6`（S-CFIE Sphere MLFMA）FEKO/Mie 超阈值
  - `F7`（V-EFIE Plate Direct）FEKO 超阈值
  - `P1/P3` PMCHW direct OOM
  - `A1/A2/A4` 解析天线门失败
- 说明：
  - `F5/F6` “缺失 sphere_600MHz.nas”已确认为过时结论，真实根因是半径提取实现问题，现已修复并有测试覆盖。
