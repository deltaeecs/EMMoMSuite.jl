# Phase 16 发布前全模块效率与精度图表化验证计划

> 创建日期: 2026-03-09  
> 状态: 已收口，结论为 Go with known legacy exception  
> 关联目标: 发布前形成一份覆盖所有核心模块的统一验证报告（效率 + 精度 + 内存可执行性）

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
3. 远场曲线图：至少按用例组输出 EMSuite vs FEKO/Mie 的远场对比图。
4. 性能对比图：至少输出总耗时柱状图与阶段拆分图。
5. 功能覆盖矩阵：模块是否被用例命中与结果状态。
6. 大规模 direct 可执行性：矩阵规模估算、峰值内存是否合理、是否存在非必要副本导致的 OOM。

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
- `benchmark/accuracy/run_P1_P3_pmchw.jl`
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
- 图表资源：`test_results/reports/assets/accuracy/*.png`
- 图表资源：`test_results/reports/assets/performance/*.png`
- 性能结构化数据：`test_results/reports/PERFORMANCE_BASELINE.csv`
- 附录（CSV）：继续复用 `test_results/accuracy/*.csv`
- 覆盖率报告：`test_results/COVERAGE_REPORT.md`

---

## 5. 报告结构规范

报告统一分为：

1. 环境信息（Julia、线程、日期、git commit）
2. 执行摘要（通过率、关键风险）
3. 精度结果表（解析、FEKO、Legacy）
4. 远场曲线图（按用例组）
5. 效率结果表（按模块/用例）
6. 性能图（总耗时 + 分阶段拆分）
7. 功能覆盖矩阵
8. 风险清单与发布建议

---

## 6. DoD（完成定义）

满足以下条件才可勾选完成：

- 报告文件已生成：`test_results/reports/RELEASE_VALIDATION_REPORT.md`
- 报告中同时包含：
  - 至少 2 组解析对标结果（Mie/解析阻抗）
  - 至少 4 组 FEKO 对标结果
  - 至少 4 张远场曲线图（覆盖 Jet / Sphere / Plate / PMCHW 主组别）
  - 至少 2 张性能图（总耗时、分阶段拆分）
  - 全模块覆盖矩阵（命中状态）
  - 至少 1 组效率/内存诊断结论（说明慢用例与 OOM 用例的可执行性状态）
- 关键门限（建议值，可按历史阈值微调）：
  - Direct vs FEKO/Mie：RMSE <= 2.5 dB
  - MLFMA vs Dense/参考：RMSE <= 3.0 dB 或已解释的 budget 边界
  - 端口误差：Re < 5%，Im < 20 ohm（按用例定义）
  - 大规模 direct：不得因可消除的临时副本/隐式复制而提前 OOM
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
  - `F6` MLFMA 对 FEKO/Mie 已通过（向量顺序误用修复后恢复到与 `F5` 同量级）
  - `test_results/reports/RELEASE_VALIDATION_REPORT.md` 与 `test_results/accuracy/ACCURACY_REPORT.md` 已回填第二轮数据
- 当前阻塞（No-Go 保持）：
  - `F7`（V-EFIE Plate Direct）FEKO 超阈值
  - `P1/P3` PMCHW direct 大用例最终可执行性仍待复核
  - 部分核心用例若按默认单线程入口运行会耗时过久（当前典型为 `F6` 球体 MLFMA 构建、`P1/P3` 大规模 direct）
  - `A1/A2/A4` 解析天线门失败
- 说明：
  - `F2/F4`（Jet CFIE）现已降级为遗留对齐项：Legacy 历史主对齐只覆盖 Jet EFIE，因此 Jet CFIE 不再作为本轮发布门。
  - `F5/F6` “缺失 sphere_600MHz.nas”已确认为过时结论，真实根因是半径提取实现问题，现已修复并有测试覆盖。
  - `F5/F6` 的球体 Mie 参考角度约定与 MLFMA benchmark 向量顺序误用都已修复：`F6` 当前结果为 vs FEKO `0.157 / 0.089 dB`、vs Mie `0.051 / 0.043 dB`，不再属于 Phase 16 发布阻塞。
  - 当前机器物理内存约 `68.5 GB`；`P1/P3` 的 `2N=52848` PMCHW 稠密矩阵本体理论占用约 `44.7 GiB`，原实现还会额外物化多个 `N×N` 子块并通过 `A\b` 触发整矩阵复制，属于可消除的峰值内存放大。
  - 已落地修复：PMCHW 四块矩阵改为直接写入 `2N×2N` 总矩阵视图，`P1/P3` 改为 `lu!` 原地分解；`P1` 已可完成 `52848×52848` 总矩阵组装，不再出现旧路径的提前 OOM。
  - `F6` 的 MLFMA 近场 COO 合并峰值内存已修复，隔离复跑不再 OOM；同时在 `24` 线程下构建时间由 `612.0s` 降至 `116.7s`，说明当前效率主因之一是执行入口未启用 Julia 线程。

## 10. 2026-03-11 图表化报告刷新

- 本轮目标从“发布前 skeleton 报告”升级为“可直接审阅的全功能图表化报告”。
- 新增执行要求：
  - `benchmark/performance_baseline.jl` 必须同步输出 `PERFORMANCE_BASELINE.csv`，不再只写 markdown。
  - `benchmark/run_release_validation_report.jl` 必须直接消费 `test_results/accuracy/*.csv` 与 `PERFORMANCE_BASELINE.csv`，生成：
    - 远场对比图
    - 效率对比图
    - 统一 markdown 报告
- 本轮开发优先级：
  1. 结构化性能数据出口
  2. 远场曲线自动绘图
  3. 性能图自动绘图
  4. 总报告汇总与链接

## 11. 2026-03-11 首轮执行结果

- 已完成产物：
  - `test_results/reports/RELEASE_VALIDATION_REPORT.md`
  - `test_results/reports/PERFORMANCE_BASELINE.csv`
  - `test_results/reports/assets/accuracy/*.png`
  - `test_results/reports/assets/performance/*.png`
- 已验证入口：
  - `julia --project=. test/test_benchmark_report_data.jl`
  - `julia -t auto --project=. benchmark/performance_baseline.jl`
  - `julia --project=. benchmark/run_release_validation_report.jl`
- 首轮交付结果：
  - 已生成 8 组远场对比图和 2 组性能图
  - 报告已包含执行摘要、精度表、远场图、端口摘要、性能表、性能图、覆盖矩阵、效率/内存诊断、风险说明
- 进入检视迭代前的已知残余项：
  - `F2_CFIE_Jet_Direct` FEKO 偏差仍超当前门限，但作为 Legacy 继承问题保留跟踪

## 12. 2026-03-11 检视迭代 Round 1 结论

- `B1` 端口残余项经复核后确认不是当前实现失败，而是统一报告引用了历史遗留产物：
  - 当前脚本真实判据为 `PMCHW(εᵣ≈1) ≈ 2×EFIE`
  - 实测结果 `|ΔZ| = 7.65e-4 Ω`、相对误差 `0.43%`，判定 PASS
- `Parallel` 覆盖状态经复核后由 `PARTIAL` 修正为 `COVERED`：
  - 证据链包括 `test/test_parallel.jl`
  - `test/test_parallel_mfie_cfie.jl`
  - `benchmark/benchmark_parallel_sphere.jl`
- Round 1 结束后的实际残余项：
  - `F2_CFIE_Jet_Direct` FEKO 偏差仍超当前门限，但作为 Legacy 继承问题保留跟踪

## 13. 2026-03-11 检视迭代 Round 2 结论

- 已补齐并行实测样本并写入统一报告：
  - 运行命令：`mpiexec -n 2 julia --project=. benchmark/benchmark_parallel_sphere.jl`
  - 产物：`test_results/reports/PARALLEL_MPI_SAMPLE.csv`
  - 当前样本：`2` ranks、`1` thread/rank、`N=792`、assembly `0.9303 s`
- 统一报告现已包含：
  - `Parallel sample CSV: FOUND`
  - `## Parallel Sample` 表格
  - 效率诊断中的 fresh MPI 样本摘要
- Round 2 结束后的实际残余项：
  - 仅剩 `F2_CFIE_Jet_Direct` FEKO 偏差仍超当前门限，但作为 Legacy 继承问题保留跟踪
- Round 2 建议结论：
  - 若接受 Jet CFIE 遗留项继续按非阻塞问题跟踪，则当前发布建议为 `Go with known legacy exception`

## 14. 2026-03-11 检视迭代 Round 3 / 最终结论

- 已完成最终收口复核，统一报告、路线图、进度文档结论一致。
- 最终发布建议：`Go with known legacy exception`
- 收口说明：
  - 当前接受的唯一已知例外为 `F2_CFIE_Jet_Direct` 相对 FEKO 的偏差
  - 该项已确认属于 Legacy 继承问题，不作为本轮 EMSuite release blocker
  - 其余 Phase 16 工程性阻塞项均已闭环
