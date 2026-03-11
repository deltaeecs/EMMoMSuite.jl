# EMSuite 重构进度

> 最后更新: 2026-03-11

## 当前总览

- 当前主线问题：MLFMA 多层 upward/downward pass 的逐层对齐与修正
- 当前次主线问题：F7 direct 路径 Round 3 性能回收
- 发布流程、统一报告、依赖瘦身、README 与理论文档修复已完成当前轮次收口

## 2026-03-11 Update：全量文档检视修复

- 已完成理论文档全量扫描，确认问题分为两类：
  - 数学公式被错误写成 `$...$` 的转义形式 `\$...\$` 或 `\$\$...\$\$`，导致渲染失败
  - `.github/REFACTORING_ROADMAP.md` 与 `.github/REFACTORING_PROGRESS.md` 含大面积编码损坏文本
- 已修复 `docs/src/theory/` 下全部受影响理论文档：
  - `basis_functions.md`
  - `electromagnetics.md`
  - `excitations.md`
  - `fast_algorithms.md`
  - `integral_equations.md`
  - `method_of_moments.md`
  - `post_processing.md`
  - `solvers.md`
- 已统一替换为标准 Markdown / Documenter 可渲染数学定界符
- 已重建本进度文档与路线图文档，移除不可维护的乱码历史块，改为“当前快照 + 历史摘要”结构

## 2026-03-11 Update：检视迭代结论（代码 + 文档）

- 已完成一轮针对 MLFMA 多层问题的代码/文档交叉核对。
- 结论一：`.github` 路线图、进度文档与 `docs/src/theory/fast_algorithms.md` 现已一致指向同一事实，即 `nLevels >= 3` 的偏差仍属于 upward/downward pass 实现问题，不能归因于“只是文档没刷新”。
- 结论二：`docs/build/` 仍包含旧构建产物，站点展示内容落后于 `docs/src/` 源文档；后续判断必须以源码文档和测试结果为准，且需要补做一次文档重建。
- 结论三：通用 `test/test_mlfma.jl` 当前仍偏向结构/链路 smoke test，尚未形成针对 `nLevels >= 3` 的通用数值 fidelity 回归，这与路线图中“补齐多层回归测试”的开放项一致。

## 2026-03-11 Update：理论文档按实现约定刷新

- 已将理论文档从“教材版简写”补充为“教材公式 + EMSuite 实现约定”的双层说明，重点覆盖以下高频误读点：
  - RWG 在仓库中采用 `edge_length + signs` 的存储形式，文档现已给出与“正负半基函数”完全等价的统一写法
  - MoM 装配中的 RWG-RWG 配对现已明确写为四个 support 子三角形的求和，并说明符号可吸收到局部几何因子中
  - MLFMA 文档现已覆盖 leaf aggregation、upward pass、translation、downward pass、leaf testing 的完整链路，而不是只有概念性流程图
  - MLFMA 文档现已明确 `sorted_ids` 只服务于八叉树内部遍历；外部 `mul!` 的输入/输出仍是物理 basis 顺序，避免 benchmark 误用
  - Solver 文档现已明确当前 GMRES/左预条件写法对应的是 `M^{-1} Z I = M^{-1} V`，并建议同时检查预条件残差与物理残差
- 本轮刷新文件：
  - `docs/src/theory/basis_functions.md`
  - `docs/src/theory/method_of_moments.md`
  - `docs/src/theory/fast_algorithms.md`
  - `docs/src/theory/solvers.md`
- 当前结论：关于“RWG 正负号如何统一写进公式”和“MLFMA / GMRES 在 EMSuite 中的实际使用口径”这两个知识缺口，文档侧已补齐。

## 2026-03-11 Update：公式渲染检视 Round 1

- 已按“公式内等号不要单独一行”的规则追加一轮理论文档渲染检视。
- 本轮定位到两个明确风险点：
  - `docs/src/theory/fast_algorithms.md` 中 3 组 MLFMA upward / translation / downward 公式把 `=` 单独放在一行
  - `docs/src/theory/method_of_moments.md` 中面-体耦合块矩阵公式把 `=` 单独放在一行
- 已全部改为单行等号写法，避免 Documenter / Markdown 数学渲染器把公式切断。
- 当前结论：本轮新增或近期重写的理论公式中，已确认的“独立等号”渲染风险已清除；文档构建剩余告警当前集中在 `electromagnetics.md`，并且 `fast_algorithms.md` 中仍存在未转义美元符号，`integral_equations.md` 的旧告警仍需在后续轮次复查。

## 2026-03-11 Update：README 刷新

- 已将 `README.md` 从“历史阶段总结页”调整为“当前可执行使用说明页”
- 已重绘架构总览，按 runtime / benchmark / release 三层组织
- 已更新快速开始、使用指引、模块结构、性能基线、精度验证、结果可视化等章节
- 已将 `benchmark/run_release_validation_report.jl` 升级为支持从根环境自动切换到 benchmark 环境

## 2026-03-11 Update：Phase 17 发布流程规范化与依赖瘦身收口

### 交付物

- 新增 `src/Accuracy/ReleaseWorkflow.jl`
- 冻结 `benchmark/configs/*.toml` 配置体系
- 新增统一入口 `benchmark/run_release_workflow.jl`
- 将 release report 拆分为 `collector.jl` / `plotting.jl` / `writer.jl`
- 引入 case status matrix、known exception registry、artifact index、run manifest
- 极坐标远场图纳入统一报告

### 依赖收敛

- `Plots` 已隔离到 `benchmark/Project.toml`
- `MoM_Visualizing` 已从主依赖移除
- 已以仓库内轻量实现替换：
  - `LoggingExtras`
  - `Primes`
  - `Roots`
  - `NearestNeighbors`
  - `ProgressMeter`

### 验证结果

- `using EMSuite` 通过
- `test/test_release_workflow.jl` 通过
- `test/test_benchmark_report_data.jl` 通过
- `benchmark/run_release_validation_report.jl` 通过
- `benchmark/run_release_workflow.jl benchmark/configs/release_quick.toml` 通过

## 2026-03-11 Update：Phase 16 统一图表化报告收口

- 统一报告入口固定为 `benchmark/run_release_validation_report.jl`
- 已形成远场曲线图、性能图、Markdown 总报告的统一产物链
- 已完成多轮检视，确认 Phase 16 以 known legacy exception 方式收口

## 2026-03-10 Update：F7 Legacy 对齐与 direct 性能回收

### Legacy 对齐

- 建立 F7 点对点对比脚本 `scripts/verification/debug_f7_legacy_pair.jl`
- 确认几何、SWG 计数、激励、等效体电流、RCS 后处理逐层一致
- near/far 门控切换为 Legacy 一致公式 `0.15 * lambda0 / sqrt(abs(eps_r))`
- 修复 tetra 11 点 singular quadrature 坐标列错排
- 切换 far 主路径到 Legacy 公式 `_ordered_swg_far_kernel(5,4)`
- F7 主用例对齐结果：
  - `phi0 RMSE = 0.097 dB`
  - `phi90 RMSE = 0.126 dB`

### 性能回收

- Round 1: `1830.4 s -> 432.0 s`
- Round 2: `432.0 s -> 82.2 s`
- Round 3: 仍在进行，目标是优于 Legacy 两倍

## MLFMA 算法分析阶段

- 已完成 Legacy MLFMA 算法报告：`.github/plans/legacy_mlfma_algorithm_report.md`
- 已确认 EMSuite 与 Legacy 的总系数链一致：`k^2 eta / 16pi^2`
- 已定位多层误差主源位于 upward/downward pass，而非整体系数、聚合或近场门控
- 已完成两轮检视，当前进入逐层节点级别排查阶段

## 历史阶段摘要

### Phase 1-6

- 完成项目骨架、核心模块、几何 / 基函数 / 方程 / 快速算法 / 求解器 / 后处理主干落地

### Phase 7-10

- 完成 Legacy 对齐、关键误差修正、性能主线优化与精度收敛

### Phase 11-12

- 完成 PWC、PWCHex、RBF 及混合体元路径支持

### Phase 13-17

- 完成 MPI / 性能探索、统一验证报告、发布流程规范化、依赖瘦身与文档刷新

## 当前未完成事项

- [ ] MLFMA 多层 upward/downward pass 修正
- [ ] F7 direct Round 3 性能回收
- [ ] Julia General Registry 发布