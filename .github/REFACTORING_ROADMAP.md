# EMSuite 重构路线图

> 最后更新: 2026-03-11（PMCHW medium 逐层定位已确认主因是 leaf source aggregation 的 3 点积分不足；source / receive 均切到 4 点后 medium gate 显著收敛）

## 当前焦点

### 1. MLFMA 主链分析与修正

- [x] 完成 Legacy MLFMA 算法拆解与系数链核对
- [x] 确认 `nLevels >= 3` 的偏差主要位于 PMCHW `M-pass` 的 leaf source / receive 链路，而非 `aggregate_upward!`
- [x] 针对 PMCHW medium `M-pass` 做逐层点对点比对，定位到首个主导失配节点
- 当前定位进展：`benchmark/compare_pmchw_upward_downward_localization_medium.jl` 已先用 exact reintegration 对 `upward` 各层父盒做点对点对照，证明旧 3 点 source aggregation 的父盒误差虽仅约 `1e-4 ~ 1e-3`，但把 exact-upward 改为 4 点 source 积分后，default `k0/k1` 的最差叶盒误差可从 `5.09e-3 / 2.33e-2` 直接降到 `3.90e-5 / 2.70e-5`，loose `k0/k1` 也分别降到 `6.69e-3 / 4.37e-4`。据此已将 `aggregate_leaf_pmchw!` 从 3 点提升到 4 点，并与既有 4 点 receive 配套；`test/test_pmchw_block_fidelity_medium.jl` 当前已收敛到 `default_m_k0 ≈ 4.02e-6`、`default_m_k1 ≈ 1.28e-6`、`loose_m_k0 ≈ 2.32e-4`、`loose_m_k1 ≈ 3.93e-5`
- [x] 在修正后补齐回归测试，覆盖 `nLevels >= 3` 多层场景（已新增 `test/test_pmchw_multilevel_quadrature_regression.jl` 与 `test/runtests_batch6.jl`）
- [x] 完成两轮检视迭代，确认无新问题后再推进下一阶段

### 2. F7 Legacy 对齐后的性能回收

- [x] F7 主用例完成 Legacy 对齐，远场误差已收敛到 dB 子量级
- [x] direct 性能回收 Round 1: `1830.4 s -> 432.0 s`
- [x] direct 性能回收 Round 2: `432.0 s -> 82.2 s`
- [ ] direct 性能回收 Round 3：继续面向“优于 Legacy 两倍”目标排查算子热点与缓存策略

### 3. 文档与发布链维护

- [x] README 已刷新为当前可执行的使用说明
- [x] 理论文档源码已恢复为标准 Markdown 数学写法 `$...$` / `$$...$$`
- [x] `docs/make.jl` 已增加构建前数学定界符转换层，兼容 Documenter / Julia Markdown
- [x] 本地文档站点已切换为 `.html` 直达链接，避免 `file` 浏览下二次点击目录页
- [x] `.github/REFACTORING_ROADMAP.md` / `.github/REFACTORING_PROGRESS.md` 已从损坏文本重建为可维护版本
- [x] 将理论文档补充为“教材公式 + EMSuite 实现约定”的双层说明，统一记录 RWG 符号折叠、MLFMA 向量顺序与预条件残差口径
- [ ] 后续若新增 Phase，计划文档需继续遵循 `copilot-instructions.md` 中的计划规范与检视要求

## 当前状态快照

### 已完成的核心里程碑

- [x] Phase 1-6：基础架构、核心模块、几何 / 基函数 / 积分方程 / 快速算法 / 求解器 / 后处理主干完成
- [x] Phase 7：与 Legacy 的关键算子与结果对齐完成
- [x] Phase 7.5：约 3 dB 系统偏差根因修复完成
- [x] Phase 8：主要性能优化阶段完成，关键 direct / MLFMA 路径已显著收敛
- [x] Phase 9：代码质量、格式化、测试覆盖、文档与发布准备主线完成
- [x] Phase 10：全路径精度因子修复与验证闭环完成
- [x] Phase 11：PWC 四面体基函数支持完成
- [x] Phase 12：PWCHex / RBF / TetraHex 混合路径支持完成
- [x] Phase 16：统一精度 / 性能图表化报告完成
- [x] Phase 17：发布流程规范化、统一入口、报告拆层、依赖瘦身完成

### 当前开放项

- [ ] General Registry 发布仍未执行
- [ ] MLFMA 多层问题已在 PMCHW medium `M-pass` 上完成主因修正；当前代码侧回归与两轮检视已收口，剩余事项主要是 `docs/build` 站点重建与后续更大 case 覆盖
- [ ] F7 direct Round 3 性能回收仍未收口

## 近期已完成工作

### 2026-03-10 F7 Legacy 对齐排查

- [x] 建立 EMSuite vs Legacy 的 F7 点对点对比脚本 `scripts/verification/debug_f7_legacy_pair.jl`
- [x] 确认几何、SWG 基函数计数、激励向量、等效体电流、RCS 后处理逐层与 Legacy 对齐
- [x] 将 SWG tetra VEFIE 主路径 near/far 门控改为 Legacy 一致公式 `0.15 * lambda0 / sqrt(abs(eps_r))`
- [x] 新增回归测试 `test/test_legacy_parity.jl` 锁定该阈值公式
- [x] 修复 tetra 11 点 singular quadrature 坐标列错排
- [x] 切换 off-diagonal far 主路径到 Legacy 公式 `_ordered_swg_far_kernel(5,4)`
- [x] F7 主用例完成对齐：`phi0 RMSE = 0.097 dB`, `phi90 RMSE = 0.126 dB`

### 2026-03-11 Phase 16 统一报告收口

- [x] 交付统一远场图、性能图、总报告入口 `benchmark/run_release_validation_report.jl`
- [x] 将结构化 CSV 明确为后续性能绘图唯一数据源
- [x] 完成多轮检视迭代，确认以已知 Legacy 遗留项收口 Phase 16

### 2026-03-11 Phase 17 发布流程规范化与依赖瘦身

- [x] 新增 `src/Accuracy/ReleaseWorkflow.jl`
- [x] 固化 `benchmark/configs/*.toml` 配置体系
- [x] 建立 `benchmark/run_release_workflow.jl` 统一编排入口
- [x] 将报告链拆分为 `collector.jl` / `plotting.jl` / `writer.jl`
- [x] 极坐标远场图纳入正式交付
- [x] `Plots` 隔离到 `benchmark/Project.toml`
- [x] 用仓库内轻量实现替换 `LoggingExtras`、`Primes`、`Roots`、`NearestNeighbors`、`ProgressMeter`
- [x] 完成 `using EMSuite`、release workflow、report data 等验证闭环

### 2026-03-11 文档修复

- [x] 重写 README，移除过时的阶段描述、脚本入口和安装说明
- [x] 修复 `docs/src/theory/*.md` 中被错误转义的数学定界符
- [x] 重建本路线图文档，清理大面积编码损坏文本
- [x] 重建进度文档，保留当前事实并压缩历史摘要

### 2026-03-11 检视迭代补记

- [x] 确认源码文档与路线图已同步到“MLFMA 多层问题仍在分析中”的当前状态
- [x] 确认 `docs/build/` 仍停留在旧构建产物，不能作为判断 MLFMA 现状的事实来源
- [x] 确认通用 `test/test_mlfma.jl` 目前以 smoke test 为主，尚未形成 `nLevels >= 3` 的通用数值回归闸门
- [x] 已补入 `test/test_pmchw_multilevel_quadrature_regression.jl`，以 medium PMCHW `M-pass` far-only 对照直接锁定 `nLevels > 3` 多层 source/receive quadrature 一致性
- [x] 检视 Round 1 已发现并修正 `benchmark/compare_pmchw_upward_downward_localization_medium.jl` 的参考语义漂移：主 `exact-upward` 现对齐当前 4 点 source 实现，旧 3 点链路明确标记为 legacy 对照
- [x] 检视 Round 2 已完成 clean review：未发现新的功能正确性 / 架构问题；仅记录 `_receive_terms` 中未使用 `η` / `normal` 的代码整洁度观察项，不影响本轮放行
- [x] 检视 Round 3：落地修复 `_receive_terms` 死代码（移除 `η` 参数与 `normal`/`v1`/`v2`/`v3`），并发现并修正 benchmark 脚本 3 处旧接口调用（旧签名传入 `eta` 会在运行时抛 `MethodError`）；Round 3 **不属于 clean round**
- [x] 检视 Round 4（clean）：确认所有接口签名一致、相位符号正确、batch 入口结构合理————第 1 个连续 clean round
- [x] 检视 Round 5：发现并删除 Phase 17 添加但从未接入的 `src/Accuracy/ReleaseWorkflow.jl`（功能已由 `benchmark/support/release_support.jl` 覆盖）；Round 5 **不属于 clean round**，继续 Round 6
- [x] 检视 Round 6（clean）：Phase 17 依赖收敛变更链一致，连续 1
- [x] 检视 Round 7（clean）：ProgressMeter/Roots/runtests 等变更一致，连续 2
- [x] 检视 Round 8（clean）：`factor_values`/`knn_bruteforce`/`find_zero_bisection` 语义等价验证通过，连续 3
- [x] **连续 3 轮检视无新问题，终止条件满足，PMCHW multilevel M-pass 检视迭代收口**
- [ ] 补做 `docs/build/` 站点重建，消除源码文档与本地构建产物之间的状态滞后

### 2026-03-11 理论文档实现约定补记

- [x] 在 `basis_functions.md` 中补齐 RWG 的统一符号写法，明确仓库采用“公共边长 + support-local sign”等价表示
- [x] 在 `method_of_moments.md` 中补齐 RWG-RWG 四子三角形配对与符号并入系数后的装配公式
- [x] 在 `fast_algorithms.md` 中补齐 leaf aggregation / upward pass / translation / downward pass / leaf testing 的完整链路说明
- [x] 在 `fast_algorithms.md` 中明确 `sorted_ids` 仅用于 MLFMA 内部遍历，外部输入输出仍以物理基函数顺序为准
- [x] 在 `solvers.md` 中补齐左预条件 GMRES 的数学形式，并记录“预条件残差 + 物理残差”双口径检查建议

### 2026-03-11 文档渲染检视补记

- [x] 针对理论文档追加“公式内等号不得单独成行”的渲染约束检查
- [x] 修正 `fast_algorithms.md` 中 upward / downward pass 三条公式的独立等号写法
- [x] 修正 `method_of_moments.md` 中面-体耦合块矩阵公式的独立等号写法
- [x] 明确“源码 Markdown 可读性优先”的修复原则：源文档保留 `$...$` / `$$...$$`，兼容逻辑下沉到 `docs/make.jl`
- [x] 在 `docs/make.jl` 中加入构建前数学定界符转换层，使 Documenter 仍接收兼容语法而源码保持标准写法
- [x] 将本地站点链接从目录式 pretty URL 调整为 `.html` 直达，修复首页左侧 tab 需要二次点击的问题
- [x] 显式配置 Documenter `edit_link = "master"`，消除分支推断环境告警
- [x] 完成 `docs/src/theory/` 全目录重建验证，确认理论文档源码与站点构建同时满足可读性与可渲染性

## 模块结构总览

```text
src/
|- Core/               # 核心类型、常量、配置、基础接口
|- Geometry/           # 网格、几何信息、输入输出与求积支撑
|- BasisFunctions/     # RWG、SWG、PWC、PWCHex、RBF
|- IntegralEquations/  # EFIE、MFIE、CFIE、VEFIE、SCFIE、PMCHW 等
|- FastAlgorithms/     # MLFMA、Lebedev、快速算子相关模块
|- Solvers/            # Direct、GMRES、BiCGSTAB、预处理
|- PostProcessing/     # RCS、远场、近场、天线后处理
|- Parallel/           # MPI、线程并行
|- IO/                 # 结果导出与数据写出
|- Utilities/          # 日志、轻量支撑工具、辅助算法
`- Driver.jl           # 统一仿真入口
```

## 历史阶段摘要

### Phase 1-6

- 完成项目骨架、核心模块划分与主求解链首版落地

### Phase 7-10

- 完成 Legacy 对齐、关键误差因子修复、性能主线优化与精度闭环

### Phase 11-12

- 完成 PWC / PWCHex / RBF 及混合体元支持

### Phase 13-17

- 完成 MPI 与性能探索、统一验证报告、发布流程规范化、依赖收敛和文档刷新

## 路线图使用原则

- 只保留当前仍有决策价值的路线信息
- 已完成且稳定的历史内容在进度文档中保留摘要，不再在此重复展开
- 每次实质进展后，同步更新本文件与 `REFACTORING_PROGRESS.md`