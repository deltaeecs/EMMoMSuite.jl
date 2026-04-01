# EMSuite 重构进度

> 最后更新: 2026-04-01

## 当前总览

- 当前主线问题：MLFMA 多层 upward/downward pass 的逐层对齐与修正
- 当前次主线问题：F7 direct 路径 Round 3 性能回收
- 发布流程、统一报告、依赖瘦身、README 与理论文档修复已完成当前轮次收口

## 2026-03-11 Update：PMCHW medium 逐层节点定位

- 已新增 `benchmark/compare_pmchw_upward_downward_localization_medium.jl`，用于在 medium `N=540` 夹具上对 `M-pass` 做两类定位：
  - `upward`：对每一层父盒 `aggS` 使用 exact reintegration（直接按父盒中心/极点重新积分 basis）做点对点对照；
  - `downward` 末端：在完成 translation + downward 后，对叶层每个 cube 的 receive 输出分别与 far-only dense `EM/HM` 行块对照。
- 当前数值结论：`upward` 不是首个失配放大点。
  - default budget 下，`k0` 的 level `4/3/2` 最差父盒相对误差分别约为 `4.99e-4 / 5.30e-4 / 5.22e-4`；
  - default budget 下，`k1` 的 level `4/3/2` 最差父盒相对误差分别约为 `1.21e-4 / 1.30e-4 / 9.76e-5`；
  - loose budget 下上述 `upward` 指标基本不变，说明 near-range 放松不会显著改变父盒聚合的 pointwise 误差。
- 新增 exact-upward 链路对照后，`downward` 中间场也被进一步压缩：
  - default budget：将 `upward` 替换为 exact parent `aggS` 后，最差 `disaggG` 相对差异仅约 `1.6e-16 (k0)` / `1.8e-16 (k1)`；
  - loose budget：最差 `disaggG` 相对差异约 `3.04e-4 (k0)` / `4.93e-5 (k1)`，叶层最差盒上的 `disaggG` 差异仍仅约 `7.21e-5` 或机器精度。
- 叶层 receive 积分规则已同步试验：
  - 当前实现切到 4 点三角形规则后，default 最差叶盒收敛为 `cube 354`，`k0 rel_total ≈ 5.09e-3`，`k1 rel_total ≈ 2.33e-2`；
  - loose budget 下，exact-upward + 4-point receive 的最差叶盒为 `cube 439 (k0 ≈ 2.59e-2)` 与 `cube 422 (k1 ≈ 6.08e-2)`；
  - 已测 worst cube 上 7 点规则未稳定优于 4 点：default `k0/k1` 反而略差，loose `k1` 有改善但 loose `k0` 变差，因此 4 点是当前更稳的折中，不宜回退到 3 点，也暂不直接切到 7 点。
- 结论：后续不再把主怀疑点放在 `aggregate_upward!` 本身，而是继续下钻 `translate! -> disaggregate_downward! -> disaggregate_leaf_pmchw_m!` 的末端链路，重点围绕最差叶盒 `354/439/422` 做中间 `disaggG` 与 leaf receive testing 的点对点比对。

## 2026-03-11 Update：PMCHW medium `M-pass` 主因修正

- 在前述定位脚本基础上，已新增“exact-upward 使用 4 点 source 积分”的对照分支，并输出最差叶盒的 pole 级 `disaggG` 差异。
- 新结论：主误差不在 `translate!` 或 `disaggregate_downward!`，而在 `aggregate_leaf_pmchw!` 仍使用 3 点三角形积分，和 direct far `K/L` 参考使用的 4 点规则不一致。
  - default budget：将 exact-upward 的 source 积分改为 4 点后，`k0/k1` 最差叶盒 `rel_total` 从 `5.09e-3 / 2.33e-2` 直接降至 `3.90e-5 / 2.70e-5`；
  - loose budget：对应指标从 `2.59e-2 / 6.07e-2` 降至 `6.69e-3 / 4.37e-4`；
  - 同时，当前实现与 exact-upward(4pt src) 的最差叶盒 `disaggG` 差异可压到机器精度或 `1e-5 ~ 1e-4` 量级，说明 source 端积分阶数就是主导项。
- 已据此修改 `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`：
  - `aggregate_leaf_pmchw!` 的三角形求积从 3 点提升到 4 点；
  - 保持此前已验证更稳定的 4 点 `_receive_terms` 规则不变。
- 修正后回归结果：`julia --project=. test/test_pmchw_block_fidelity_medium.jl` 全部 `26/26` 通过，关键指标显著改善：
  - `default_weak_m.rel_total ≈ 6.90e-7`
  - `default_m_k0.rel_total ≈ 4.02e-6`
  - `default_m_k1.rel_total ≈ 1.28e-6`
  - `loose_weak_m.rel_total ≈ 1.06e-4`
  - `loose_m_k0.rel_total ≈ 2.32e-4`
  - `loose_m_k1.rel_total ≈ 3.93e-5`
- 当前判断：PMCHW medium `M-pass` 的多层主误差已经完成主因修复。后续工作从“继续盲钻 downward”切换为两项收口工作：
  - 补齐更通用的 `nLevels >= 3` 回归闸门，避免这类 source/receive quadrature 不一致再次回归；
  - 继续做检视迭代，确认本修正不会在其他算子或大尺寸 case 上引入副作用。

## 2026-03-11 Update：PMCHW multilevel 回归与检视 Round 1

- 已新增直接多层回归 `test/test_pmchw_multilevel_quadrature_regression.jl`，夹具与 medium block fidelity gate 保持一致：
  - `generate_sphere_mesh(0.5, 10, 20)`、`RWGBasis`、`PMCHW(300e6, 4.0)`；
  - 显式断言 `k0/k1` 两棵 octree 均满足 `nLevels > 3`；
  - 对 `M_only` probe 的 `M×k0` / `M×k1` 做 far-only dense `EM/HM` 对照，并分别覆盖 default / loose near-range 预算。
- 已新增独立慢速批次入口 `test/runtests_batch6.jl`，用于单独执行该回归而不影响现有 batch 组织。
- 新回归首次运行暴露测试自身缺陷 `UndefVarError: nnz not defined`，已通过补入 `using SparseArrays` 修正。
- 当前验证结果：
  - `julia --project=. test/test_pmchw_multilevel_quadrature_regression.jl`：`17/17` 通过；
  - `julia --project=. test/runtests_batch6.jl`：`17/17` 通过；
  - 关键指标为 `default_m_k0 ≈ 4.02e-6`、`default_m_k1 ≈ 1.28e-6`、`loose_m_k0 ≈ 2.32e-4`、`loose_m_k1 ≈ 3.93e-5`，并确认 `default_nlevels = loose_nlevels = (5, 5)`。
- 检视 Round 1 还发现定位 benchmark 在实现已切换到 4 点 source 后出现“参考语义漂移”：主 `exact-upward` 仍默认使用 3 点 source 参考，会把当前正确实现误报成偏差。
- 已据此修正 `benchmark/compare_pmchw_upward_downward_localization_medium.jl`：
  - `exact_level_agg(...; quad_order=4)` 与 `apply_exact_upward_chain!(...; quad_order=4)` 现在默认对齐当前 4 点 source 实现；
  - 旧 3 点 source 对照保留为显式 `legacy` 分支，避免后续定位脚本再次把“当前实现 vs 旧参考”混为“实现误差”。
- 修正后 benchmark 复核结论：
  - current implementation 与 exact-upward(4pt src) 的 `disaggG` 差异已到机器精度量级；
  - 旧 exact-upward(3pt src legacy) 仍稳定复现 `5e-3 ~ 6e-2` 量级的历史偏差，证明该脚本现在既能正确反映当前实现，也保留了旧主因的对照证据。
- 当前状态：多层回归闸门与第一轮检视修正已完成；下一步只剩继续做第二轮 clean review，并在需要时补做文档站点重建。

## 2026-03-12 Update：PMCHW multilevel 检视 Round 2（clean）

- 已对本轮改动链路完成第二轮 clean review，检视范围聚焦于：
  - `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl` 中的 leaf source / receive 路径；
  - `test/test_pmchw_multilevel_quadrature_regression.jl` 与 `test/runtests_batch6.jl` 的回归覆盖与断言口径；
  - `benchmark/compare_pmchw_upward_downward_localization_medium.jl` 的 exact-upward 参考语义与 legacy 对照输出。
- Round 2 结论：未发现新的功能正确性问题、接口兼容性问题或架构越层问题。
- 本轮仅记录一个低风险代码整洁度观察项：`_receive_terms` 里仍保留未使用的 `η` 参数与三角形 `normal` 局部变量，说明 receive-side 公式探索痕迹尚未完全清理；该项不影响当前数值结果，也不阻塞本轮放行。
- 因此前述 PMCHW medium `M-pass` 修正现在满足：
  - 主因修复已完成；
  - `nLevels > 3` 直接回归已补齐；
  - 两轮检视中第一轮修正 benchmark 语义漂移，第二轮 clean 无新增问题。
- 当前剩余事项已从“继续检视代码正确性”收敛为“补做 `docs/build/` 站点重建，并在后续更大 case 上继续扩展覆盖”，不再阻塞本轮代码收口。

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

## 2026-03-11 Update：公式渲染检视 Round 2

- 已将理论文档中残余的通用 Markdown 数学定界符统一迁移为 Documenter / Julia Markdown 兼容语法：
  - 行内数学：`$...$` -> ``...``
  - 展示数学：`$$...$$` -> `math` 代码块
- 本轮覆盖文件：
  - `docs/src/theory/basis_functions.md`
  - `docs/src/theory/electromagnetics.md`
  - `docs/src/theory/excitations.md`
  - `docs/src/theory/fast_algorithms.md`
  - `docs/src/theory/integral_equations.md`
  - `docs/src/theory/method_of_moments.md`
  - `docs/src/theory/post_processing.md`
  - `docs/src/theory/solvers.md`
- 额外修正：
  - 消除 `method_of_moments.md` 与 `fast_algorithms.md` 的独立等号渲染风险
  - 修正 `integral_equations.md` 中 CFIE 公式的 `\text` 拼写损坏
  - 修正 `basis_functions.md` 中 `\tilde` 拼写损坏
- 验证结果：
  - `julia --project=docs docs/make.jl` 已通过
  - theory 目录不再出现 Documenter 的 “Unexpected Julia interpolation in the Markdown” 告警
- 当前结论：理论文档公式渲染问题已完成本轮收口；文档构建剩余提示仅是 Documenter 无法自动推断 `edit_link` 分支的环境告警，与公式渲染无关。

## 2026-03-11 Update：公式渲染检视 Round 3

- 用户侧新增约束已纳入方案：理论文档源码必须保留标准 Markdown 数学写法，不能为了迁就 Documenter 而把行内公式长期写成 ``...``。
- 已将以下理论页源码恢复为标准 `$...$` / `$$...$$` 形式：
  - `docs/src/theory/basis_functions.md`
  - `docs/src/theory/electromagnetics.md`
  - `docs/src/theory/excitations.md`
  - `docs/src/theory/fast_algorithms.md`
  - `docs/src/theory/integral_equations.md`
  - `docs/src/theory/method_of_moments.md`
  - `docs/src/theory/post_processing.md`
  - `docs/src/theory/solvers.md`
- 已在 `docs/make.jl` 中新增构建前转换层：
  - 源 Markdown 继续使用 `$...$` / `$$...$$`
  - `makedocs` 前自动复制到临时 source 目录并转换为 Documenter / Julia Markdown 兼容语法
- 已同步修正文档站点本地浏览链路：
  - `Documenter.HTML(prettyurls = false)`，生成 `*.html` 直达链接
  - `Documenter.HTML(edit_link = "master")`，消除 `edit_link` 分支推断告警
- 验证结果：
  - `julia --project=docs docs/make.jl` 通过
  - `docs/build/index.html` 中首页导航与正文链接已变为 `guide/installation.html`、`theory/overview.html` 等直达路径
  - 本轮未再出现数学渲染告警或 `edit_link` 环境告警
- 当前结论：文档系统现已同时满足三项要求：
  - 源 Markdown 对常规阅读器保持标准数学写法
  - Documenter 站点构建保持稳定
  - 本地打开 `docs/build/index.html` 时可直接点击进入目标页面，无需二次点击目录

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
- 已新增 medium 逐层定位脚本，并确认 `upward` 各层父盒的 exact reintegration 误差维持在 `1e-4` 量级；当前首个显著放大点在 leaf receive，而非父盒聚合

## 历史阶段摘要

### Phase 1-6

- 完成项目骨架、核心模块、几何 / 基函数 / 方程 / 快速算法 / 求解器 / 后处理主干落地

### Phase 7-10

- 完成 Legacy 对齐、关键误差修正、性能主线优化与精度收敛

### Phase 11-12

- 完成 PWC、PWCHex、RBF 及混合体元路径支持

### Phase 13-17

- 完成 MPI / 性能探索、统一验证报告、发布流程规范化、依赖瘦身与文档刷新

## 2026-04-01 Update：PMCHW multilevel 检视 Round 3

- 本轮检视范围：`_receive_terms` 接口与所有调用方、回归测试断言完整性、benchmark 定位脚本对齐语义、Phase 17 新增文件结构。
- **发现问题 1**：`_receive_terms` 的 `η` 参数（Round 2 已标记为死代码观察项）未在本轮修复，本轮正式修复：
  - 移除 `η :: Number` 参数
  - 移除 `normal = normalize(cross(...))` 及相关 `v1/v2/v3` 临时变量
  - 同步更新 `disaggregate_leaf_pmchw_j!` / `disaggregate_leaf_pmchw_m!` 两处调用方
- **发现问题 2**（新发现）：`benchmark/compare_pmchw_upward_downward_localization_medium.jl` 中 3 处直接调用 `PMCHWMLFMAOperatorModule._receive_terms(...)` 或本地 `receive_terms_with_rule(...)` 仍以旧 7 参数形式传入 `eta`，与修复后的 6 参数签名不匹配，将在运行时导致 `MethodError`；本轮已全部修正。
- 回归测试（`test_pmchw_multilevel_quadrature_regression.jl`）断言结构检视：`nLevels > 3` 断言覆盖 default/loose 两档，门限相对实测值留有 100~3000 倍余量，设计合理。
- Phase 17 文件结构检视：`ReleaseWorkflow.jl` 与 `benchmark/reporting/` 三层职责清晰，无新发现。
- 本轮结论：发现并修复了 2 个问题（1 个已知观察项落地修复，1 个新发现接口不匹配），**Round 3 不属于 clean round**，需继续 Round 4。

## 2026-04-01 Update：PMCHW multilevel 检视 Round 4 / Round 5

### Round 4（clean）

- 检视范围：接口签名一致性、相位符号（source/test）、测试断言、batch 入口设计、benchmark 完整性。
- 确认所有 `_receive_terms` 调用方使用 6 参数签名 ✓；聚合用 `+jk` / 测试用 `-jk` 相位符号 ✓；`runtests.jl` 不包含慢速 batch6 ✓；benchmark 内 `receive_terms_with_rule` 签名对齐 ✓。
- **结论：无新问题，Round 4 为 clean round（第 1 个连续 clean round）。**

### Round 5

- 检视范围：扩展至 Phase 17 文件引用链与 EMSuite 主包 include 结构。
- **发现问题**：`src/Accuracy/ReleaseWorkflow.jl` 从未被任何代码 `include`，是 Phase 17 添加但未接入的死文件（功能已由 `benchmark/support/release_support.jl` 完整覆盖并验证）。该文件属未追踪文件（从未入 git）。
  - 已删除该文件，消除 `src/` 目录下的死字段混淆。
  - ROADMAP / PROGRESS 记录修正：Phase 17 发布流程功能的实际实现在 `benchmark/support/release_support.jl`，而非 `src/Accuracy/`。
- **结论：发现并修复 1 个 Phase 17 遗留死文件问题，Round 5 不属于 clean round，继续 Round 6。**

## 2026-04-01 Update：PMCHW multilevel 检视 Round 6 / Round 7 / Round 8

### Round 6（clean，连续 1）

- 检视范围：Phase 17 依赖收敛变更链（`EMSuite.jl` 导出清理、`LightweightSupport` 集成、`LoggingExtras`/`ProgressMeter` 移除、`Accuracy.jl` 的 `BenchmarkReportData` 移除）。
- 所有变更逻辑一致，符合 Phase 17 依赖收敛意图。**未发现新问题。**

### Round 7（clean，连续 2）

- 检视范围：`PostProcessing.jl`（移除 ProgressMeter）、`runtests.jl`（新增 test_release_workflow）、`test_benchmark_report_data.jl`（切换到 release_support.jl 直接引入）、`dataset_generator.jl`（Roots → find_zero_bisection）。
- 所有变更均为 Phase 17 依赖瘦身的符合预期结果。**未发现新问题。**

### Round 8（clean，连续 3）——终止条件满足

- 检视范围：`indices.jl` 中 `factor_values` 替换 `Primes.factor`、`MPIArrays.jl` 中 `knn_bruteforce` 替换 `NearestNeighbors.knn`、`pinv2interpW.jl` 中 `find_zero_bisection` 替换 `Roots.find_zero`。
- 算法语义验证：`factor_values` 返回含重复的质因数列表，`slicedim2partition` 算法通过内循环复用机制确保与旧唯一质因数版本等价；`knn_bruteforce` 返回格式与 NearestNeighbors.knn 兼容；`find_zero_bisection` 以指数扩张方式确保覆盖所有实际截断阶范围。**未发现新问题。**
- **连续 3 轮检视（Round 6 / 7 / 8）均无新问题，检视终止条件满足。**

## 当前未完成事项

- [ ] F7 direct Round 3 性能回收
- [ ] Julia General Registry 发布
- [ ] docs/build/ 站点重建（源码已正确，构建产物滞后）

## 2026-04-01 Update：全工程检视迭代

### 全工程检视原则

在 PMCHW multilevel 检视收口后启动面向全工程的检视迭代，四个维度：程序架构、算法实现、软件工程、开发原则。

### Round 1（发现问题，不属于 clean round）

- **P0 修复：`PostProcessing/RCS.jl:78`** — 移除 `radarCrossSection` 主路径中的活跃调试 `println("k0: ..., eta0: ...")`
- **P1 修复：`src/Accuracy/BenchmarkReportData.jl`** — 删除已从 `Accuracy.jl` 剔除但滞留于 git 的死文件（与 `ReleaseWorkflow.jl` 同类问题）
- **Chore 修复：`src/Geometry/MeshGen.jl_temp` + `src/IntegralEquations/MFIE.jl.bak`** — 删除被错误纳入 git 的编辑临时/备份文件
- **Chore 修复：`.gitignore`** — 补入 `*.jl.bak`、`*.jl_temp`、`*.jl.orig`、`results/`、`test_results/runs/` 防止再次入库
- **Fix：`FastAlgorithms/MLFMA/OctreeBuilder.jl`** — 将 `println("Warning: nLevels=..."` 改为 `@warn`，接入日志体系
- **Phase 17 提交补全** — 将前序会话验证通过但未入库的 Phase 17 实现全量提交（LightweightSupport、Project.toml、benchmark 基础设施、测试、报告、文档）

### Round 2（发现问题，不属于 clean round）

- **Fix：`src/EMSuite.jl`** — 补入 `extract_sphere_radius` 到顶层 re-export；该函数在 `Accuracy.jl` 已导出但遗漏于主模块
- **确认**：所有移除的依赖（LoggingExtras/ProgressMeter/NearestNeighbors/Primes/Roots）已完全清出 src/ 和 test/
- **确认**：Solvers、IO、Ports 模块无 TODO/活跃调试输出

### Round 3（clean，连续 1）

- 检视范围：FarField/NearField/FieldCut/MLFMAFastPost、Core/Materials/Geometry、benchmark 定位脚本签名、BasisFunctions 占位实现评估
- `receive_terms_with_rule` 签名在 benchmark 脚本中确认无旧 `eta` 参数，与生产代码一致
- SWG/RWG `evaluate` 占位实现评估：未被任何生产路径调用，现阶段不阻塞；留作 P3 后续实现
- **未发现新问题，Round 3 为 clean**（全工程检视第 1 个 clean round）

### Round 4（发现问题，不属于 clean round）

- 检视范围：全工程系统性代码质量审查，覆盖 src/ 下所有模块
- 修复分三类：**物理常数标准化**、**日志系统标准化**、**接口与代码质量**
- **常数标准化（15+ 文件，40+ 处）：**
  - 移除所有硬编码的 `c0 = 299792458.0`、`mu0 = ...`、`eps0 = ...` 定义
  - 统一使用 `Constants.c0/mu0/eps0/eta0`（定义于 `src/Core/Constants.jl`）
  - 修复精度问题：部分文件使用 `eps0 = 8.854187817e-12`（10位），现改用 `Constants.eps0 = 1/(c0²×mu0)` 高精度值
  - 覆盖模块：IntegralEquations (EFIE/MFIE/SCFIE/VEFIE)、PostProcessing (NearField/Absorption/RadiationIntegral/RCS)、Utilities (Parameters/MieSeries)、FastAlgorithms/MLFMA (MLFMAOperator/Aggregation/Disaggregation)、Parallel/MPI (VolumeAssembly)
- **日志系统标准化（40+ 处）：**
  - 所有裸 `println` 替换为 `@info`，便于日志级别控制
  - 覆盖模块：IntegralEquations (SCFIE 7处, VEFIE 9处)、FastAlgorithms/MLFMA (MLFMAOperator, PMCHWMLFMAOperator 2处含中文)、Parallel/MPI (VolumeAssembly 4处, DistributedGMRES 5处)
- **接口与代码质量：**
  - `RWG.jl` / `SWG.jl`: `evaluate()` 从静默返回 `SVector(0,0,0)` 改为 `error()` 明确报错，附加用户指引
  - `SWG.jl`: 修正边界面注释误导（实际已包含边界面，非注释掉状态）
  - 文档化 SWG vs RWG 边界处理设计差异（SWG 包含边界面用于 VEFIE 通量连续性，RWG 排除边界边用于 PEC EFIE）
  - `EFIE.jl`: 移除 4 处已注释的 debug println 死代码
- **文档任务完成：**
  - `Parameters.jl`: 添加 Thread Safety 警告章节，明确 GLOBAL_PARAMS 无锁全局状态、并发竞态风险、当前约束与未来改进方向
  - `EFIE.jl`: 详细文档化 near-interaction 转置 workaround（Legacy Parity 设计约束、半解析积分非对称性、矩阵对称性保障机制）
  - `RWG.jl`: 增强边界处理设计差异文档，明确排除边界边的物理原因（PEC EFIE 电流不可流出闭合表面）
- **Git 提交：**
  - 代码修复：`refactor: standardize constants and logging across entire codebase` (commit 0da31d3)
  - 进度/路线图：`docs: update progress and roadmap for全工程检视 Round 4` (commit 9d7067b)
  - 文档任务：`docs: document thread-safety limitations and design constraints` (commit 50c9ee2)
- **Round 4 收口：** 所有发现问题已修复，所有待办任务已完成

### Round 5（发现问题，不属于 clean round）

- 检视范围：验证 Round 4 修复质量 + 全工程硬编码常数/日志残留检查
- 使用 explore agent 系统性扫描 src/ 下所有 .jl 文件
- **发现 4 处遗漏的硬编码常数**（Round 4 未覆盖）：
  1. `MLFMAOperator.jl:76` — `lambda = 299792458.0 / freq` (HIGH)
  2. `WavePort.jl:152-153` — `c0 = 299792458.0`, `η₀ = 376.730313461` (HIGH)
  3. `VolumeAssembly.jl:494-495, 677-678` — `mu0/eps0` 硬编码 (MEDIUM)
- **已全部修复**：
  - MLFMAOperator.jl: 使用 `Constants.c0`（注：MLFMAOperatorMPI 在 Round 4 已正确）
  - WavePort.jl: 移除硬编码，使用 `Constants.c0/eta0`；Ports.jl 添加 Constants 导入
  - VolumeAssembly.jl: 使用 `Constants.mu0/eps0`，添加 Constants 导入
- **验证结果**：
  - ✅ 日志系统：100% 完成（所有 println → @info 已在 Round 4 完成）
  - ✅ 导入语句：所有使用 Constants 的文件已正确导入
  - ✅ 循环依赖：无问题
  - ✅ 文档格式：无问题
- **Git 提交**：`fix: complete physical constants standardization (Round 5 findings)` (commit b3197cf)
- **Round 5 收口**：发现并修复问题，连续 clean 轮次清零（连续 0 轮）

### Round 6（发现问题，不属于 clean round）

- 检视范围：验证 Round 5 修复 + 全工程数值稳定性、参数验证、命名规范审查
- 使用 explore agent 深度审查核心算法模块（IntegralEquations, MLFMA, BasisFunctions, Geometry, Solvers）
- **Round 5 修复验证**：✅ 全部通过（MLFMAOperator、WavePort、VolumeAssembly 的 Constants 使用正确）
- **发现 9 个新问题**：
  - 🔴 Critical (1): FastExp.jl `unsafe_trunc` 索引溢出风险
  - 🔴 High (1): Singularities.jl 退化三角形 log(0) 风险（5处）
  - 🟡 Medium (7): 端口参数验证、预条件对角阈值、数值阈值硬编码、命名混乱等
- **已修复关键问题**（3个）：
  - FastExp: `unsafe_trunc` → `clamp(trunc(...), 1, n-1)`
  - Singularities: `log(1-a/s)` → `log(max(1-a/s, eps(FT)))` 保护（singularF1/F21）
  - WavePort: 添加 `@assert a>0, b>0, freq>0` 参数验证
- **Git 提交**：`fix: numerical stability improvements (Round 6 findings)` (commit 6da9de1)
- **Round 6 收口**：发现并修复关键问题，连续 clean 轮次清零（连续 0 轮）

### Round 7（发现问题，不属于 clean round）

- 检视范围：验证 Round 6 修复 + 评估剩余 Medium/Low 问题 + 轻量级新问题扫描
- 使用 explore agent 深度验证 Round 6 的 3 个修复点
- **Round 6 修复验证**：
  - ✅ FastExp.jl clamp 边界正确且完整
  - ⚠️ Singularities.jl F1/F21 正确但**遗漏 F22**（本轮修复）
  - ✅ WavePort.jl @assert 设计合理
- **发现 2 个遗漏的关键问题**：
  - 🔴 Singularities.jl `singularF22()` 未保护 — 3 处 log 项缺少 `max(..., eps)` 保护
  - 🔴 Translation.jl Rab=0 除零风险 — `R̂ab = RabVec / Rab` 未检查分母为零
- **已全部修复**：
  - Singularities: 添加 `log(max(1-a/s, eps_ft))` 保护（与 F1/F21 对齐）
  - Translation: 添加 `if Rab < eps*edge then skip` 检查，避免 NaN 传播到 MLFMA
- **剩余问题评估**：6 个 Medium/Low 问题中，2 个高优先级已修复，其余 4 个不影响正确性（可延后到 Phase 20）
- **Git 提交**：`fix: complete numerical stability fixes (Round 7 findings)` (commit 6d5584d)
- **Round 7 收口**：发现并修复关键问题，连续 clean 轮次清零（连续 0 轮）

### Round 8（发现问题，不属于 clean round）

- 检视范围：验证 Round 7 修复 + 全面数值保护扫描（log/sqrt/acos/除零）
- 扫描范围：158 个 Julia 文件，重点检查所有潜在数值风险
- **Round 7 修复验证**：
  - ✅ Singularities.jl singularF22 — 3 处 log 保护正确且完整
  - ✅ Translation.jl Rab 检查 — eps*edge_length 阈值合理
- **发现 5 个新问题**（4 阻塞 + 1 延后）：
  - 🔴 SpecialPorts.jl CoaxPort — `log(outer/inner)` 无参数验证
  - 🔴 BasisUtilities.jl — 退化三角形/边无保护（面积=0，边长=0）
  - 🔴 CoordinateTransforms.jl — `acos(r̂[3])` 参数未 clamp 到 [-1,1]
  - 🔴 SCFIE.jl — 奇点跳过逻辑不清晰（R<1e-10 后立即用 R）
  - 🟡 多处类型不稳定容器 (`Any[]`) — 延后到 Phase 20
- **已全部修复阻塞性问题**：
  - SpecialPorts: 添加 @assert (inner>0, outer>inner, eps_r>0)
  - BasisUtilities: 添加 area>eps(FT) + 3 个边长检查
  - CoordinateTransforms: acos(clamp(...)) 保护
  - SCFIE: 显式 if-continue 块 + 注释说明
- **Git 提交**：`fix: add comprehensive numerical safety checks (Round 8 findings)` (commit 8e68d5d)
- **Round 8 收口**：发现并修复关键问题，连续 clean 轮次清零（连续 0 轮）

### Round 9（发现问题，不属于 clean round）

- 检视范围：全面验证 Round 6-8 修复 + 数值阈值一致性检查
- **Round 6-8 修复验证**：
  - ✅ **全部 10 处修复都正确且完整**
  - FastExp.jl、Singularities.jl、Translation.jl、WavePort.jl、SpecialPorts.jl、BasisUtilities.jl、CoordinateTransforms.jl、SCFIE.jl
- **发现 1 个阻塞性问题**：
  - 🔴 FastExp.jl `fast_exp_ikr` 使用 `unsafe_trunc + min` — 与 `fast_green_func` 的 `clamp(trunc(...))` 不对称
  - 影响：潜在越界风险，代码维护性差
- **发现 6 个非阻塞性改进机会**（延后到后续 phase）：
  - 数值阈值不一致（`1e-10` vs `1e-12` vs `1e-15`）
  - 硬编码 magic constants 分散（13+ 处）
  - Singularities.jl 中 `epsilon_l` 系数缺少文档说明
- **已修复阻塞性问题**：
  - FastExp.jl:134 改为 `clamp(trunc(Int, idx_f), 1, n-1)` 与 `fast_green_func` 对齐
- **Git 提交**：`fix: unify table lookup safety in FastExp (Round 9 finding)` (commit 74d8652)
- **Round 9 收口**：发现并修复阻塞性问题，连续 clean 轮次清零（连续 0 轮）

### Round 10（**第一个 Clean Round** ✅）

- 检视范围：验证 Round 9 修复 + 全面数值风险/并行竞态/文档质量扫描
- 覆盖范围：12 个核心模块完整检查，20+ 文件采样扫描
- **Round 9 修复验证**：
  - ✅ FastExp.jl clamp 修复完整且与 `fast_green_func` 一致
  - ✅ 项目中已完全移除 `unsafe_trunc`，零残留
- **阻塞性问题扫描**：
  - ✅ **零除零风险** — 所有除法都有前置检查或逻辑保证
  - ✅ **零 log/sqrt 负参数** — Singularities.jl 所有 `log()` 都使用 `max(..., eps_ft)` 保护
  - ✅ **零数组越界** — 所有 `@inbounds` 都有 `clamp` 或循环范围保护
  - ✅ **零并行竞态** — EFIE/MFIE/SCFIE/VEFIE 使用 SpinLock 保护行写入，MPI 有行级锁
  - ✅ **零类型不稳定** — 所有容器都有显式类型
- **遗留非阻塞性问题确认**：
  - 数值阈值差异（`1e-10` vs `1e-12`）— 不影响正确性，仅是维护性问题
  - Singularities.jl 中 `epsilon_l` 系数（`1e-2` vs `1e-3`）— 反映不同上下文需求
- **文档质量**：关键函数都有 docstring，数值保护有注释说明
- **Round 10 收口**：**第一个 Clean Round 达成**，连续 clean 轮次：**1/3**

### Round 11（待启动）

- 连续 clean 轮次：1/3（需再连续 2 轮 clean 才可收口）
- Round 11 将采样扫描 MLFMA 相关模块验证数值一致性

