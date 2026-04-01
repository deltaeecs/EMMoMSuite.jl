# Phase 17 发布流程规范化与报告产物整备计划

> 创建日期: 2026-03-11
> 状态: 已完成
> 关联目标: 在 Phase 16 图表化统一报告闭环基础上，补齐发布前流程规范、配置输入、日志与产物目录、报告脚本分层，以及远场绘图风格统一，形成可重复执行的 release workflow。

> 2026-03-11 决议补记：
> - 配置格式固定为 TOML
> - 允许引入 CairoMakie，但必须遵守“新增依赖最小化”约束，不作为无条件核心运行时依赖扩散
> - 极坐标远场图纳入正式交付，而不是仅保留可选接口
> - 已知遗留项采用独立 registry 维护
>
> 2026-03-11 收口补记：
> - `Plots` 已隔离到 `benchmark/Project.toml`，不再作为主包运行时依赖
> - release/report 支撑逻辑已迁出主包 API 面，统一收敛到 `benchmark/support/release_support.jl`
> - `LoggingExtras`、`Primes`、`Roots`、`NearestNeighbors`、`ProgressMeter` 已由仓库内轻量实现替代
> - 验证闭环：`using EMSuite`、两组发布链测试、release report 与 quick workflow 全部通过

---

## 0. 开发原则遵循声明

本 Phase 适用并遵循以下原则（来自 `.github/copilot-instructions.md`）：

- 原则 1（TDD 工作流）：先冻结配置/产物/日志 schema，再改脚本与测试。
- 原则 2（严格 Legacy 对齐）：数值结果的真值仍只来自 `MoM_Kernels` / `MoM_Basics` / `MoM_AllinOne`、FEKO CSV 与解析基准；流程规范化不得改变物理口径。
- 原则 3（差异排查流程）：若规范化后结果变化，必须按几何-常数-积分-组装顺序排查，不得用经验系数补偿。
- 原则 4（进度同步）：本计划执行后必须同步 `REFACTORING_ROADMAP.md` 与 `REFACTORING_PROGRESS.md`。
- 原则 5（提交规范）：配置 schema、日志/产物目录、报告脚本重构、绘图美化分别独立闭环。
- 原则 6（终端输出规范）：不做 shell 重定向；持久化日志由 Julia 脚本内部写入。
- 原则 7（检视迭代）：Phase 完成后至少进行 2 轮 clean review。
- 原则 9（计划文档规范）：本计划明确 Legacy/参考锚点、DoD、检视与回链。
- 原则 10（新增约束）：新增第三方依赖必须从严控制，优先使用 Julia 标准库与现有依赖，避免为发布流程引入大面积兼容风险。

---

## 1. 现状诊断

当前发布验证链已经能生成 `test_results/reports/RELEASE_VALIDATION_REPORT.md` 与图表，但流程仍存在以下工程缺口：

1. 执行入口分散：`benchmark/run_full_benchmark.jl`、`benchmark/run_full_accuracy_benchmark.jl`、`benchmark/performance_baseline.jl`、各组 `benchmark/accuracy/run_*.jl`、`benchmark/run_release_validation_report.jl` 各自管理输入与输出，缺少统一的顶层 profile。
2. 配置不可复用：大多数脚本仍将频率、求解器、门限、输出目录硬编码在脚本内部，难以复跑同一 release profile。
3. 日志结构不统一：多数 benchmark 仍以 `println/@printf` 直接输出到终端，只有部分脚本沉淀 CSV 或 markdown，缺少统一的 run manifest、状态摘要和错误记录。
4. 产物目录未按 run 隔离：当前主要写入 `test_results/accuracy/` 与 `test_results/reports/`，适合单次人工运行，不适合多轮 release 候选并存。
5. 绘图风格尚未统一：Phase 16 现用 `Plots.jl + GR` 可出图，但视觉风格与 `MoM_Visualizing` 中的 2D/极坐标风格不一致。
6. 报告脚本职责偏重：`benchmark/run_release_validation_report.jl` 同时承担数据发现、图表绘制、markdown 写作与发布结论拼装，后续继续扩展会越来越脆。

---

## 2. Legacy / 参考锚点

### 2.1 数值真值锚点

- Legacy 网格/夹具：`MoM_AllinOne/meshfiles/`
- Legacy/历史 CSV：`test_results/legacy_baseline/`
- FEKO 基线：`../MoM_AllinOne/deps/compare_feko/`
- 解析基线：`src/Utilities/MieSeries.jl` 与已有 A1-A4/B1-B5 端口链路

### 2.2 流程与图形风格锚点

- 统一报告入口现状：`benchmark/run_release_validation_report.jl`
- 性能结构化出口现状：`benchmark/performance_baseline.jl`
- 精度单项脚本现状：`benchmark/accuracy/run_F1_F4_jet.jl`、`run_F5_F6_sphere.jl`、`run_F7_F9_plate.jl`、`run_P1_P3_pmchw.jl`、`run_B1_B5_antenna.jl`
- 可视化风格参考：`../MoM_Visualizing/src/far_field_viz.jl`、`../MoM_Visualizing/src/lines_plot.jl`、`../MoM_Visualizing/src/polar_plot.jl`

### 2.3 从 MoM_Visualizing 抽取的风格约束

本 Phase 远场绘图样式以以下特征为参考，而不是继续沿用默认 GR 风格：

1. 2D 线图采用明确的 theme 与统一轴样式，而不是每张图临时拼参数。
2. 远场对比图优先使用“实线 = EMSuite，散点/空心标记 = FEKO/Mie/Legacy”这一比较语义，减少双实线重叠后的可读性问题。
3. 轴标签使用物理量命名：`theta (deg)`、`RCS (dBsm)`，并保持 5 个主刻度左右的紧凑出版风格。
4. 图例背景透明、无边框，避免遮挡主瓣/零陷附近信息。
5. 导出尺寸固定为论文/报告友好的小图规格，避免当前大图在 markdown 中缩放后字太小或太密。
6. 对需要方向图呈现的用例，保留升级到极坐标图的接口，而不仅限于笛卡尔切面图。

---

## 3. 目标产物

本 Phase 交付物分为五层：

1. 顶层 release profile 配置文件：定义本轮发布要跑哪些用例、门限、绘图风格、输出根目录、是否允许复用已有 CSV。
2. 标准化 run artifact 目录：一次运行对应一个唯一 run id，包含 manifest、logs、metrics、plots、report 五类产物。
3. 统一日志与状态摘要：至少有 machine-readable 的 `run_manifest.toml/json` 与 `run_status.csv/json`，而不是只看终端文本。
4. 分层报告生成脚本：拆分为“收集/汇总/绘图/写报告”四段，减少单脚本职责堆叠。
5. 美化后的远场绘图输出：默认风格对齐 `MoM_Visualizing` 的 2D 比较图语义，并正式交付极坐标视图。
6. 依赖收敛清单：梳理当前发布链额外依赖，优先移除不必要第三方库、避免把报告链依赖扩散到核心求解链。

---

## 4. 目录与文件规范（拟定）

### 4.1 配置目录

新增目录：`benchmark/configs/`

建议内容：

- `release_default.toml`：默认发布验证 profile
- `release_quick.toml`：轻量 smoke profile
- `plot_style.toml`：统一绘图风格参数
- `known_exceptions.toml`：已知遗留项与 release waiver
- `thresholds.toml`：精度/性能/端口门限

### 4.2 标准 run 目录

新增根目录：`test_results/runs/<run_id>/`

建议子目录：

- `manifest/`
  - `run_manifest.toml`
  - `environment.json`
- `logs/`
  - `orchestrator.log`
  - `<case_id>.log`
- `metrics/`
  - `accuracy/*.csv`
  - `performance/*.csv`
  - `summary/run_status.csv`
- `plots/`
  - `accuracy/*.png`
  - `accuracy_polar/*.png`
  - `performance/*.png`
- `report/`
  - `RELEASE_VALIDATION_REPORT.md`

### 4.3 稳定别名目录

保留当前对外稳定路径，但改为“最新一次成功 run 的镜像/复制目标”：

- `test_results/accuracy/`
- `test_results/reports/`

这样既保留现有脚本兼容性，也避免多轮 release 产物互相覆盖。

---

## 5. 脚本分层规范（拟定）

### 5.1 顶层编排

- 新增单一入口，例如：`benchmark/run_release_workflow.jl`
- 职责：读取 profile、生成 run id、调度 accuracy/performance/report 三条链、写 manifest 和状态摘要

### 5.2 用例执行层

- 保留现有 `benchmark/accuracy/run_*.jl` 与性能脚本的数值主体
- 逐步将硬编码参数替换为从 profile 读入
- 每个脚本的职责限定为：执行用例并输出标准 metrics，不直接决定最终报告结构

### 5.3 产物汇总层

- 抽离统一 artifact collector，例如 `src/Accuracy` 下继续扩展 report-data / manifest-data 读取接口
- 统一解析 accuracy/performance/antenna/parallel artifacts

### 5.4 绘图层

- 将远场图与性能图从 report writer 中拆出为独立 plotting module / helper
- 绘图层只消费标准化 metrics，不直接扫描整个目录做猜测式发现
- 远场绘图默认同时生成笛卡尔切面图与极坐标图，两者共享同一 style schema

### 5.6 依赖边界

- 发布流程新增能力优先依赖 Julia 标准库与仓库既有依赖，避免引入 YAML、额外 CLI 桥接库等非必要第三方组件
- CairoMakie 若用于发布图，应限制在“报告/绘图层”，不得把图形依赖扩散到核心求解主路径
- 若现有依赖只为旧报告链或可替换工具服务，需在本 Phase 中列出可移除/可降级清单

### 5.5 报告写作层

- markdown writer 只负责：章节布局、指标填表、图链插入、release recommendation
- 不再直接负责 CSV 搜索、分类和绘图细节

---

## 6. 日志与状态规范（拟定）

### 6.1 每次 run 必须记录

- 启动时间、git commit、Julia 版本、线程数、MPI rank 数
- 使用的 profile 文件与阈值文件版本
- 每个 case 的开始时间、结束时间、状态（PASS/FAIL/SKIP/ERROR）
- 每个 case 的关键指标摘要（RMSE、MaxErr、Zin gap、总耗时）

### 6.2 日志分级

- `INFO`：阶段开始/结束、路径、耗时摘要
- `WARN`：历史产物缺失、某项仅复用旧数据、已知遗留项
- `ERROR`：脚本异常、CSV 结构不匹配、关键产物缺失

### 6.3 最低 machine-readable 产物

- `run_manifest.toml`：输入配置、环境、选项
- `run_status.csv`：每个 case 的状态与摘要指标
- `artifact_index.csv`：本次 run 产物清单与相对路径
- `known_exceptions.toml`：本次 release 接受的 waiver 及依据

---

## 7. 绘图规范（拟定）

### 7.1 远场对比图

- 默认输出两类：
  - 笛卡尔切面图（主交付）
  - 极坐标图（正式交付）
- 默认语义：
  - EMSuite：实线
  - FEKO / Mie / Legacy：采样散点或空心标记
- 默认视觉参数：
  - 统一字体大小与线宽
  - 透明无边框图例
  - dotted grid
  - 统一导出尺寸与 dpi
- 默认信息标注：
  - 曲线标题保留 case、cut、reference、RMSE
  - 误差摘要放在副标题或 caption，避免覆盖主图

### 7.2 性能图

- 总耗时柱状图保持现有信息密度
- 阶段拆分图改为更清晰的分组/堆叠规范，并补全单位、线程/MPI 说明
- 若数据量增加，优先引入按 solver family 分组，而不是继续单图硬塞全部 case

### 7.3 绘图库迁移策略

- 首选方向：在不扩大核心依赖面的前提下，把发布图迁移到 CairoMakie，以便更接近 `MoM_Visualizing` 风格与导出质量
- 约束条件：CairoMakie 只服务于 release plotting path，不作为求解主路径的通用硬依赖扩散
- 保守过渡：若一次性迁移成本过高，则先在现有结构上抽象 plot style schema，第二步再替换 backend

---

## 8. 缺失功能补项清单

除用户显式要求外，本 Phase 认为还缺以下能力：

1. run manifest：否则报告无法追溯“这张图是哪次执行产出的”。
2. artifact index：否则 markdown 之外的 CSV/PNG 仍是散落文件。
3. profile-based rerun：否则每次换机器或换线程都要手改脚本。
4. latest vs archived run 机制：否则下一轮复测会覆盖上一轮发布候选证据。
5. case status matrix：明确 PASS/FAIL/SKIP/KNOWN_EXCEPTION，而不是只在报告正文散落描述。
6. known-exception registry：把 `F2_CFIE_Jet_Direct` 这种已接受遗留项从正文中抽离为明确的 release waiver 列表。
7. plot caption/metadata：在报告中说明参考源、门限、线程/MPI 条件，避免图表脱离上下文。
8. 依赖收敛审计：识别哪些第三方库只用于旧图表/旧脚本，是否可迁出主依赖或延后加载。

---

## 9. 分步实施计划

### 9.1 Step A：冻结 schema 与目录

- 定义 profile schema
- 定义 run manifest / run status / artifact index schema
- 定义 run 目录结构与 latest 镜像策略

### 9.2 Step B：接通顶层 orchestrator

- 新增统一入口脚本
- 接入现有 accuracy/performance/report 脚本
- 生成 run id 与基础 manifest/log/status

### 9.3 Step C：抽离 plotting/report helpers

- 拆分 `benchmark/run_release_validation_report.jl`
- 让 collector / plotting / writer 三层独立

> 2026-03-11 进展补记：以上目标已落地，当前已形成：
> - `benchmark/reporting/collector.jl`
> - `benchmark/reporting/plotting.jl`
> - `benchmark/reporting/writer.jl`
> - `benchmark/run_release_validation_report.jl` thin wrapper

### 9.4 Step D：latest 镜像策略收敛

- 先将 stable alias 同步进 run 工件，保证兼容历史脚本
- 再将本次 run 工件回灌到 stable alias，明确 stable 目录只表示“最近一次成功 run”
- 避免 stable 目录继续承担唯一源职责

> 2026-03-11 进展补记：以上镜像双向链路已在 `benchmark/run_release_workflow.jl` 落地。

### 9.4 Step D：远场图美化

- 对齐 `MoM_Visualizing` 的 line/scatter compare 语义
- 验证 Jet / Sphere / Plate / PMCHW 四组图的可读性
- 同步产出极坐标版并验证其在 markdown/报告中的版式

### 9.5 Step E：依赖收敛

- 盘点当前发布链新增/现有第三方依赖
- 优先移除或隔离非必要依赖，特别是只服务于报告链但当前挂在核心依赖面的组件
- 明确 CairoMakie、Plots、MoM_Visualizing 在最终发布流程中的边界

### 9.6 Step F：检视与发布前确认

- 做两轮 clean review
- 汇总待确认项
- 获得确认后再作为 release pipeline 默认入口

---

## 10. 已确认决议与剩余待确认项

### 10.1 已确认决议

1. 配置格式固定为 TOML，优先使用标准库 `TOML`。
2. 允许引入 CairoMakie，但必须限制在发布图路径，不能无约束扩散依赖面。
3. 极坐标远场图纳入正式交付。
4. 已知遗留项单独维护 `known_exceptions.toml/md`。

### 10.2 剩余待确认项

1. `test_results/reports/` 是否继续作为“最新结果稳定入口”？当前建议保留，但改为从 `test_results/runs/<run_id>/report/` 镜像生成。
2. CairoMakie 的接入方式采用直接依赖、延迟加载还是扩展式隔离；当前建议优先选择隔离式接入。
3. 现有 `Plots` / `MoM_Visualizing` 依赖在 Phase 17 完成后是否还能进一步收敛出核心包主依赖面。

---

## 11. DoD（完成定义）

满足以下条件才可标记本 Phase 完成：

- 至少 1 个 release profile 能从单一入口完整驱动 accuracy/performance/report 全链路
- 已形成标准 run 目录，包含 manifest、logs、metrics、plots、report
- 当前 `test_results/reports/` 与 `test_results/accuracy/` 仍保持兼容输出
- 远场图默认风格已完成一次统一升级，并通过 Jet / Sphere / Plate / PMCHW 四组图的笛卡尔版与极坐标版实测验证
- 统一报告已接入 run manifest、artifact index、case status matrix、known exception 摘要
- 已完成一次发布链依赖审计，并明确非必要依赖的移除/隔离结论
- 至少 2 轮检视迭代完成，且连续两轮无新增阻塞问题

---

## 12. 检视迭代计划

### Round 1（工程结构）

- 检查配置 schema 是否足以表达当前所有主用例
- 检查 run 目录与 latest 镜像是否会破坏现有脚本兼容性
- 检查 report/plot/data collector 是否已经完成职责解耦

### Round 2（发布可用性）

- 抽样复核 manifest、status、artifact index 与实际产物是否一致
- 复核美化后的笛卡尔图与极坐标图在 markdown 中的可读性
- 复核 known exception 机制是否能清晰表达 release waiver
- 复核依赖收敛后，发布图链与核心求解链是否已保持边界清晰

连续两轮无新增阻塞问题后，方可把该流程作为发布默认入口。

---

## 13. 回链

- Roadmap 勾选位置：新增 `Phase 17` 条目
- Progress 更新点：新增“Phase 17 发布流程规范化与报告整备”日志
- 前置成果：`.github/plans/phase_16_release_validation_report_plan.md`
- 现有统一报告：`test_results/reports/RELEASE_VALIDATION_REPORT.md`

---

## 14. 2026-03-11 首轮实现结果

### 14.1 已落地产物

- `src/Accuracy/ReleaseWorkflow.jl`
- `benchmark/run_release_workflow.jl`
- `benchmark/configs/release_default.toml`
- `benchmark/configs/release_quick.toml`
- `benchmark/configs/plot_style.toml`
- `benchmark/configs/thresholds.toml`
- `benchmark/configs/known_exceptions.toml`
- `test/test_release_workflow.jl`

### 14.2 已落地能力

1. 单一 TOML profile 入口已可生成标准 `test_results/runs/<run_id>/` 目录。
2. 已可自动写出：
  - `manifest/run_manifest.toml`
  - `metrics/summary/run_status.csv`
  - `metrics/summary/artifact_index.csv`
3. 统一报告已补齐：
  - 极坐标远场图
  - `Case Status Matrix`
  - `Known Exceptions`
4. 发布建议已从硬编码 F2 逻辑切换到 `known_exceptions.toml` 驱动。
5. 已完成一项确定性的依赖收敛：`MoM_Visualizing` 已从 `Project.toml` 主依赖移除。

### 14.3 本轮验证

- `julia --project=. test/test_release_workflow.jl` 通过（13/13）
- `julia --project=. test/test_benchmark_report_data.jl` 通过（19/19）
- `julia --project=benchmark benchmark/run_release_validation_report.jl` 通过
- `julia --project=. benchmark/run_release_workflow.jl benchmark/configs/release_quick.toml` 通过

### 14.4 当前未闭环项

1. `release_default.toml` 目前默认仍以“复用既有 accuracy 产物 + 重跑性能/报告”为主，尚未完成所有 accuracy 子脚本的 profile 参数化接入。
2. 发布图链现已要求通过 benchmark 环境执行，后续如需进一步统一入口体验，需要决定是否增加包装脚本或自动环境提示。
