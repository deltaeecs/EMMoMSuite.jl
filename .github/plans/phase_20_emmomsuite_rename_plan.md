# Phase 20 EMSuite -> EMMoMSuite 全工程改名计划

> 创建日期: 2026-04-15
> 状态: 规划中
> 关联目标: 将当前包从 `EMSuite.jl` 完整改名为 `EMMoMSuite.jl`，新名称既明确体现 Method of Moments，又与重构前 Legacy 体系 (`MoM_Kernels` / `MoM_Basics` / `MoM_AllinOne`) 以及当前过于泛化的 `EMSuite` 做出清晰区分，同时满足 Julia 包命名与模块命名的一致性要求。

---

## 0. 开发原则遵循声明

本 Phase 适用并遵循 `.github/copilot-instructions.md` 中以下原则：

- 原则 1（TDD 工作流）：先建立改名前失败、改名后通过的导入/文档/发布 smoke gate，再进行批量重命名。
- 原则 4（进度同步）：本计划执行后必须同步 `REFACTORING_ROADMAP.md` 与 `REFACTORING_PROGRESS.md`。
- 原则 5（提交规范）：包元数据、源码入口、脚本迁移、文档迁移、发布/仓库层更名分别独立闭环提交。
- 原则 6（终端输出规范）：验证通过终端直接采集，不做 shell 重定向。
- 原则 7（检视迭代）：改名完成后必须做至少 2 轮 clean review，重点检查遗漏引用、文档失真与环境失配。
- 原则 9（计划文档规范）：本计划明确给出命名基准、影响边界、DoD、检视与回链。

---

## 1. 命名决议

### 1.1 采用名称

- Julia 包名：`EMMoMSuite`
- 顶层模块名：`EMMoMSuite`
- 仓库与对外文档展示名：`EMMoMSuite.jl`

### 1.2 采用原因

1. `MoM` 被直接嵌入名称，能明确表达该工程不是泛化的电磁工具箱，而是以矩量法为主线的 CEM 框架。
2. `EM` + `MoM` + `Suite` 保留当前“工程套件/工作流集合”的定位，不会误导为单一算子或教材式最小实现。
3. 与重构前 `MoM_Kernels`、`MoM_Basics`、`MoM_AllinOne` 有连续性，但不与任一 Legacy 包撞名。
4. 相比 `MethodOfMoments.jl` 这类完全通用名，更容易体现本仓库是“工程化、可发布、含 benchmark/release 链”的重构后主包。
5. 作为 Julia 包名可合法使用 PascalCase，且与顶层模块保持一一对应。

### 1.3 明确不采用的命名方向

- 不采用 `EMSuite.jl`：语义过宽，无法突出 MoM 主线，也不利于与早期重构阶段区分。
- 不采用 `MoMSuite.jl`：虽然更短，但弱化了电磁语境，过于接近其他潜在数值/统计领域命名。
- 不采用 `MethodOfMoments.jl`：名字虽直白，但过于通用，不利于与本仓库既有工程资产和历史路线图衔接。

---

## 2. 重命名边界与影响面

本次改名不是单点替换，而是覆盖五个层面：

### 2.1 包与源码入口层

- `Project.toml` 的 `name`
- `src/EMSuite.jl` 文件名与 `module EMSuite`
- `src/Driver.jl` 等内部 `using ..EMSuite` / `using EMSuite.*` 引用
- 任何依赖顶层模块名的 `include`、docstring、导出说明

### 2.2 子环境与构建层

- `docs/Project.toml` 中对主包的依赖名
- `docs/make.jl` 中 `using EMSuite`、`modules = [EMSuite]`、`sitename = "EMSuite.jl"`
- `docs/Manifest.toml` 与其他 path-based manifest 不手工硬改，而是在执行阶段通过重新 instantiate / resolve 刷新
- `benchmark/` 子环境中所有 `using EMSuite` 与限定模块导入

### 2.3 测试、脚本、基准与发布链

- `test/` 中所有入口与子模块限定路径
- `benchmark/` 中所有脚本、report、配置说明、生成文件列名中对包名的显式引用
- `run_tests.jl`、`scripts/`、`docs/` 构建脚本
- `.github/RELEASE_CHECKLIST.md` 中的 registry/publish/install 文案

### 2.4 文档与历史工程记录层

- `README.md`、`CHANGELOG.md`、`docs/src/**`、`docs/INTEGRATION_COMPARISON.md`
- `.github/REFACTORING_ROADMAP.md`、`.github/REFACTORING_PROGRESS.md`
- `.github/plans/**` 与老 Phase 计划文档

处理原则：

1. 当前使用说明、可执行指令、API 文档、发布清单必须全部切换到 `EMMoMSuite`。
2. 历史计划/进度中若是“引用当前包名作为现行事实”，应同步改名；若是“描述历史阶段当时名为 EMSuite 的事实”，允许保留原文，但需要在本 Phase 文档中说明 `EMSuite` 为旧名。
3. 所有 `using EMSuite` 形式的可执行入口必须清零，避免仓库内部继续依赖旧名。

### 2.5 仓库与对外发布层

- 本地工程目录名、远程仓库名、README clone 指令、Documenter 站点标题
- General Registry 发布清单与后续 `Pkg.add("EMMoMSuite")` 安装说明

说明：仓库目录与远程仓库改名不在本仓库源码 patch 内一步完成，但必须纳入本 Phase 计划的最终 rollout 步骤。

---

## 3. 兼容性决策

### 3.1 不提供同包双顶层名兼容

本 Phase 不计划在同一个 Julia 包内长期维持 `using EMSuite` 与 `using EMMoMSuite` 两套顶层入口并行，原因如下：

1. Julia 包注册、`Project.toml` 名称、`src/<PackageName>.jl` 入口和用户导入名天然绑定，双名并行会带来额外维护成本与歧义。
2. 仓库内部已有大量脚本和文档可统一迁移，不需要为内部调用保留旧名。
3. 本次改名的目标本就包含“与旧名明确切割”，保留旧导入名会削弱这一目标。

### 3.2 可接受的过渡手段

- 在 `README.md`、发布说明、CHANGELOG 中显式声明：`EMSuite.jl` 已更名为 `EMMoMSuite.jl`。
- 在本 Phase 结束时新增一节“rename migration notes”，统一列出外部用户需修改的命令与导入方式。
- 如后续确有外部兼容需求，可另建极薄迁移包或在 release note 中给出映射，不在本 Phase 主线内实现。

---

## 4. TDD 与执行分阶段方案

### Phase 20.A 基线冻结与失败门建立

目标：先把“什么算改名完成”冻结成可验证 gate。

建议 gate：

1. `julia --project=. -e "using EMMoMSuite"` 在改名前应失败，改名后必须通过。
2. `julia --project=docs docs/make.jl` 在改名后必须通过。
3. release/benchmark quick smoke 至少保留一条，例如 `benchmark/run_release_workflow.jl benchmark/configs/release_quick.toml`。
4. 全仓搜索中，不应再存在可执行路径上的 `using EMSuite`。

建议新增/调整测试：

- 若当前没有 package-name smoke test，可在 `test/` 新增一个轻量导入测试，显式断言顶层模块名和关键导出可访问。
- 对 docs / benchmark 采用命令级验证，不强行把构建流程塞入单元测试。

### Phase 20.B 元数据与顶层源码改名

目标：先完成最核心的包入口切换。

涉及文件：

- `Project.toml`
- `src/EMSuite.jl` -> `src/EMMoMSuite.jl`
- `src/Driver.jl`
- 其他直接 `using ..EMSuite` 的源码文件

完成标准：

- 主环境下 `using EMMoMSuite` 成功。
- 所有源码内部限定路径均指向 `EMMoMSuite`。

### Phase 20.C 仓库内脚本与子环境迁移

目标：让 docs、test、benchmark、scripts 全部跟上新包名。

涉及范围：

- `docs/Project.toml`
- `docs/make.jl`
- `docs/src/**`
- `benchmark/**`
- `test/**`
- `run_tests.jl`
- `scripts/**`

执行重点：

1. 先替换导入语句与模块限定路径。
2. 再替换当前面向用户的标题、说明、命令示例。
3. path-based manifests 统一通过 resolve/instantiate 重新生成，不手写大量锁文件差异。

### Phase 20.D 文档、路线图与历史记录治理

目标：把仓库的“当前身份”统一到 `EMMoMSuite.jl`，同时避免篡改历史语义。

执行策略：

1. `README.md`、`CHANGELOG.md`、`docs/src/**`、`.github/RELEASE_CHECKLIST.md` 完全切到新名。
2. `.github/REFACTORING_ROADMAP.md`、`.github/REFACTORING_PROGRESS.md` 顶部与当前状态描述切到新名，并新增“旧名 EMSuite”说明。
3. `.github/plans/**` 中仅对“仍被当前工作流引用、且表达当前事实”的内容改名；对历史叙述可保留，并在本计划中统一解释。

### Phase 20.E 发布与仓库层 rollout

目标：完成仓库外部入口与用户迁移说明。

应包含：

1. 仓库目录/远端仓库改名为 `EMMoMSuite.jl`
2. 文档站点、clone 指令、README badge/链接更新
3. General Registry 发布说明改为 `Pkg.add("EMMoMSuite")`
4. release note / changelog 中新增 rename migration section

---

## 5. 文件分组与建议提交粒度

### Commit 1: 包元数据与入口改名

- `Project.toml`
- `src/EMMoMSuite.jl`
- `src/Driver.jl`
- 相关顶层源码导入修正

建议提交信息：`refactor: rename package entry to EMMoMSuite`

### Commit 2: 测试、benchmark、脚本迁移

- `test/**`
- `benchmark/**`
- `run_tests.jl`
- `scripts/**`

建议提交信息：`refactor: migrate test and benchmark imports to EMMoMSuite`

### Commit 3: docs 与用户文档迁移

- `README.md`
- `CHANGELOG.md`
- `docs/**`
- `.github/RELEASE_CHECKLIST.md`

建议提交信息：`docs: rename public documentation to EMMoMSuite`

### Commit 4: 治理文档与历史引用收口

- `.github/REFACTORING_ROADMAP.md`
- `.github/REFACTORING_PROGRESS.md`
- `.github/plans/**` 中需要跟随当前事实的文档

建议提交信息：`docs: register EMMoMSuite rename across planning docs`

### Commit 5: 仓库/发布层 rollout

- 非源码层操作、最终说明和回归结果补记

建议提交信息：`chore: finalize EMMoMSuite repository rollout`

---

## 6. 风险清单与排查顺序

### 6.1 高风险项

1. `docs/Project.toml` 与主包名不一致，导致 `docs/make.jl` 无法解析主模块。
2. `benchmark/` 与 `test/` 中存在深层 `using EMSuite.Submodule` 漏改，导致仅局部脚本失败。
3. `Manifest.toml` / `docs/Manifest.toml` 中 path 和 name 不一致，造成环境解析异常。
4. `.github/plans/**` 和 `README.md` 中的命令示例仍保留旧名，形成新的错误文档。
5. 仓库目录和远端仓库未同步改名，导致 clone / dev path 文档与实际不一致。

### 6.2 推荐排查顺序

1. 先修包入口：确保 `using EMMoMSuite`。
2. 再修 test/benchmark/docs 的导入语句。
3. 然后重建 docs / 子环境 manifests。
4. 再清理 README / release checklist / current-facing docs。
5. 最后审查历史文档，区分“现行事实”与“历史叙述”。

---

## 7. 检视迭代计划

Phase 完成后至少进行以下两轮检视；若发现新问题，继续迭代直到连续 2 轮 clean：

### Round 1: 运行面检视

- `using EMMoMSuite`
- 主测试批次 smoke
- docs build
- release quick workflow
- grep 检查残留 `using EMSuite`

### Round 2: 文档与发布面检视

- README 安装/开发/测试命令逐条核对
- `.github/RELEASE_CHECKLIST.md` 与 registry 名称核对
- `.github/REFACTORING_ROADMAP.md` / `REFACTORING_PROGRESS.md` 当前状态段落核对
- `.github/plans/**` 是否存在误导性“当前包名仍为 EMSuite”的现行描述

若连续两轮均未发现新问题，才允许宣布本 Phase 收口。

---

## 8. DoD（完成定义）

本 Phase 完成需同时满足：

1. 主包可通过 `using EMMoMSuite` 导入。
2. 仓库内部所有可执行脚本、测试、benchmark、docs build 均切换到 `EMMoMSuite`。
3. 当前对外文档中不再把 `EMSuite.jl` 当作现行包名。
4. `README.md`、`CHANGELOG.md`、`.github/RELEASE_CHECKLIST.md` 明确声明旧名到新名的迁移关系。
5. `.github/REFACTORING_ROADMAP.md` 与 `.github/REFACTORING_PROGRESS.md` 已同步记录本 Phase 的启动与完成。
6. 至少完成 2 轮 clean review。
7. 完成仓库与发布层 rollout 清单，包含 `Pkg.add("EMMoMSuite")` 口径。

---

## 9. 待机链接

- 路线图更新位置：`REFACTORING_ROADMAP.md` -> `3. 文档与发布链维护` 下新增 `Phase 20` 改名任务
- 进度更新位置：`REFACTORING_PROGRESS.md` -> `2026-04-15 Update` 改名规划启动记录
- 本计划文件：`.github/plans/phase_20_emmomsuite_rename_plan.md`
