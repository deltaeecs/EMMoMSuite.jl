# Phase 15 N-Muller Dense Baseline Plan

> 创建日期: 2026-03-07  
> 状态: dense baseline + 小球/中尺度对照 + DeltaGap 最小链路已完成  
> 关联主计划: `.github/PHASE_15_DIELECTRIC_ANTENNA_PLAN.md`

---

## 0. 开发原则遵循声明

本子计划遵循 `copilot-instructions.md` 中以下原则：

- 原则 1（TDD）：先落 `test/test_nmuller.jl`，再补 formulation 实现。
- 原则 2（Legacy / 文献对齐）：不使用经验系数；连续式与离散块组合以 scuff-em 文档和 Yla-Oijala/Taskinen 2005 为对照。
- 原则 3（差异排查）：先锁 dense block 组合，再比较 PMCHW vs N-Muller conditioning / solve 行为。
- 原则 7（Phase 结束检视）：本子流完成后必须做至少 2 轮 clean review。
- 原则 9（计划文档规范）：本计划明确给出文献锚点、DoD 与治理回链。

## 1. Legacy / 文献对齐基准

- 外部连续式锚点：scuff-em `lsInnardsSections/Formulations.tex` 中 N-Muller `NMuller0 / NMuller / M^(1) - M^(2)` 描述。
- formulation 文献锚点：Yla-Oijala, Taskinen 2005, Well-conditioned Muller formulation。
- 代码内对照锚点：现有 `src/IntegralEquations/PMCHW.jl` 的 L/K 子块装配与 `src/IntegralEquations/PMCHWBlockOperators.jl` 的 RWG surface Gram。

## 2. 本轮目标

本轮只做最小 dense baseline：

1. 新增 `NMuller` formulation 类型。
2. 复用现有 PMCHW 子块装配，形成 `M^(1) - M^(2)` 的 dense 2N×2N 组装。
3. 新增正式测试，锁定：
   - 尺寸/有限性
   - 非平凡 block 结构
   - 与 PMCHW 不是同一矩阵
   - 可执行 direct solve
4. 在同一 dielectric sphere 上补一个正式对照入口，统一记录 `cond`、GMRES residual 与 `rel_vs_LU`。

## 3. 非目标

- 本轮不接 MLFMA backend。
- 本轮不提前宣告 conditioning 改善已经验收。
- 本轮不把 N-Muller 接入天线 B1-B5 脚本。

## 4. DoD

- `test/test_nmuller.jl` 通过。
- `test/test_nmuller_comparison.jl` 通过。
- `NMuller` 已从顶层 `EMSuite` 导出。
- roadmap / progress / Phase 15 主计划已同步指向本子流。
- 小球 dense 对照脚本 `benchmark/compare_pmchw_nmuller_sphere.jl` 已能稳定输出 PMCHW / N-Muller conditioning 与 GMRES 指标。

## 5. 2026-03-07 实施结果

- `test/test_nmuller.jl`：`15/15 pass`
- `test/test_nmuller_comparison.jl`：`9/9 pass`
- `test/test_nmuller_excitation.jl`：`10/10 pass`
- 共享球夹具（`freq=120 MHz`, `radius=0.1 m`, `N=54`）当前观测：
   - `cond(PMCHW)   = 1.711099e7`
   - `cond(NMuller) = 2.538759e5`
   - `cond ratio    = 1.483700e-2`
   - `GMRES PMCHW   : iters=200, res=4.901509e-1, rel_vs_LU=1.019453`
   - `GMRES NMuller : iters=200, res=1.236587e-1, rel_vs_LU=2.729159e-1`
- 当前解释边界：
   - 这证明当前 dense N-Muller 基线并非仅仅“可组装”，而是在共享球夹具上已表现出显著优于 PMCHW 的线性系统性质。
   - 本结论当前仅覆盖 dense / 小球 / 共享 RHS 对照，不替代后续 medium-scale 或 backend 级验收。

- `benchmark/compare_pmchw_nmuller_sphere.jl medium` 当前也已给出第二组探针结果：
   - 对应专门回归入口：`test/test_nmuller_comparison_medium.jl`
   - `N=150`
   - `cond(PMCHW)   = 6.102107e7`
   - `cond(NMuller) = 2.254102e5`
   - `cond ratio    = 3.693974e-3`
   - `GMRES PMCHW   : res=7.078038e-1, rel_vs_LU=1.000066`
   - `GMRES NMuller : res=4.174158e-2, rel_vs_LU=7.808047e-2`
- 当前解释边界补充：
   - 这说明 N-Muller 相对 PMCHW 的优势并不局限于最小夹具。
   - 该现象现已被 `test/test_nmuller_comparison_medium.jl` 固化为正式中尺度 dense 专门门禁，但仍未并入默认 `runtests`。

- 作为对“larger fixture + 输入阻抗路径”分支的补充，当前已新增：
   - `excitation_vector(op::NMuller, source::DeltaGapSource, basis::RWGBasis)`
   - `input_impedance(op::NMuller, source::DeltaGapSource, I_2N, basis::RWGBasis)`
   - 对应回归：`test/test_nmuller_excitation.jl` 当前 `10/10 pass`
- 当前解释边界再补充：
   - 这一步补齐的是 N-Muller 在 Phase 15 中进入最小馈电/输入阻抗工作流所需的接口，并未把 N-Muller 直接接入 B1–B5 正式脚本。
   - 现阶段它更适合作为 formulation 对照与输入阻抗链路探针，而不是替代 PMCHW 主线天线基准。
   - 最新 `benchmark/compare_pmchw_nmuller_impedance.jl` 诊断显示：当前 DeltaGap / `input_impedance` 语义仍未校准，不能把该链路视为已验证的物理端口模型。
   - 后续若要升级为正式阻抗门，必须先明确 `NMuller` 的 port RHS 与 feed-current 提取应对应哪一组未知量或哪种组合缩放。

## 6. 回链

- Roadmap: Phase 15 alternative formulation 段
- Progress: 2026-03-07 N-Muller dense baseline 启动日志
- 主计划: `.github/PHASE_15_DIELECTRIC_ANTENNA_PLAN.md`