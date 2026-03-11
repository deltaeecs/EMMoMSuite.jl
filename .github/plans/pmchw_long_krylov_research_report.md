# PMCHW MLFMA 长 Krylov 保真边界调研报告

> 日期: 2026-03-07  
> 适用阶段: Phase 15  
> 关联文档: `.github/plans/phase_15_theory_impl_test_refresh.md`, `.github/REFACTORING_ROADMAP.md`, `.github/REFACTORING_PROGRESS.md`

> 2026-03-07 状态补记：本报告中的若干“建议动作”现已部分落地。PMCHW 顶层已迁移到 `PMCHWBlockOperator` shell，默认 `strong_form` 已接入真实 RWG block pairing，`PMCHWMLFMAOperator` 已可通过 `MatrixFreePMCHWBackend` 收编为 backend；中尺度 `N=540` 的四路 Gate S 专门回归 `dense weak / dense strong / MLFMA weak / MLFMA strong` 已正式跑通，结果为 `15.S2 | 9/9 pass | 9m05.1s`。因此本报告当前用途已从“提出是否应补 strong-form/shell”转为“说明为什么这些动作是必要的，以及剩余边界仍在哪里”。

> 2026-03-07 补记：alternative formulation 的第一轮 dense 基线也已落地。当前已有 `NMuller` formulation、`test/test_nmuller.jl`、`test/test_nmuller_comparison.jl` 与 `benchmark/compare_pmchw_nmuller_sphere.jl`。在共享小球夹具上，`cond(NMuller)/cond(PMCHW) = 1.4837e-2`，且在同一 GMRES 预算下 `NMuller` 的 residual 与 `rel_vs_LU` 均优于 `PMCHW`。这意味着“先补 alternative formulation 基线”已不再是建议，而是已完成的事实；本报告里的下一步应相应前推到 medium-scale / backend 误差预算层。

> 2026-03-07 新补记：预算子流现已与长 Krylov 子流正式接通。`benchmark/compare_pmchw_mlfma_budget_krylov.jl` 在 `N=540` medium 夹具上固定比较 `default / loose_near / fixed_leaf_0p04_nr9` 三组预算。结果显示：`short` 档下三组预算相对 dense 的 `Z_in` gap 仍只有亚欧姆量级；而 `long` 档下该 gap 被放大到 `default = 7.28Ω / 2.01Ω`、`loose_near = 15.93Ω / 2.04Ω`、`fixed_leaf_0p04_nr9 = 14.76Ω / 2.13Ω`。因此本报告当前剩余边界已进一步收敛为：如何治理 MLFMA backend fidelity / budget，而不是继续笼统讨论 solver 截止。

> 2026-03-07 再补记：`test/test_pmchw_mlfma_budget_krylov_medium.jl` 已把其中最有信息量的 `default / loose_near` long-Krylov 对照升级为正式专门门禁，并在当前环境下通过 `16/16 pass | 31m06.3s`。这意味着“长 Krylov 下 budget 会放大成主导项”已经不再只是 benchmark 观察，而是带门限的正式验收事实。

## 1. 调研目标

本报告针对 EMSuite 当前已确认的问题进行外部调研与工程对照：

- 小尺度 PMCHW 主线已经 green，说明 formal correctness 主干已打通。
- 中尺度 `N=540` 上，`Gate C` 与 `GMRES(200)` parity 说明算子在常规残差水平下可工作。
- 但当 Krylov 子空间做大到足以逼近 `LU` 时，dense GMRES 可以继续向 `LU` 收敛，而 PMCHW MLFMA 仍保留非平凡阻抗偏差。

因此本次调研要回答的问题不是“有没有另一个简单预条件器”，而是：

1. 文献和开源程序通常如何处理 transmission / dielectric PMCHW 一类问题的高精度求解。
2. 它们是优先修 formulation、operator preconditioning、space pairing，还是只在代数层堆 Krylov + 预条件器。
3. 对 EMSuite 当前“长 Krylov 下的算子保真边界”最有参考价值的路线是什么。

## 2. EMSuite 当前已知事实

### 2.1 已被正式验证的事实

- 小尺度正式门禁 `GD0`、`GD1`、`GD2`、`GB`、`Gate C`、`B2` 均已通过。
- 中尺度 `N=540` 上，`Gate C` 绿色，说明随机向量 matvec 一致性已到可接受水平。
- `dense+GMRES(200)` 与 `MLFMA+GMRES(200)` 接近，说明在较低精度阶段两者行为一致。
- `benchmark/compare_pmchw_mlfma_budget_krylov.jl` 已进一步证明这一点具有预算条件：在 `short` 档下，代表性三组预算相对 dense 的 `Z_in` gap 仅 `0.03`–`0.27Ω`（实部）与 `0.18`–`0.54Ω`（虚部），预算影响仍被 Krylov 截止误差淹没。
- PMCHW shell + strong-form 主线已经不是假设，而是正式实现：Dense shell、小夹具 MLFMA shell、中尺度 dense Gate S 与中尺度四路 Gate S 都已有正式测试入口。
- `test/test_pmchw_gate_s_mlfma_medium.jl` 已证明在统一 shell 语义、统一 RHS、统一 GMRES 参数下，`dense strong` 与 `MLFMA strong` 的残差几乎重合；当前关键指标为 `relw = 3.45e-3`、`rels = 3.18e-3`、`zgw = 2.34e-5`、`zgs = 5.75e-5`。
- Phase 15 天线基准脚本中的 B2 路径也已从占位变为实跑：`benchmark/accuracy/run_B1_B5_antenna.jl B2` 现可直接给出 PMCHW MLFMA 对 Direct 的输入阻抗对照，并已在当前环境下通过。
- `NMuller` dense baseline 已完成正式接线，并已有共享球夹具对照结果：`cond(PMCHW)=1.711099e7`、`cond(NMuller)=2.538759e5`、`cond ratio = 1.483700e-2`；同一 `GMRES(200)` 预算下，`NMuller` 的 `res=1.236587e-1`、`rel_vs_LU=2.729159e-1`，均优于 `PMCHW` 的 `res=4.901509e-1`、`rel_vs_LU=1.019453`。
- `benchmark/compare_pmchw_nmuller_sphere.jl medium` 也已跑通，说明这种 formulation 优势不只停留在最小夹具：`N=150` 时 `cond(PMCHW)=6.102107e7`、`cond(NMuller)=2.254102e5`、`cond ratio = 3.693974e-3`；同一 `GMRES(250)` 预算下，`NMuller` 的 `res=4.174158e-2`、`rel_vs_LU=7.808047e-2`，而 `PMCHW` 仍为 `res=7.078038e-1`、`rel_vs_LU=1.000066`。
- 该 medium 结果现在已经由 `test/test_nmuller_comparison_medium.jl` 固化为专门中尺度 dense 门禁，因此 N-Muller 的优势不再只是一次性 benchmark 观测。
- `benchmark/compare_pmchw_mlfma_budget_krylov.jl medium` 已显示：当 Krylov 深度继续增大时，MLFMA budget/fidelity 会从“可忽略的次要项”变成显著误差源。具体来说，`long` 档下代表性三组预算相对 dense 的 `Z_in` gap 已增至 `7.28Ω`、`15.93Ω`、`14.76Ω`（实部），这说明 medium 长 Krylov 剩余边界已不能再归因于 solver 截止本身。

### 2.2 已被正式否定的路线

- `Diagonal`、`ILU(Z_near)`、`BlockJacobi(op)` 在当前 PMCHW 中尺度夹具上都没有改善结果，反而更差。
- 继续扩大 `near_range` 会在当前夹具上触发 OOM，因此“把更多项塞回近场”不是主线修复方案。
- 因而当前问题不能再表述为“还差一个预条件器接口”。

### 2.3 当前最准确的问题表述

当前主线问题应定义为：

> 在 transmission / dielectric PMCHW 问题上，现有 MLFMA 压缩算子在随机向量和低精度 GMRES 阶段已足够，但在长 Krylov、高精度逼近 `LU` 时，其算子误差或 formulation-level  conditioning 问题会暴露出来，形成可观的阻抗偏差下界。

这一定义与公开文献的主流结论一致：对 Maxwell transmission 问题，单靠“快速 matvec + 朴素 GMRES”通常不够，必须联合 formulation、pairing、operator preconditioning 或更严格的压缩误差控制。

## 3. 文献结论

### 3.1 Calderon 体系的核心结论

Bempp 团队论文《Software frameworks for integral equations in electromagnetic scattering based on Calderón identities》给出的关键信息是：

- 过去十多年，Maxwell 边界积分方程的重要进展不只是更快的乘法，而是更稳定的离散空间配对与 Calderon 型算子框架。
- 论文强调 stable dual pairing、preconditioned electric field / magnetic field / combined field formulations，以及 product algebra 对复杂 block operator 的封装。
- 这类方法的重点不是“再试一个代数预条件器”，而是通过算子恒等式和合适的 trial-test pairing 改善谱性质。

对 EMSuite 的直接启示是：

- 当 dense 算子能继续逼近 `LU` 而 fast 算子在长 Krylov 阶段失真时，应优先怀疑 formulation + operator scaling + discretisation pairing，而不是只看 Krylov 参数。
- Transmission 问题中的高精度求解，公开主流路线更接近“算子级预条件/正则化”，而不是“近场块 ILU”。

### 3.2 强形式 strong form 的意义

Bempp handbook 明确写到：

- `use_strong_form=true` 等价于对系统施加 mass-matrix preconditioning。
- 在很多情形下它会显著降低迭代次数。

这说明公开实现里，连最基础的迭代入口都默认把“离散算子应该先映射到正确 range / dual pairing 下”当成第一层预条件，而不是把 weak form 生系统原样交给 GMRES。

对 EMSuite 的启示是：

- 这条路线对 EMSuite 是正确的，且现已落地：PMCHW shell 已具备 `weak_form/strong_form/strong_form_rhs/recover_trial_coefficients`，因此 weak/strong 不再停留在概念讨论层。
- 当前证据表明 strong-form 确实起作用，但它没有把问题“直接消失”成纯 solver 话题；它的价值在于先把“谱条件差”与“快速算子逼近误差”分离开。

### 3.3 文献对 dielectric transmission formulation 的一般趋势

scuff-em 内部文档在 dielectric object 路径中明确同时保留：

- PMCHWT formulation
- N-Muller formulation

并直接引用 Yla-Oijala 与 Taskinen 2005 的《Well-conditioned Muller formulation for electromagnetic scattering by dielectric objects》。

这说明在实际工程系统中，PMCHW 不是唯一默认答案。公开软件会保留更好条件数的 transmission formulation 作为替代路线。

对 EMSuite 的含义非常直接：

- 如果目标是高精度 dielectric scattering，不应把 PMCHW 当成唯一 formulation 主线。
- 至少需要把 N-Muller 或另一种 well-conditioned transmission formulation 纳入对照。

## 4. 开源实现对照

## 4.1 Bempp

### 实现特征

- Maxwell 路径基于 blocked multitrace operator。
- 提供 `strong_form()`，并在 GMRES 接口中显式支持 `use_strong_form=true`。
- Maxwell 多迹算子在代码层面就是 block operator，而不是把 transmission 系统拍平成一个无语义的大矩阵。
- 提供 OSRC 类近似边界算子，用于把连续问题的非局部映射替换为更好的近似边界算子。

### 八叉树 / FMM 树结构如何处理

- Bempp 的用户接口不直接暴露“PMCHW 专用双八叉树”这一概念，而是把树结构下沉到 ExaFMM 后端。
- 公开参数显示其 FMM 主控量是 `depth`、`expansion_order`、`ncrit`、`near_field_representation`，这意味着树的深度、每盒临界粒子数、展开阶数由通用 FMM 后端统一管理，而不是由 PMCHW 单独维护一套显式八叉树调度。
- 对 Maxwell transmission，Bempp 不是先造一个“2N 未知量的 PMCHW MLFMA 专用树”，而是先构造 electric-field、magnetic-field 等 Maxwell block operator，再把每个非局部块交给 FMM evaluator。
- 换句话说，Bempp 的树是“kernel/operator evaluator 级别”的公共基础设施，不是 EMSuite 当前这种“PMCHW 算子自己拥有两棵八叉树并显式调四遍”的架构。

### 矩阵向量乘积如何处理

- Bempp 的 FMM matvec 不是黑盒 `A*x` 一步到位，而是“变换到点值/散度表示 + FMM evaluator + singular sparse correction”的组合。
- 对 Maxwell electric boundary operator，公开实现的核心结构是：
	- 先把 RWG 系数通过 `compute_rwg_basis_transform` / `compute_rwg_div_transform` 映射成若干源分量；
	- 再对每个分量调用 `fmm_interface.evaluate(...)`；
	- 最后用 dual-space 映射把结果投回离散未知量，并显式加上 `singular_part @ x`。
- Maxwell magnetic boundary operator 同样不是单个“PMCHW 磁场专用八叉树核”，而是通过多个 FMM 评估结果重组 curl，再加 `singular_part @ x`。
- 这和 Bempp 对 scalar hypersingular 的处理一致：奇异/近奇异部分单独稀疏装配，非局部部分走 FMM evaluator，然后在线性算子层面组合。

### 对 PMCHW / transmission 的直接含义

- Bempp 的 transmission 主线是 multitrace block algebra：`[[M, E], [-E, M]]` 这一类 block 结构在算子层持续保留。
- 因此它的 PMCHW 风格求解更接近“4 个 block operator 的代数组合，各块内部再各自调用 FMM matvec”，而不是“一个 PMCHW 专用算子内部手写 4 遍远场流程”。
- 这类设计天然更容易接 strong form、mass-matrix preconditioning、Calderon product algebra，因为 block 语义始终没有丢。

### 公开示例反映的工程取舍

- dielectric Maxwell 示例 `examples/maxwell/maxwell_dielectric.py` 默认直接使用 `lu`，并明确写明“大问题应使用 iterative solvers with preconditioning”。
- 也就是说，公开示例没有暗示“裸 GMRES + 快速乘法”足以成为 transmission 问题的统一方案。

### 对 EMSuite 的可借鉴点

- 用 block operator 语义保留物理分块，而不是把一切都压成裸 `mul!`。
- 强制提供 strong-form 主线作为 transmission solve 的正式入口。
- 若要继续做 Krylov 主线，应优先做 Calderon / multitrace 风格的 operator-level regularization 对照，而不是继续堆局部代数预条件器。

## 4.2 FreeFEM + BEMTool + Htool

### 实现特征

- FreeFEM 文档明确其 BEM 路线采用 BEMTool + Htool，也就是 H-matrix / low-rank admissible block 压缩，不是 MLFMA。
- Maxwell BEM 路径使用 RT0 型空间生成器，说明它把空间和 kernel 类型作为实现核心。
- 代码中同时存在 GMRES 与 FGMRES，以及 `HMatVirtPrecon` 一类专门针对 H-matrix 的预条件器/求解器包装。
- H-matrix 构建显式暴露 `eps`、`eta`、cluster size、target/source depth 等压缩控制参数。

### 树结构如何处理

- FreeFEM/BEMTool/Htool 的核心几何层不是 MLFMA 八叉树，而是 cluster tree。
- 它们通过 admissibility 条件把块划分为近场稠密块和远场低秩块，树深、块大小、 admissibility 参数都直接暴露给压缩器。
- 因此它的“树”不是围绕某个特定 Maxwell transmission block 专门设计，而是围绕矩阵块压缩与 ACA/低秩逼近统一构建。

### 矩阵向量乘积如何处理

- H-matrix matvec 的本质是：
	- 近场块直接稠密乘法；
	- 可容许远场块用低秩因子乘法；
	- 整体由 cluster-tree 管理块遍历。
- 这和 MLFMA 的 aggregation -> translation -> disaggregation 链路完全不同。FreeFEM 的重点不是球谐/平面波转移，而是 admissible block 的低秩近似质量。
- 因而它更自然地把 `eps`、`eta`、树深、块阈值纳入“压缩误差预算”，这也是公开实现里 FGMRES 与 H-matrix 预条件能并列出现的原因。

### 反映出的工程策略

- 公开实现会把“压缩误差控制”作为 solver 主线的一部分，而不是只把快速算法当黑盒 matvec。
- 对复杂 BEM 问题，他们愿意用 FGMRES 这类更稳健的 Krylov 包装去适配非平凡预条件和近似算子。

### 对 EMSuite 的启示

- 当前 EMSuite 的问题本质上也不是“GMRES 不会用”，而是缺少对 fast operator 误差预算的正式控制接口。
- MLFMA 若继续作为主线，必须增加“压缩误差对求解目标精度的约束”这一层治理，至少要有面向求解容差的 `L / near-far / translation accuracy` 预算，而不是只做固定经验值。
- 如果 Phase 15 的目标是可靠高精度介质结果，那么 H-matrix 类替代路线应被视为验证/兜底基线，而不是无关路线。

## 4.3 scuff-em

### 实现特征

- 明确实现 PMCHW 与 N-Muller 两条 dielectric formulation。
- 文档层面对 dielectric object 单列 PMCHWT 与 N-Muller 两种连续方程。
- 在公开 scattering/transmission 应用中，常见主线是组装系统矩阵后直接 `LUFactorize()` / `LUSolve()`。

### 反映出的工程策略

- 对 dielectric scattering，scuff-em 的公开路径更强调 formulation 选择和稳健 direct solve，而不是证明“快速迭代一定够用”。
- 这恰好说明当前 EMSuite 遇到的问题并不反常。公开成熟代码也没有把“PMCHW + fast iterative”当成天然稳定的默认配置。

### 对 EMSuite 的启示

- 若 Phase 15 的优先级是 correctness 与可验收精度，dense/direct 与 alternative formulation 应继续保留为第一层基准。
- 不应为了“统一使用 MLFMA + GMRES”而牺牲 transmission 问题的高精度可控性。

## 4.4 与 EMSuite 当前实现的直接差异

### EMSuite 当前本地结构

- EMSuite 当前 `PMCHWMLFMAOperator` 是显式 PMCHW 专用实现：
  - 自己持有 `octree0` 与 `octree1` 两棵 N 点八叉树；
  - `leaf_size` 先被统一收敛到 `λ_min/10` 标尺；
  - `near_range` 再按该标尺显式放大并传入本地 octree builder；
  - `mul!` 明确执行 `Z_near*x + y_far`，其中 `y_far` 由 `J×k0`、`J×k1`、`M×k0`、`M×k1` 四遍手工累加得到。
- 本地 translation 也由 EMSuite 自己预计算，`near_range` 会直接改变 far-offset 枚举与 `αTrans` 索引表。
- 这意味着本地 PMCHW 的树结构、近远场划分、translation stencil、四遍调度，全都直接属于 PMCHWMLFMAOperator 的实现责任。

### 与 Bempp 的关键不同

1. 树的归属不同。
	- EMSuite: PMCHW 算子自己拥有并调度两棵八叉树。
	- Bempp: transmission 仍是 block operator，树/FMM 由通用 ExaFMM evaluator 负责。
2. matvec 分解层次不同。
	- EMSuite: 先按 PMCHW 物理块手写四遍远场，再在每遍内部做 aggregation/translation/disaggregation。
	- Bempp: 先保持 electric/magnetic block operator 语义，再让每个 block 的非局部部分调用 FMM evaluator，最后在 block algebra 层组合。
3. 近场修正的组织方式不同。
	- EMSuite: 一个整体 `Z_near` 稀疏 2N×2N 矩阵与四遍远场叠加。
	- Bempp: 每个 operator block 都显式保留 `singular_part`，形成“非局部 evaluator + 奇异稀疏修正”的结构。
4. 参数控制粒度不同。
	- EMSuite: 当前主要通过 `leaf_size`、`near_range`、本地 translation 参数间接控制误差。
	- Bempp: 公开参数直指 `depth`、`expansion_order`、`ncrit`、`near_field_representation`，即直接把 FMM 误差治理接口暴露出来。

### 与 FreeFEM/Htool 的关键不同

1. 树结构不同。
	- EMSuite 是八叉树 + translation stencil。
	- FreeFEM/Htool 是 cluster tree + admissible blocks。
2. 远场近似机制不同。
	- EMSuite 使用 MLFMA 的多层聚合/平移/解聚。
	- FreeFEM/Htool 使用 H-matrix 的低秩块因子化。
3. 误差预算接口不同。
	- EMSuite 当前仍偏向“固定经验参数后直接求解”。
	- FreeFEM/Htool 明确把压缩精度、块阈值、树深并入求解控制面板。

### 对当前长 Krylov 问题的具体启示

- Bempp 式结构提示：应把 PMCHW 的快速乘法进一步拆回到 block/operator 层，而不是长期把所有误差都压在一个 PMCHW 专用四遍 `mul!` 里排查。
- Bempp 的 `singular_part + FMM evaluator` 结构说明，EMSuite 后续应分别量化：
  - 近场稀疏块误差；
  - 每个 Maxwell 子块的远场 evaluator 误差；
  - block 组合后对 Krylov 子空间方向的放大效应。
- FreeFEM/Htool 式结构提示：即便坚持 MLFMA，也必须把 `leaf_size / near_range / translation order` 从固定经验值提升为“面向目标残差与目标阻抗误差”的正式预算接口。
- scuff-em 的路线则再次提醒：对 dielectric transmission，alternative formulation 与 dense/direct 基线不是过渡品，而是必须长期保留的真值约束。

## 5. 对比结论

| 系统 | transmission 组织方式 | 树/层次结构 | matvec 主体 | 预条件/正则化重点 | 对 EMSuite 的结论 |
|------|----------------------|-------------|-------------|-------------------|-------------------|
| EMSuite 当前 | 顶层已迁移为 PMCHW block/operator shell，MLFMA 仍保留本地专用 backend 内核 | backend 内部仍为两棵显式八叉树 `octree0/octree1` | shell 负责 block algebra；MLFMA backend 仍是 `Z_near*x +` 四遍远场手工调度 | strong-form 已有正式入口，但仍缺压缩误差预算接口 | 长 Krylov 偏差已被收束到 backend fidelity / formulation regularization 边界 |
| Bempp | multitrace / Calderon 风格 block operator | 通用 ExaFMM 树，参数为 `depth/ncrit/order` | 每个 Maxwell block = `FMM evaluator + singular_part`，再做 block algebra | strong form, Calderon, dual pairing | 优先补 operator-level 路线 |
| FreeFEM+BEMTool+Htool | H-matrix BEM | cluster tree + admissible blocks | 近场稠密块 + 远场低秩块 | GMRES/FGMRES + H-matrix 误差控制 | 需要正式误差预算，不是只调 restart |
| scuff-em | PMCHWT + N-Muller | 公开主线不强调 PMCHW MLFMA 树 | 多为组装后 direct solve | formulation 选择优先于裸迭代 | EMSuite 不应只押注 PMCHW + MLFMA |

## 6. 对 EMSuite 的建议

### 6.0 架构路线决策

建议采用“**先迁移到 Bempp 式 block/operator 架构，再补 FreeFEM 风格 H-matrix 后端**”的路线，但要明确迁移目标是**组织形式与职责分层**，不是一次性照搬 Bempp 的具体实现或外部依赖。

更具体地说，当前最值得迁移的是下面三点：

1. 保留 PMCHW 的物理 block 语义。
	- 不再把 transmission 主线长期收束为一个专用 `2N` 单体 `mul!`。
	- 应把 `EJ/EM/HJ/HM` 四块重新提升为一等公民。
2. 把“非局部快速乘法”降为 backend 责任。
	- operator algebra 负责 block 组合、weak/strong form、mass matrix preconditioning。
	- MLFMA 只负责某个 block 的 nonlocal evaluator。
3. 把 FreeFEM 的 H-matrix 路线放到下一阶段。
	- 在 block/operator 架构还没稳定之前就并入 H-matrix，只会把 formulation、operator、compression 三类问题再次耦合在一起。

因此，若问题是“现在要不要直接跳去 H-matrix”，我的建议是否定的；若问题是“现在要不要先把本地架构整理成 Bempp 那种 block/operator 形态”，我的建议是肯定的。

### 6.1 立即执行

1. 把当前问题正式归类为“formulation + operator-accuracy”双因素问题，不再归类为“缺预条件器”。
2. 已完成：PMCHW 主线已迁移成 Bempp 风格的 block/operator 组织，四块算子与显式 `weak_form/strong_form` 已落在 shell 层。
3. 已完成：PMCHW strong-form / mass-matrix-aware 对照主线已落地，且中尺度四路 Gate S 已正式执行。
4. 把 medium fixture 上的 Krylov 诊断从“随机向量 Gate C”升级为“Arnoldi / 迭代子空间方向上的 matvec 误差诊断”。
5. 在 PMCHW 路径内再做一次算子拆分：至少分别度量 `E` block 与 `H` block 的 fast matvec 误差，而不是只看最终 2N 合成误差。

### 6.2 优先级最高的新对照

1. 已完成：dense 路径上的 PMCHW block operator 显式化已经落地，并验证了 block 组合与现有 Direct PMCHW 一致。
2. 已完成：在 block/operator 架构上，dense 与 MLFMA 的相同强形式 Gate S 对照已经具备正式中尺度回归。
3. 已完成：N-Muller formulation 的小球 dense 对照，当前已正式证明它在共享球夹具上优于 PMCHW 的 conditioning 与 GMRES 行为。
4. 下一优先项应转到更高阶 medium-scale dielectric sphere，对 PMCHW / N-Muller 继续比较 conditioning、输入阻抗与 GMRES 收敛行为是否保持同一趋势；当前 `N=150` 预设已给出正向证据，但还不是正式 medium gate。
5. 为 MLFMA backend 增加“目标求解容差 -> 压缩参数”的误差预算诊断，而不是固定 `near_range` 和 leaf 标尺后直接求解。
6. 把当前 B2 输入阻抗基准从“已可执行”继续推进到“参数预算敏感度可记录”，避免它只给出单点通过结论。

### 6.3 中期路线

1. 若 PMCHW strong-form 仍无法收回长 Krylov 偏差，则主线应继续向 Bempp 的 Calderon / multitrace 风格 regularization 演进，而不是继续试 `ILU`、`SPAI`、`BlockJacobi` 的变体。
2. 只有在 block/operator 架构与 strong-form 分界都稳定后，才引入 FreeFEM 风格 H-matrix backend；此时 H-matrix 应作为与 MLFMA 平行的 backend，而不是新的单体求解器分支。
3. 若项目需要一个更稳健的中尺度高精度介质解法，应优先让 `DenseBackend / MLFMABackend / HMatrixBackend` 共享同一 block/operator 外壳，以便做真正可比的 backend 级误差评估。

### 6.4 推荐外部基线用例库

下面这些不是“阅读材料”，而是适合作为 EMSuite 精度对照的外部基线来源。选择原则是：

- 几何简单，可在 EMSuite 现有网格/后处理能力上复刻；
- 观测量明确，不依赖对方内部私有状态；
- 能分别覆盖 dense operator、compressed matvec、far-field、solver、transmission 五个层次。

| 来源 | 开源用例 | 公开观测量 | 适合映射到 EMSuite 的哪一层 | 建议优先级 |
|------|----------|-----------|-----------------------------|------------|
| Bempp | `test/validation/operators/boundary/test_maxwell_boundary.py` 中 sphere 上的 electric/magnetic boundary operator | dense 离散矩阵或 dense matvec | Maxwell 子块离散正确性；适合做 `E/H` block 级 dense 基线 | 高 |
| Bempp | `test/validation/fmm/test_fmm.py` 中 Maxwell boundary FMM，`regular_sphere(2)`、`k=1.5` 与 `k=1.5+0.5i` | `dense @ vec` 与 `fmm @ vec` 对齐，`rtol=2e-3` | compressed matvec 容差基线；适合给 EMSuite 的 fast-vs-dense block 误差定量找量级 | 最高 |
| Bempp | `test/validation/fmm/test_fmm.py` 中 Maxwell potential FMM | 指定点集上的 electric/magnetic potential | 远场/势场后处理正确性；适合约束 EMSuite far-field/potential evaluator | 中 |
| Bempp | `examples/maxwell/maxwell_dielectric.py` | 两个半径 `0.4` 球体、`300 MHz`、`eps_r=2.1` 的 transmission 解与散射截面曲线 | multi-body dielectric direct baseline；适合做 PMCHW transmission 结果的几何扩展对照 | 中 |
| Bempp | `bempp_cl/api/linalg/iterative_solvers.py` + handbook `use_strong_form=true` | weak/strong form 使用同一 GMRES 接口、同一 RHS 的迭代数与残差曲线 | Gate S 的方法学基线；适合约束 EMSuite strong-form 对照流程 | 最高 |
| FreeFEM | `examples/bem/Maxwell_EFIE_sphere.edp` 与文档中的 `radius=1, 600 MHz, RT0S, eta=10, eps=1e-3` | H-matrix EFIE 解与重建场 | 压缩误差预算基线；虽非 PMCHW，但适合给“压缩参数必须进入误差预算”提供外部参照 | 中 |
| scuff-em | `doc/docs/tests/MieScattering/MieScattering.md` | 球散射与解析 Mie 级数的 cross section / PFT 对照 | dielectric sphere 真值基线；最适合做 EMSuite 单球 transmission 的第一层正确性约束 | 最高 |
| scuff-em | `unitTests/unit-test-BEMMatrix.cc` | `Dielectric sphere, low/medium real frequency`, `imag frequency`, `Two dielectric spheres` 的矩阵回归 | dense dielectric assembly 回归；适合作为 PMCHW / N-Muller dense 组装的外部结构性基线 | 高 |
| scuff-em | `doc/docs/applications/scuff-tmatrix/scuff-tmatrix.md` 中 `E10Sphere_327` / `E10Sphere_1362` | 介质球 T-matrix 元与解析球解对照 | 模态散射基线；适合作为更高精度球散射校核而不是第一优先实现 | 中 |

对 EMSuite 最有直接价值的，不是全部一起做，而是分三批落地：

1. 第一批，先锁真值与 fast-vs-dense 量级。
	- scuff-em 单球 Mie scattering。
	- Bempp Maxwell boundary FMM sphere real/complex wavenumber。
2. 第二批，锁 solver/formulation 分界。
	- Bempp weak/strong GMRES 对照流程。
	- EMSuite 自己的 dense weak/strong 与 MLFMA weak/strong 四路对照。
3. 第三批，再扩展到更复杂几何与替代压缩路线。
	- Bempp 双介质球 transmission。
	- FreeFEM Maxwell EFIE sphere H-matrix。
	- scuff-em dielectric sphere / two-sphere matrix regression。

这些外部基线在 EMSuite 中的推荐对应关系如下：

| EMSuite 目标 | 最优先外部基线 | 原因 |
|-------------|----------------|------|
| 单球 dielectric scattering 真值 | scuff-em Mie scattering | 有解析级数对照，最容易把 PMCHW 正确性与 solver 误差分开 |
| fast matvec 误差量级 | Bempp Maxwell boundary FMM | 公开给出了 dense vs compressed 的直接 matvec 对照与容差 |
| strong-form 方法学 | Bempp GMRES `use_strong_form` | weak/strong 切换接口清晰，最适合作为 Gate S 的外部范式 |
| 压缩误差预算思路 | FreeFEM Maxwell EFIE sphere | 明确暴露 `eta/eps/minclustersize` 等控制量，适合对照 EMSuite 未来预算接口 |
| alternative formulation 压力测试 | scuff-em PMCHWT / N-Muller 文档与矩阵回归 | 可用于证明 PMCHW 不是唯一高精度主线 |

## 7. 不建议继续投入的方向

- 不建议继续把主要精力放在 `Diagonal`、`ILU(Z_near)`、leaf Block Jacobi 这类局部代数预条件器上，除非先有文献或谱分析证据支持。
- 不建议再用“扩大近场直到结果变好”作为主线，因为当前夹具已证明该路线会撞到内存墙。
- 不建议再把随机向量 Gate C 绿色解读为“高精度求解已无算子问题”。

## 8. 推荐执行顺序

1. 已完成：Bempp 风格 PMCHW block/operator 外壳。
2. 已完成：Dense backend 接入该外壳，并完成 weak/strong form 对照。
3. 已完成：现有 PMCHW MLFMA 已收编为 backend，而不是继续维持单体主线。
4. 已完成：在统一外壳下比较 dense 与 MLFMA 的 strong-form 行为。
5. 已完成：Dense N-Muller 基线实现与小球球算例验证。
	- 对应执行文档：`.github/plans/phase_15_nmuller_dense_baseline.md`。
	- 当前正式对照入口：`test/test_nmuller.jl`、`test/test_nmuller_comparison.jl`、`benchmark/compare_pmchw_nmuller_sphere.jl`。
6. 下一步：把 PMCHW / N-Muller 对照扩展到 medium-scale sphere，确认小球上观察到的 formulation 优势是否在更大离散系统中保持。
7. 若 block/operator + strong-form + N-Muller dense 对照后仍失败，再进入 Calderon / multitrace regularization 设计。
8. 在有了误差预算接口之后，再补 FreeFEM 风格 H-matrix backend。

## 9. 最终结论

外部文献与开源实现给出的共识非常一致：

- Maxwell transmission / dielectric scattering 的难点首先是 formulation 与 operator conditioning，其次才是线性代数层面的 Krylov 加速。
- 公开成熟实现普遍不会把“PMCHW + 近场块预条件 + 裸 GMRES”当成最终高精度路线。
- EMSuite 已经沿着这条路线完成了第一轮关键迁移：shell、strong-form、MLFMA backend 收编与中尺度四路 Gate S 都已落地。
- 对 EMSuite 当前问题，最合理的下一步不再是“是否要做 shell/strong-form”，也不再是“是否需要先补 N-Muller baseline”；这两步现在都已完成。
- 当前真正尚未关闭的是两条边界：
	1. formulation 侧：把已在小球上观察到的 N-Muller 优势扩展到更大离散系统，并继续绑定外部真值；
	2. backend 侧：把 MLFMA 的压缩误差控制从固定经验参数提升为正式预算接口。

## 10. 参考资料

- Bempp publications: https://bempp.com/publications.html
- Bempp handbook solvers: https://bempp.com/handbook/api/solvers.html
- Bempp Maxwell dielectric example: https://github.com/bempp/bempp-cl/blob/main/examples/maxwell/maxwell_dielectric.py
- Bempp Maxwell operators: https://github.com/bempp/bempp-cl/blob/main/bempp_cl/api/operators/boundary/maxwell.py
- FreeFEM BEM documentation: https://doc.freefem.org/documentation/BEM.html
- FreeFEM Maxwell BEM section: https://doc.freefem.org/documentation/BEM.html#maxwell-bem-problem-in-freefem
- FreeFEM BEM/H-matrix implementation: https://github.com/FreeFem/FreeFem-sources/tree/main/plugin/mpi/bem.hpp
- scuff-em implementation notes: https://github.com/HomerReid/scuff-em/tree/main/doc/docs/forDevelopers/Implementation.md
- scuff-em dielectric formulations reference: https://github.com/HomerReid/scuff-em/tree/main/doc/docs/tex/Formulations.tex
- Scroggs, Betcke, Burman, Smigaj, van 't Wout, Software frameworks for integral equations in electromagnetic scattering based on Calderon identities, arXiv:1703.10900, https://arxiv.org/abs/1703.10900