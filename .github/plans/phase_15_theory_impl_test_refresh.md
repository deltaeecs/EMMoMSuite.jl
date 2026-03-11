# Phase 15 Theory-Implementation-Test 治理闭环

> 目标: 对 PMCHW transmission 主线建立“先理论、后实现、再验收”的强约束开发合同；当前架构目标已从单体 `PMCHWMLFMAOperator` 刷新为 **Bempp 风格 block/operator 外壳 + 可插拔 backend**。  
> 状态: 生效中。未经本文档定义的理论前提、函数合同、门禁测试，不允许继续实现或修改。  
> 通过条件: 文档检视连续 3 轮无新发现；正式 PMCHW 主线 Gate A/B/C/D/S/O 全部通过。共享 MLFMA core 的诊断性 `@test_broken` 可保留，但必须被明确标注为“机制回归”而非主线阻塞。

## 当前状态

- 正式 PMCHW 主线已恢复 green：`GD0`、`GD1`、`GD2`、`GB`、`Gate C`、`B2` 当前均已通过。
- 当前保留红灯仅用于锁定 shared core 机制：较小 `near_range` 下，leaf far-path 会把过近交互误留在远场。
- 因此后续 Phase 15 工作重点不再是“修正式 PMCHW operator 正确性”，而是“在保持正式 green 的前提下回收 near-field 膨胀带来的性能代价”。
- 另一个已确认的后续问题是**中尺度求解阶段**：当 `N=540` 时，`dense+GMRES(200)` 与 `MLFMA+GMRES(200)` 一致，但两者对 `LU` 仍有明显阻抗偏差，因此不能再把该问题简单记为“MLFMA matvec 正确性已完成”。
- 最新修正：上述现象已被重新定义为**谱条件问题与 MLFMA 压缩误差的耦合问题**。长 Krylov 实验表明，dense GMRES 在大 `restart` 下已基本恢复 `LU`，而 PMCHW MLFMA 仍保留可观阻抗偏差，因此本 Phase 必须新增一条 **strong-form 对照主线**，把 conditioning 与 fast-operator fidelity 显式拆开。
- 当前架构路线也已刷新：后续主线不再建议继续堆叠 `PMCHWMLFMAOperator` 单体特化逻辑，而应先迁移到 **Bempp 风格 block/operator 外壳**，把 `EJ/EM/HJ/HM` 四块重新提升为一等公民，再把 Dense、MLFMA、未来 H-matrix 作为 backend 并行接入。
- 同时，文档示例使用的 `BlockJacobiPreconditioner(op, basis)` 接口已正式补齐；但在当前 `N=540` 夹具上，`Diagonal` / `ILU(Z_near)` / `BlockJacobi(op)` 均未改善 PMCHW 收敛，因此后续工作不得再默认“补一个预条件器接口即可解决中尺度偏差”。

## 0. 适用原则

- 原则 1: 严格 TDD，顺序必须是 RED -> GREEN -> REFACTOR。
- 原则 2: 可对齐路径必须以 Legacy 为唯一真相源；PMCHW 特有路径以 EMSuite Direct PMCHW 为内部真相源。
- 原则 3: 排查顺序固定为几何 -> 常数 -> 积分 -> 组装 -> MLFMA 四步链路。
- 原则 5: 每个逻辑闭环验证成功后立即提交，禁止积压大批未提交试验。
- 原则 7: Phase 结束前必须做检视迭代；本 Phase 改为连续 3 轮无新发现才允许通过。
- 原则 9: 计划文档必须精确到函数、测试、门禁、完成定义。

## 1. 为什么之前会失控

之前的问题不是“没有写文档”，而是文档不可执行：

- 文档只描述现象，没有把函数输入、输出、理论前提、验收条件绑定到一起。
- 对 near_range、leaf_size、L_min、kernel medium convention 的假设没有冻结，导致修改跨层扩散。
- 调试脚本和正式 mul! 流程不一致，结论不具备可追溯性。
- Gate B 只看结果，没有拆到 pass 级别的必要中间量，导致错误归因依然靠猜。
- 检视迭代没有闭环标准，出现“写了计划但没有证据表明计划已收敛”。

## 2. 理论基线与真相源

### 2.1 Legacy 基线

- MoM_Basics: 几何、基函数、符号约定、局部边/面方向。
- MoM_Kernels: Green 函数、奇异积分、L/K 算子核、常数链。
- MoM_AllinOne: 端到端矩阵组装与求解流程。

### 2.2 PMCHW 特有基线

Legacy 中没有完整 PMCHW MLFMA 实现，因此 PMCHW 的真相源分两层：

- 块矩阵真相源: EMSuite Direct PMCHW 组装结果。
- MLFMA 真相源: Direct PMCHW 的 far-field block decomposition。

### 2.3 冻结的统一约定

未经单独更新本文档，以下约定不得擅自改变：

- 相位约定: 源端聚合使用 $e^{+jk \hat{r} \cdot (r-r_c)}$，接收端解聚使用 $e^{-jk \hat{r} \cdot (r-r_c)}$。
- Translation 定义: 相对位移固定为 target cube ID3D 减 source cube ID3D。
- near/far 互补定义: 叶层 near 由 searchNearCubes 决定；far 必须是其严格补集。
- kernel medium convention: Direct 若采用某一 $k$ 与 $\eta$ 的近似，MLFMA pass 必须使用同一 convention，不允许一边 real-k 一边 complex-k。
- 诊断流一致性: 任何用于判断正确性的脚本，其聚合/上行/平移/下行/解聚顺序必须与正式 mul! 完全一致。

### 2.4 冻结的架构目标

未经单独更新本文档，后续 PMCHW transmission 主线的架构目标固定为：

- 顶层不再以单个 `PMCHWMLFMAOperator` 作为唯一演化中心，而是以 **PMCHW block/operator 外壳** 为中心。
- 顶层外壳必须显式保留四块语义：`EJ`、`EM`、`HJ`、`HM`。
- `weak_form`、`strong_form`、mass-matrix-aware solve、block algebra 属于外壳责任。
- Dense、MLFMA、H-matrix 仅作为 backend；backend 不得反向定义 transmission 的物理分块语义。
- 现有 `PMCHWMLFMAOperator` 允许作为兼容 facade 保留，但不得再承担新的主线架构责任。

## 3. 开发对象分层

本 Phase 的核心不是“一个算子”，而是一串分层函数。任何问题必须先定位到层，再定位到函数。

### 3.1 几何与八叉树层

- build_octree
- setLevelInfo! 两个重载
- searchNearCubes
- setKidLevelFarNeighbors!

### 3.2 插值与平移预计算层

- truncationLCal
- levelIntegralInfoCal
- compute_translation_factors!
- cal_alpha_trans_on_level!

### 3.3 远场链路层

- aggregate_leaf_pmchw!
- aggregate_upward!
- translate!
- disaggregate_downward!
- disaggregate_leaf_pmchw_j!
- disaggregate_leaf_pmchw_m!

### 3.4 近场与总算子层

- assemble_near_field_pmchw
- PMCHWMLFMAOperator constructor
- mul!

### 3.5 Operator Algebra 层

- PMCHW block/operator shell
- EJ/EM/HJ/HM block composition
- weak_form / strong_form
- mass matrix / dual pairing

### 3.6 Backend 层

- Dense backend
- MLFMA backend
- 未来 H-matrix backend

### 3.7 验收层

- test/test_pmchw_mlfma_operator.jl
- 共享 MLFMA smoke tests

### 3.8 Strong-Form 对照层

- Dense weak-form solve
- Dense strong-form solve
- MLFMA weak-form solve
- MLFMA strong-form solve

## 4. 函数级开发合同

本节是核心。每个函数必须回答四个问题：理论职责、禁止猜测的点、开发步骤、验证方式。

### 4.1 build_octree

理论职责:
- 给定 basis center、leafCubeEdgel、介质波长参数，构造多层 cube 拓扑。
- 只负责拓扑与层级信息，不允许偷偷引入“精度修正策略”。

禁止猜测:
- 不允许因为 Gate B/C 失败，直接在 build_octree 内拍脑袋修改 near_range 或 leaf_size_eff。
- 如果需要自适应参数，必须先在本文档新增理论推导和门禁矩阵。

开发步骤:
1. 固定输入参数语义: leafCubeEdgel 是叶层几何长度，不是“目标精度旋钮”。
2. 明确每一层 cubeEdgel、nLevels、kidsInterval、bfInterval 的构造公式。
3. 明确 build_octree 对 setLevelInfo! 与 compute_translation_factors! 的参数透传关系。

补充说明:
- 对常规 PEC MLFMA，lowest-level leaf size 取约 `0.2-0.25 λ` 是常见经验值。
- 本 Phase 的 PMCHW 双介质路径不能机械套用该经验值；这里的 `λ` 必须理解为最短介质波长，而且还必须与 `near_range` 联动验证。
- 若仅凭“0.25 λ 经验值”压低 leaf 尺寸或 near span，而未通过 Gate D2/GB/C/B2 复验，则视为违反开发合同。

验证方式:
- 单测: cube center 覆盖完整、kidsInterval 连续、bfInterval 无重叠。
- Gate A: 两棵 octree 在要求一致的拓扑前提下，sorted_ids 与 inv_sorted_ids 的语义必须可证明。

### 4.2 setLevelInfo! 两个重载

理论职责:
- 叶层版本负责 basis center 到 cube 的分桶。
- 非叶层版本负责 child cube 到 parent cube 的聚合。

禁止猜测:
- 不允许把某层的 L、near_range、neighbors 修正写成“只对某个 level 生效”的隐式分支，除非文档明示。

开发步骤:
1. 先证明 leaf 与 non-leaf 的 cubesID3D 生成规则一致。
2. 再接入 searchNearCubes 结果。
3. 最后接入 levelIntegralInfoCal 生成 L/poles。

验证方式:
- 单测: 每层 nCubes、ID3D 唯一性、kidsInterval 闭区间连续性。
- 诊断: parent.kidsInterval 展开后的 child 集与实际 child cube 集完全一致。

### 4.3 searchNearCubes

理论职责:
- 定义某层 cube 的 near 邻域，这是 near/far 划分唯一入口。

禁止猜测:
- 不允许为了“让 Translation 更稳定”临时改 near_range，而不同时证明对 leaf-level near matrix 和 non-leaf farneighbors 的影响。

开发步骤:
1. 固定 near 的数学定义: 三维 offset 立方体内的所有 cube。
2. 明确边界裁剪规则。
3. 明确输出 neighbors 包含自身。

验证方式:
- 单测: 内部 cube、边界 cube、角点 cube 三类用例。
- 不变量: far = 全部同层 cube 减 near。

### 4.4 setKidLevelFarNeighbors!

理论职责:
- 把 parent near 关系映射为 kid level farneighbors，是多层 MLFMA 的核心定义转换。

禁止猜测:
- 不允许只看 leaf far pair 数量就判断正确；必须验证“覆盖无缺口、无重叠”。

开发步骤:
1. 先证明 neighborsKids = parent near neighbors 的孩子全集。
2. 再证明 kidCube.farneighbors = neighborsKids \ kidCube.neighbors。

验证方式:
- Gate A: 对任意 kidCube，nearneighbors 与 farneighbors 互斥且并集等于 parent near neighborhood 的 child 展开。

### 4.5 truncationLCal / levelIntegralInfoCal

理论职责:
- 给定 cubeEdgel 和波长，生成截断数 L 与极点集合。

禁止猜测:
- 不允许用经验常数强行把 L 往上调，除非先证明当前 L 违反截断准则并把新准则写入文档。

开发步骤:
1. 冻结 truncationLCal 公式。
2. 明确 L_min 仅是下界，不是替代理论公式。
3. 明确 poles 与 quadrature 权重和必须满足的积分精度目标。

验证方式:
- 单测: cubeEdgel 增大时 L 单调不减。
- 单测: sum(Wθϕs) 逼近 $4\pi$。
- 诊断: 记录每层 $kR_{min}/L$，只作为分析指标，不作为自动调参触发器。

### 4.6 compute_translation_factors! / cal_alpha_trans_on_level!

理论职责:
- 为某一层所有 far offset 预计算 translation 因子 $\alpha$。

禁止猜测:
- 不允许在这里“补偿”别处的错误，比如通过改常数因子去掩盖聚合或解聚问题。

开发步骤:
1. 先固定 far offset 索引范围与 αTransIndex 的定义。
2. 再验证单个 offset 的解析级数与数值积分关系。
3. 最后验证 αTrans 的列和与目标 Green 函数之间的理论对应。

验证方式:
- 单测: αTransIndex 对 near offset 不给出 far 有效索引。
- Gate B 前置: 单 offset translation 精度测试，固定 $R$、固定 $L$、固定 pole set。
- 诊断: 对 level、offset、kRab 三元组输出误差，不允许只看 norm。

### 4.7 translate!

理论职责:
- 对每个 target cube，从其 farneighbors 汇总所有 source cube 的 aggS 贡献。

禁止猜测:
- 不允许因为输出偏小/偏大就直接怀疑 αTrans；先核对 farneighbors 覆盖，再核对 aggS 源，再核对 factor。

开发步骤:
1. 明确 disaggG 在每一轮 translate! 前必须清零。
2. 明确 relative3DID 的方向是 target-source。
3. 明确每个 farneighbor 的贡献是逐极点逐极化累加。

验证方式:
- 单测: 单 source cube、单 farneighbor 场景下，disaggG 等于 factor .* aggS。
- 诊断: worst column 对应 cube 的 farneighbor 贡献排序必须可导出。

### 4.8 aggregate_leaf_pmchw!

理论职责:
- 将 PMCHW 的 J 或 M 未知量在叶层投影到球面极点场。

禁止猜测:
- 不允许把 J/M、k0/k1、E/H 接收逻辑混到一个“先跑通再说”的大函数里。

开发步骤:
1. 先实现 J-pass，固定 x_range 与 kernel convention。
2. 通过 Gate B 的 EJ-k0、EJ-k1 后，再实现 M-pass。
3. 明确每一遍开始前 aggS 必须清零。

验证方式:
- 单测: x 仅一个 basis 激励时，leaf aggS 非零 cube 只出现在该 basis 所在 cube。
- Gate B: pass-level 对齐必须分 J/k0、J/k1、M/k0、M/k1 四遍单独验收。

### 4.9 disaggregate_leaf_pmchw_j! / disaggregate_leaf_pmchw_m!

理论职责:
- 将叶层 disaggG 回投到 RWG testing function，形成 EJ/HJ/EM/HM 四块输出。

禁止猜测:
- 不允许通过改末端 block scaling 去掩盖前面 pass 的错误。

开发步骤:
1. 先把 J receive kernel 和 M receive kernel 理论分清。
2. 再分别实现 J receive 与 M receive，不允许共用含糊中间式。
3. 最后核对与 PMCHW Direct block 的外部常数链。

验证方式:
- 单测: 只激励 J-pass 时，HJ/EM 的符号关系必须满足理论不变量。
- Gate A: near-field block 中 HJ + EM 的抵消关系。
- Gate B: 每块输出独立和 direct far block 比较。

### 4.10 assemble_near_field_pmchw

理论职责:
- 以 leaf-level near pattern 为掩膜，直接精确组装 2N x 2N 的 EJ/EM/HJ/HM 四块近场矩阵。

禁止猜测:
- 不允许为了让总误差变小而擅自扩大近场掩膜；任何 near pattern 修改都必须先在 searchNearCubes 层被证明。

开发步骤:
1. 先验证稀疏模式来源唯一且来自 leaf near。
2. 再分别验证四块填充位置和符号。
3. 先在 Gate D0 重新确认 dense PMCHW 本身对 Mie 基线成立。
4. 最后验证与 dense PMCHW 在同一 sparsity pattern 上逐元素一致。

验证方式:
- Gate D0: Dense PMCHW 球散射必须先对 Mie 基线复验通过，E-plane RMSE < 1.5 dB，H-plane RMSE < 2.0 dB，后向散射误差 < 0.2 dB。
- Gate D1: 在 `Z_near` 的非零位置上逐元素比较 dense PMCHW，最大元素相对误差 < 0.1%。
- Gate A: HJ + EM near-field block infinity norm ratio < 1e-8。

### 4.11 PMCHWMLFMAOperator constructor

理论职责:
- 只负责把 Direct PMCHW 所需的 near field、octree、permutation、共享缓冲区组装起来。

禁止猜测:
- 不允许在 constructor 里夹带“临时 debug 逻辑”或“只对某 εr 生效”的试验参数。

开发步骤:
1. 先完成纯构造，不加自动调优。
2. 确认 sorted_ids0 / sorted_ids1 / inv_sorted_ids0 / inv_sorted_ids1 语义一致。
3. 确认 Z_near 的索引空间和 mul! 的输入空间一致。

验证方式:
- 单测: permutation 为双射。
- Gate A: non-trivial near/far split，且近场矩阵维度与 2N 一致。

### 4.12 mul!

理论职责:
- 组织完整四遍远场和近场叠加，输出最终 2N 向量。

禁止猜测:
- 不允许在 mul! 里加“经验校正因子”“临时比例系数”“只修某一块的补偿项”。

开发步骤:
1. 先写出 pass 顺序与每遍输入输出缓冲区定义。
2. 再分别实现四遍链路，不做合并简写。
3. 最后做近场叠加和 permutation 还原。

验证方式:
- Gate B: 四遍分别与对应 direct far block 对齐。
- Gate C: `15.11 MLFMA mul! vs Direct` 相对误差 < 0.10。
- Gate D: 共享 MLFMA 核心的 EFIE/SCFIE smoke tests 仍通过。

### 4.13 strong-form 对照主线

理论职责:
- 用统一右端、统一 Krylov 参数、统一验收指标，比较 weak-form 与 strong-form 两类离散系统，显式分离“谱条件恶化”与“MLFMA 压缩误差”。

禁止猜测:
- 不允许只因为 `strong-form` 迭代次数更少，就直接宣布 MLFMA fidelity 问题已解决。
- 不允许在 strong-form 对照中更换 RHS、容差、restart、停止准则，从而把 solver 参数差异伪装成 formulation 结论。
- 不允许跳过 dense strong-form 基线，直接只做 MLFMA strong-form。

开发步骤:
1. 先定义 PMCHW 的 strong-form/ mass-matrix-aware 离散系统，并冻结 weak/strong 的未知量与右端映射关系。
2. 先在 dense 路径上做 weak vs strong 对照，比较与 `LU` 的阻抗误差、解向量差异、残差历史。
3. 再在 MLFMA 路径上做 weak vs strong 对照，保持同一 Krylov 参数与停止准则。
4. 最后比较 `dense strong` 与 `MLFMA strong` 的差异，把剩余偏差归入 fast-operator fidelity，而不是继续混在 conditioning 里。

验证方式:
- Gate S1: `dense weak` 与 `dense strong` 必须对同一 direct/LU 真值给出可比较的阻抗与残差曲线。
- Gate S2: 若 `dense strong` 明显优于 `dense weak`，则正式记录“conditioning 是显著因素”；若两者接近，则不得把主因归到 strong-form 缺失。
- Gate S3: 在 strong-form 下比较 `dense` 与 `MLFMA`；若 `dense strong` 已接近 `LU` 而 `MLFMA strong` 仍有显著阻抗偏差，则正式记录“剩余主因属于 MLFMA 压缩/算子保真度”。
- Gate S4: 记录指标必须至少包括输入阻抗误差、相对残差、迭代步数、相对解误差或 Krylov 子空间方向误差。

基线用例约束:
- Gate S 不允许只依赖 EMSuite 自身历史结果，必须至少绑定两类外部基线：
	- 单球 dielectric 真值基线：优先采用 scuff-em 的 Mie scattering / dielectric sphere 系列；
	- compressed-vs-dense matvec 基线：优先采用 Bempp Maxwell boundary FMM sphere real/complex wavenumber 系列。
- 若做 strong-form 对照，优先参照 Bempp 的 `use_strong_form` 方法学：同一 RHS、同一 GMRES 参数、同一停止准则，只切换 weak/strong 离散。
- 若做压缩误差预算试验，优先参照 FreeFEM Maxwell EFIE sphere 的 `eta/eps/minclustersize` 风格参数治理，而不是只改单个经验 `near_range`。

### 4.14 Bempp 风格 block/operator 外壳

理论职责:
- 以四块 `EJ/EM/HJ/HM` 为核心组织 PMCHW transmission 系统，显式承载 block algebra、weak/strong form、mass-matrix-aware solve，并向下调用 Dense/MLFMA/H-matrix backend。

禁止猜测:
- 不允许把 block/operator 外壳再次退化成“先拼一个 2N 大矩阵或大 matvec，再在上面解释物理含义”。
- 不允许让 backend 决定 transmission 的物理分块、符号关系或 strong-form 语义。
- 不允许在 block/operator 外壳还未稳定前并入 H-matrix backend。

开发步骤:
1. 先定义 block/operator 抽象：四块语义、domain/range/dual 语义、`weak_form/strong_form` 接口。
2. 先接入 Dense backend，验证 block 组合与 Direct PMCHW 完全一致。
3. 再把现有 PMCHW MLFMA 收编为 backend，而不是继续扩写单体 `PMCHWMLFMAOperator`。
4. 等 Gate O 与 Gate S 稳定后，再预留 H-matrix backend 接口。

验证方式:
- Gate O1: block/operator 外壳在 dense backend 下必须与现有 Direct PMCHW 在四块级别逐一一致。
- Gate O2: `weak_form` 与 `strong_form` 必须在同一外壳下可切换，且切换不改变 RHS/未知量的物理语义。
- Gate O3: 现有 `PMCHWMLFMAOperator` 若保留，必须仅作为兼容 facade，且其输出与 block/operator + MLFMA backend 路径一致。

## 5. 开发顺序约束

必须按下面顺序推进，跳步视为违规：

1. 几何层: build_octree -> setLevelInfo! -> searchNearCubes -> setKidLevelFarNeighbors!
2. 预计算层: truncationLCal -> levelIntegralInfoCal -> compute_translation_factors!
3. 近场层: assemble_near_field_pmchw
4. 单 pass 层: aggregate_leaf_pmchw! J/k0 -> J/k1 -> M/k0 -> M/k1
5. 末端接收层: disaggregate_leaf_pmchw_j! -> disaggregate_leaf_pmchw_m!
6. 总算子层: constructor -> mul!
7. operator/block 外壳层: block shell -> weak_form/strong_form -> backend binding
8. strong-form 对照层: dense weak/strong -> MLFMA weak/strong
9. 端到端层: Gate C + Gate S + Gate O + impedance gate + regression gate

任何阶段失败，都必须退回当前函数层修正，不允许跳到更上层继续猜。

## 6. 强制门禁矩阵

### Gate A: 结构不变量

- A1: `HJ + EM` near-field block infinity norm ratio < 1e-8。
- A2: `nnz(Z_near) < (2N)^2`，必须存在真实 far field。
- A3: octree permutation、neighbors、farneighbors 的互补关系成立。

### Gate B: pass 级对齐

使用固定随机种子和单位向量两类输入。

- B1: EJ k0 pass vs direct-far EJ-k0。
- B2: EJ k1 pass vs direct-far EJ-k1。
- B3: EM k0 / k1 pass 独立对齐。
- B4: HM k0 / k1 pass 独立对齐。
- B5: 同一 block 的 pass sum 与 direct-far block sum 对齐。

记录指标:

- relative error
- norm ratio
- correlation
- worst column id
- worst row id

### Gate C: 端到端

- C1: `test/test_pmchw_mlfma_operator.jl` 中 `15.11 MLFMA mul! vs Direct` 必须 < 0.10。
- C2: 输入阻抗实部误差 < 5%，虚部误差 < 20Ω。
- C3: 结果满足物理约束，辐射阻抗实部为正。

### Gate D: 回归与可追溯性

- D0: 进入近场逐元素对比前，必须重新验证 Dense PMCHW -> Mie 基线。
- D1: 在 `Z_near` 的非零位置上逐元素比较 dense PMCHW，最大元素相对误差 < 0.1%。
- D2: EFIE/SCFIE 共享 MLFMA smoke tests 通过。
- D3: 用于决策的诊断必须沉淀为正式测试或删除，不允许继续堆积 scratch scripts。

### Gate O: operator/block 架构分界

- O1: PMCHW transmission 顶层必须以四块 `EJ/EM/HJ/HM` 组织，而不是以单体 backend operator 直接暴露主线语义。
- O2: Dense backend 必须先在该外壳下跑通并通过四块级别对照，之后才允许接入 MLFMA backend。
- O3: `strong_form` 必须定义在 block/operator 外壳层，而不是定义在某个特定 backend 的私有实现里。
- O4: FreeFEM 风格 H-matrix backend 在 Gate O 与 Gate S 未稳定前不得进入主线实现。

### Gate S: strong-form 分界

> 2026-03-07 状态补记：`test/test_pmchw_gate_s_dense.jl` 已锁定 dense 半边，`test/test_pmchw_gate_s_mlfma_medium.jl` 已在 `N=540` 夹具上正式跑通 `dense weak / dense strong / MLFMA weak / MLFMA strong` 四路专门回归并写入进度文档；当前已形成 S1-S5 的可执行主对照链。S6-S8 相关的外部真值/外部 fast baseline/压缩预算对照仍保留为后续独立验收项。

- S1: 对同一 `N=540` 夹具，必须同时记录 `dense weak`、`dense strong`、`MLFMA weak`、`MLFMA strong` 四条求解结果。
- S2: 四条结果必须使用同一 RHS、同一 restart/maxiter/停止准则，禁止只改单一分支参数。
- S3: 若 `dense strong` 收回而 `MLFMA strong` 仍不收回，则后续主线问题正式归类为“MLFMA 压缩误差/算子保真度”。
- S4: 若 `dense weak` 与 `dense strong` 都不收回，则后续主线问题不得先归因到 MLFMA，必须先回到 formulation/离散系统定义。
- S5: strong-form 对照结论必须同步写入进度文档，不允许只保留在一次性实验输出中。
- S6: 至少一项单球 dielectric 结果必须与外部真值基线对照；优先使用 scuff-em Mie scattering 或等价解析球结果。
- S7: 至少一项 fast-vs-dense block matvec 结果必须与外部 compressed-matvec 基线量级对照；优先使用 Bempp Maxwell boundary FMM sphere real/complex wavenumber。
- S8: 若引入压缩参数预算试验，必须记录与 FreeFEM/Htool 风格参数 `eta/eps/minclustersize` 的类比关系，禁止只给出“调大了 near_range 更准”这类无预算结论。

## 7. 禁止事项

- 禁止新增一次性诊断脚本留在仓库根目录或 scripts 目录中。
- 禁止把“调参数”写成“修算法”。
- 禁止通过扩大近场、改变 leaf_size_eff、抬高 L_min 来掩盖未定位的错误。
- 禁止在未通过 Gate D0 前把 dense PMCHW 直接当作无条件真相源。
- 禁止在未通过 Gate D1 前讨论 Gate C 结果归因。
- 禁止在未通过单 pass 对齐前合并四遍输出谈总算子正确性。
- 禁止在未完成 Gate S 前，把长 Krylov 偏差草率归类为“纯 solver 问题”或“纯 MLFMA 压缩问题”。
- 禁止在 block/operator 外壳未稳定前，把 FreeFEM 风格 H-matrix backend 直接并入主线。
- 禁止继续把 `PMCHWMLFMAOperator` 单体扩写为后续所有 transmission 能力的唯一承载点。

## 8. 文档检视迭代记录

本次检视对象是本文档本身，不是算法实现结果。目标是确认“文档已经可执行、可审计、可验收”。

### Round 1

发现的问题:

- 原文没有把函数层级拆开，无法知道错误该回退到哪一层修。
- 原文没有把 near/far 划分入口锁死，导致 near_range 修改仍可能跨层失控。

修正:

- 新增第 3 节分层对象。
- 新增第 4 节函数级合同。
- 新增第 6 节 Gate A/B/C/D/S 对应的层级映射。

结论:

- 本轮有新发现，不能通过。

### Round 2

发现的问题:

- 原文虽然有 Gate，但没有把“禁止猜测”的边界写死，仍可能以经验常数或参数漂移绕过理论定位。
- 原文没有明确诊断产物的处置规则，容易再次堆积废物脚本。

修正:

- 在每个关键函数合同中加入“禁止猜测”条款。
- 新增第 7 节禁止事项。
- 在 Gate D 中新增“诊断必须沉淀为正式测试或删除”的要求。

结论:

- 本轮有新发现，不能通过。

### Round 3

检查结果:

- 分层完整，已覆盖几何、预计算、近场、远场、总算子、验收六层。
- 每个关键函数均已具备理论职责、禁止猜测、开发步骤、验证方式四项信息。
- 门禁、开发顺序、禁止事项、退出条件已经闭合，没有发现新的文档级缺口。

结论:

- 本轮无新发现。
- 文档级治理合同通过，可作为后续实现唯一执行文档。

## 9. 完成定义

一个实现任务完成，必须同时满足：

- 有对应函数级合同。
- 有 RED/GREEN/REFACTOR 的测试与提交证据。
- Gate A/B/C/D/S/O 全部通过并记录结果。
- Gate S 对照已完成并给出明确归因。
- 本文档未被违反。
- `.github/REFACTORING_ROADMAP.md` 与 `.github/REFACTORING_PROGRESS.md` 已同步更新。

## 10. 关联位置

- Roadmap: `REFACTORING_ROADMAP.md` -> Phase 15 -> `15.G1`
- Progress: `REFACTORING_PROGRESS.md` -> 2026-03-06 文档治理刷新记录
