# Findings: MoM_Lebedev

（随调研持续补充）

## 代码结构（2026-08-08 快照）
- `LebedevSortedPoints.jl`: 读取 `deps/sphere_lebedev/nodesSorted/` 下 `p.n.txt` Lebedev 节点/权重；
  `modiTgetFileName` 对不存在的阶数向上取最近文件；`get_t_nodes(t)` 在 p=2t+1 超出上限时退回
  GL(θ)×uniform(ϕ) 张量积网格。
- `dataset_generator.jl`: 用随机 EFIE 双电流源聚合生成 (pk→pt) 数据集，k=2π 固定，无 RNG 种子；
  `pt` 参数被重新赋值（`pt = 2τt + 1`）。
- `pinv2interpW.jl`: IDW(kNN) 初始化 → 逐行 `pinv` 最小二乘拟合 → 80/20 测试误差 → 精度更好才存 h5。
- `LVI.jl`: `levelIntegralInfoCal(::Val{:LbTrained1Step})` 选择 Lebedev 阶数；
  `interpolationCSCMatCal` 按 n2pDict 反查阶数并加载/生成 W；interpolate/anterpolate 为稀疏 matvec。

## 关键疑点（待验证）
1. `LVI.jl:39` 向 `truncation_kernel` 传 `levelCubeEdgel*2π/λ`（ka），基实现 `Interpolation.jl:91` 传 `cubel/λ`（a/λ），
   疑似把 L 放大 ~2π 倍 → Lebedev 阶数/点数被系统性高估。
2. `LVI.jl:49` 条件 `2truncL+1 < maximum(keys(p2nDict))`，当需要 p=131 时 131<131 为 false → 最大 Lebedev 集永不使用。
3. `deps/InterpolationWeights/` 为空 → 首次运行 MLFMA 会现场训练权重（慢、占内存、不可复现）。
4. 训练数据与真实 MLFMA 数据分布不一致（固定 k、特定 EFIE 聚合、随机几何），且 2×2 块引入了 θ↔ϕ 交叉耦合。

## 数值验证结果（2026-08-08，julia 1.12.3，脚本 scripts/verify_lebedev_interp.jl / verify_lebedev_fullsize.jl）

### A. truncation_kernel 参数不一致（已确认）
| a/λ | L_correct(a/λ) | L_LVI(ka) | 放大 |
|-----|---------------|-----------|------|
| 0.125 | 9.98 | 24.46 | 2.45x |
| 0.25  | 13.59 | 37.14 | 2.73x |
| 0.5   | 19.13 | 59.45 | 3.11x |

`levelIntegralInfoCal(0.5, Val(:LbTrained1Step))` 返回 truncL=60（应 19）。
实际取到 p=125 -> 5294 点；按正确 L=19 -> p=39 -> 最近文件 p=41 -> 590 点，点数虚高约 9 倍。

### B. p=131 死代码（已确认）
`2truncL+1 < 131`：truncL=65 时 131<131=false → 回退 Lagrange。最大 Lebedev 集（131.5810.txt）永不使用。

### C. 训练式插值权重的实测表现（pk=13 -> pt=27）
| 方法 | EFIE 留出集平均相对误差 | 限带函数(deg<=6) 最大相对误差 | 向量场最大相对误差 | 行和 max|W*1-1| |
|------|----------------------|------------------------------|--------------------|------------------|
| IDW 初始化 | 1.129 | - | - | - |
| 现行管线（仅实部约束） | 1.475（完整数据集 0.79-0.89） | 3.605e+02 | 6.011e+02 | 2.785e+02 |
| 修正版（实+虚约束） | 0.550 | 4.137e+01 | - | - |
| 球谐精确插值 W=Y_fine*(Y_coarse)^+ | - | 3.983e-15 | 4.823e-15 | 1.776e-15 |

- SH 基正交性验证: 3.55e-14（Lebedev 权重下近似单位阵）→ 球谐实现正确。
- 现行 W: nnz=8692（每行 16.3）；球谐 W 稠密 266x74=19684 非零。
- 训练时耗（完整 50ρx500pos，最小一对 13->27）：数据集 13.2 s + 逐行 pinv 11.0 s。
- 可复现性：两套随机种子的 W 输出差异 12.6（seed A 误差 8.97，seed B 11.4）。
- 内存：pk=65->131 时 pArray 预分配 4.65 GB（不可行）。

### D. 只使用实部约束（代码级确认）
`pinv2W!` 中 `w` 只有 2np 行，`xx2D/yy2D = vcat(real, imag)` 有 4np/4nt 行 → 拟合只用实部，
虚部行永远不被索引。修正版（实+虚）误差从 1.475 降到 0.550，但局部 9 点支撑仍是根本瓶颈。

## 更好的插值权重方式（结论）
1. **球谐最小二乘全局插值（推荐）**：W = Y_fine (Y_coarse)^+，degree <= (p_coarse-1)/2；
   θ/ϕ 分量独立标量插值（分块对角），无交叉块；确定性、无训练数据；限带函数机器精度（实测 3.98e-15）。
2. **局部支持 + 多项式精确性约束（保持稀疏的一步插值）**：支撑取球冠内全部点（或 >= 2(L+1)^2），
   约束 W*Y_coarse = Y_fine（对 degree<=L 精确）→ 行和自动为 1。
3. **工程兜底**：主链路已用 GL 张量积 + 两段 Lagrange（现有 LagrangeInterpInfo），对限带函数可控精度，
   无训练成本；Lebedev 仅在需要避免极区奇异性/省点数时启用 1/2。

## 文献
- X.-Y. He, D.-H. Kong, "Efficient Vector Interpolation Method Based on Lebedev Quadrature for MLFMA",
  IEEE Antennas Wireless Propag. Lett., 21(5):1056-1059, 2022, DOI 10.1109/LAWP.2022.3158397
  （本代码的理论来源，第一作者 xyhe 即本工程作者）。
- 球谐展开 MLFMA（SE-MLFMA）：Eibert 等，"Improving the Spherical Harmonics Expansion-Based MLFMA"。
- Eibert, "Performing interpolation and anterpolation entirely by FFT in 3-D MLFMA"。
- Ergül & Gürel, "Enhancing the accuracy of the interpolations and anterpolations in MLFMA", 2006。

---

# 新目标进展（2026-08-08）：复现 -> 结构分析 -> 优化到极限

## R1: 原始实现确认（github.com/deltaeecs/MoM_Lebedev.jl 归档仓库）
- 原始方法与合并版一致：`interpWeightsInitial`(KDTree kNN 反距离) +
  `pinv2W!` 逐行 pinv + `vcat(real,imag)` 布局（只实部约束）+ 2x2 交叉块。
- 原始 LVI 用 `truncationLCal`（正确 a/λ）；合并版引入 `*2π` 回归 bug。
- 原始包每次构造 `LbTrainedInterp1tepInfo` 都重新 `runpinvCal`（!ispath 被注释）。
- 旧权重（jld2）从未提交到 git；需要从零复现。

## R2: 复现结果（与论文数值对齐需论文 PDF，已向用户索取）
| 对 (pk->pt) | EFIE 留出平均相对误差 | 限带函数最大相对误差 | nnz/行 |
|---|---|---|---|
| 13->27 | 1.475（完整数据集 0.79-0.89） | 3.6e2 | 16.3 |
| 27->53 | 0.558 | 1.1e2 | 15.8 |
| 65->131 | 不可行（pArray 4.65 GB + 逐行 pinv） | - | - |

重要发现：EFIE 训练数据的实际带宽（源间距最大 0.12λ -> 阶数 ~11）超出粗网格
p=13 的可表示范围（degree<=6-7，74 点），球谐精确插值在该数据上也只有 ~0.98 误差
——原始方法在拟合粗网格根本表示不了的内容，这是"不理想"的结构性根因之一。

## R3: Lebedev 点结构
- 每阶点集 = 24 元素八面体群（真旋转）的轨道并集：轨道大小 6/8/12/24。
  p=13: 5 轨道 [6,8,12,24,24]; p=27: 13 轨道; p=41: 26; p=65: 62; p=131: 244。
- 轨道数约 = 点数/24 -> 权重构造与存储可压缩约 24x。
- 纬度环的 φ 非等间距（含 0 间距重复）-> 无直接 ring-FFT；旋转群等变是主要可用结构。
- 网格对 24 旋转群精确不变；权重和 = 4π；球谐正交性在精确度范围内成立。

## R4: 优化实现（src/FastAlgorithms/Lebedev/SHInterp.jl + LVI.jl 修复）
| 方法 | 13->27 限带误差 | 27->53 限带误差 | nnz/行 | 备注 |
|---|---|---|---|---|
| 原始管线（复现） | 3.6e2 | 1.1e2 | 16 | 训练式、不可复现 |
| 局部约束 L_loc=3/4/6 | 1.4 / 5.0 / 2.9e-13 | 0.59 / 0.50 / 0.87 | 18-61 | 行和=1 精确 |
| 局部约束全度 | - | 4.5e-5 | 234 | 近全局支撑 |
| 球谐精确（稠密） | 5.1e-15 | 3.6e-15 | 74-266 | 机器精度、确定性 |
| 轨道压缩 | 与朴素一致 (1e-12/7e-14) | 同左 | 同左 | 构造/存储 ~1/24 |

- LVI 修复：truncation_kernel(levelCubeEdgel/λ)（原 *2π 放大 2.45-3.1x）；
  `2L+1 <= max(p)`（原 off-by-one 使 p=131 死代码）。修复后 0.5λ 盒 590 点（原 5294）。
- `LbTrainedInterp1tepInfo` 新增 method=:sh_exact / :sh_local / :sh_local_orbit；
  默认 :trained 保持原行为。test_lebedev.jl 9/9 通过。

## R5: 结论与极限
- 精度极限：球谐精确 W（限带函数机器精度），但稠密；高阶层应用成本高。
- 稀疏与精确不可兼得：一步插值核本质全局（要精确到 degree L 需 ~(L+1)^2 支撑）；
  局部约束方案是"精度-稀疏"的可调折中，L_loc<=6 时 55-90 点/行。
- 性能极限组合：正确阶数选择（省 ~9x 点数） + 球谐权重（省训练、机器精度）
  + 轨道压缩（省 ~24x 构造/存储） + :sh_local 稀疏版按需折中。

---

# 权重矩阵修复记录（2026-08-08）

## 根因链（用户："现有的权重矩阵全都是坏的"）
- 直接证据：修复前所有方法（原始训练式、修正复约束、球谐精确、全支撑 LS）在
  EFIE 留出集上的论文指标 εi 都是 O(1)，且点数越多误差越大（16->24 点达 1e2~1e3）。
- **主因**：`generate_dataset_on_pkpt` 对粗/细层分别调 `generate_dataset_on_poles`，
  后者内部每次重新生成随机源几何 -> 两层样本来自不同辐射函数（噪声拟合）。
  判别实验：内联"共享几何"生成器 k=8 得 εi=2.15e-4，真实生成器（不共享）得 1.04。
- 次因：`vcat(real,imag)` 布局只让实部参与逐行 pinv；固定臂长 0.12λ 与低阶盒子失配。
- 频谱实验：原版生成器数据 99.9% 能量阶数远超粗网格容量（p=13 容量 degree<=6-7）。

## 修复内容（源码）
- dataset_generator.jl：新增 random_source_geometry / evaluate_poles!；
  generate_dataset_on_pkpt 每样本生成一次几何并供粗/细层共用；几何按 rel_l 缩放；
  兼容旧调用（geom 关键字、默认 arm_max/off_max）。
- pinv2interpW.jl：runpinvCal 中 vcat(real,imag) -> hcat(real,imag)（全复数约束）。
- SHInterp.jl：interp_weights_exact 默认 degree=Lb（粗网格精确表示极限，非 Lb+1）；
  新增 interp_weights_auto（按规模选精确/局部稀疏）。
- LVI.jl：interpolationCSCMatCal 默认 method=:sh_auto；truncation/off-by-one 已修。

## 修复后验证
- 论文 Fig.2 复现：(13->27) k=8 εi=2.46e-4, k=18=2.7e-6；(27->53) k=8=1.03e-3；
  完整 25000 样本 k=9：εi 均值 1.51e-4 / 最大 4.92e-4 —— 与论文 ~6e-4 同量级且趋势一致。
- 球谐精确：确定性、行和 2e-15、限带矢量场 4.8e-15。
- 局部约束/轨道压缩：行和精确、orbit==naive（1e-12）。
- test_lebedev.jl 9/9 通过。

## 已知残差
- EFIE 矢量场的 L+1 切向投影耦合超出低阶粗网格精确表示极限；标量球谐在该数据类
  εi~0.6；局部数据拟合更优。理论完备解：矢量球谐（VSH）插值（后续工作）。

# MLFMA 集成与真实基准（2026-08-08）

## 集成点
- interp_type 映射 + OctreeBuilder(interp_method) + Aggregation/Disaggregation 一步分支
  + MLFMAOperator(interp_method, near_range)。

## 基准结果（PEC 球 EFIE）
- r=1λ（N=792）：GL≡Lebedev（全算子差 5.6e-9），εq=1.344e-2（该规模插值不进远场）。
- r=2λ（N=1440）：修复后训练式 k=8 εq=1.160e-1 == GL Lagrange，插值内存 1.5 MB vs GL 0.64 MB，
  极点数 -29%，MVM -35%（0.225s vs 0.349s）；原版坏权重 εq=6.975（60x 差）；零矩阵 7.01。
- 稠密球谐精确 εq=0.372：标量 SH 无法吸收矢量 L+1 耦合；数据拟合（交叉块）更优。
- EMMoMSuite MLFMA 远场校准按 near_range=4；near_range=1 有既有远邻居空隙。
