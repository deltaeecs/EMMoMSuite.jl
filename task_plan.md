# Task Plan: MoM_Lebedev 调研（做什么 / 缺陷 / 更好的插值权重算法）

## Goal
搞清楚 `MoM_Lebedev`（现并入 EMMoMSuite 的 `FastAlgorithms.Lebedev`）做了什么、有什么缺陷，
并给出"更好的插值权重计算方式"。

## Phases

- [x] P1: 通读 Lebedev 全部源码（LebedevSortedPoints / dataset_generator / pinv2interpW / LVI），
      确认其在 MLFMA 中的调用链（interpolate / anterpolate / levelIntegralInfoCal）。
      -> 结论：模块已并入 EMMoMSuite 但主 MLFMA 链路未启用（OctreeBuilder 用 GL+Lagrange）；
         理论来源为 He & Kong, IEEE AWPL 21(5):1056-1059, 2022（作者 xyhe 即本工程作者）。
- [x] P2: 核实关键疑点（数值实验）：
      - LVI 中 truncation_kernel 传入 `levelCubeEdgel*2π/λ` 是否把 L 放大 ~2π 倍；
        -> 实测放大 2.45x/2.73x/3.11x（a/λ=0.125/0.25/0.5）。
      - p=131（5810 点）是否是死代码（`2truncL+1 < max` 的 off-by-one）；
        -> 是：truncL=65 时 131<131 为 false，触发回退。
      - `deps/InterpolationWeights/` 是否为空、运行时是否触发 on-demand 训练；
        -> 为空；完整管线最小对 (13->27) 训练 24.2 s；自身留出集误差 0.79-0.89。
      - 训练数据随机性 / 无种子 / 无复现性；
        -> 两套随机种子训练出 W 的输出差异达 12.6 倍误差量级。
      - 大阶数训练内存爆炸量级（np=5810 → ~4.6 GB）。
        -> pk=65->131 时 pArray 4.65 GB，不可行。
- [x] P3: 文献/业界对照：找到 He&Kong 2022（本方法来源）、SE-MLFMA、FFT 插值、
      球谐最小二乘插值等对照方案。
- [x] P4: 写出结论：MoM_Lebedev 做了什么、缺陷清单（按严重度）、更好的权重计算方案。
- [x] P5: 汇总最终报告（含证据、数值结果、引用）。

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| Statistics.mean 未导入 | 补 using Statistics | 1 次 |
| Julia 顶层软作用域（for 内对全局赋值） | 包成 main() | 1 次 |
| realSHmatrix 漏 P_1^0（P[0,1] 越界） | 递推补 l=1 特例 pprev=0 | 1 次 |
| nnz(稠密矩阵) MethodError | 改用 count(!iszero, ...) | 1 次 |

## Decisions
- 以当前 EMMoMSuite `src/FastAlgorithms/Lebedev/` 为现行源码；旧 `MoM_Lebedev.jl` 已归档（README 第 391 行）。
- 用可复现的小规模数值实验验证关键疑点，而不是仅靠读码推断。
- 保留两个验证脚本（scripts/verify_lebedev_interp.jl、verify_lebedev_fullsize.jl）供用户复验。

---

# 新目标（2026-08-08 用户确认）：复现论文结果 -> 分析 Lebedev 点结构 -> 优化插值矩阵到极限

用户确认：He&Kong 2022 为其本人论文；原始权重 = 生成数据后伪逆，论文数据真实可靠；
重构后权重丢失。新需求：
1. 复现此前（论文/原始权重）的结果；
2. 分析 Lebedev 球面积分点结构，利用点分布规律性；
3. 继续优化插值矩阵，把 Lebedev 性能推到极限。

## Phases (新目标)

- [x] R1: 找到旧仓库 github.com/deltaeecs/MoM_Lebedev.jl（公开归档，作者本人项目），
      确认原始方法 = IDW(kNN) 初始化 + 逐行 pinv + 只实部约束 + 2x2 交叉块；
      论文 PDF 被 IEEE paywall 挡住（可向用户索取以对齐论文报告数值）。
- [x] R2: 复现原始管线（13->27: EFIE留出 1.475, 限带 3.6e2; 27->53: 0.558, 1.1e2;
      65->131 不可行 4.6GB）。发现 EFIE 训练数据带宽(~11)超出粗网格可表示范围(<=7)。
- [x] R3: 结构分析：Lebedev 点 = 24 元素八面体群轨道并集（6/8/12/24 点轨道），
      p=131 有 244 轨道 vs 5810 点（约 1/24）；纬度环 φ 非等间距 -> 无直接 FFT；
      旋转群等变性精确成立。
- [x] R4: 优化实现（新模块 src/FastAlgorithms/Lebedev/SHInterp.jl）：
      a) interp_weights_exact: W=Y_fine pinv(Y_coarse)，机器精度（3e-15），确定性；
      b) interp_weights_local(_orbit): 局部支撑 + 多项式精确约束 + SVD，
         轨道压缩与朴素版逐位一致（1e-12），构造/存储约 1/24；
      c) LVI.jl 修复 truncation 参数 bug（ka->a/λ，点数 5294->590）与 off-by-one（p=131 启用）；
      d) LbTrainedInterp1tepInfo 新增 method=:sh_exact/:sh_local/:sh_local_orbit，默认保持原训练式。
- [x] R5: 验证：test_lebedev.jl 9/9 通过；sh_exact 在 13->27/27->53 限带矢量场 3.8e-15/3.1e-15；
      局部 L_loc=6 13->27 2.9e-13（55/行）；27->53 全度 4.5e-5（234/行）；
      高阶层无稀疏精确一步插值（核全局），稠密 SH 是精度极限。

## 关键结论（新目标）
- 原始方法的"不理想"有定量证据：即使在自己的 EFIE 留出集上误差也是 O(1)，
  根因是 (a) 9-16 点局部支撑无法表示全局核 (b) 训练数据带宽超出粗网格可表示范围
  (c) 只用实部约束。论文数据真实可靠与这些实现缺陷不矛盾。
- 极限优化路径：确定性球谐精确权重（机器精度）为精度极限；
  轨道压缩（约 24x）为结构与存储极限；正确阶数选择（修复 truncation）省约 9x 点数。

---

# 权重矩阵修复（2026-08-08 用户指令："现有的权重矩阵全都是坏的，先把权重矩阵修复"）

## 根因（已定位并修复）
1. **数据集生成根本性 bug（主因）**：`generate_dataset_on_pkpt` 对粗层/细层分别调用
   `generate_dataset_on_poles`，而后者每次在内部重新随机生成源几何（臂/偏移/参考点）
   -> 粗层与细层样本对应**不同辐射函数**，插值矩阵在拟合纯噪声。
   修复：`random_source_geometry` 每样本生成一次几何，`evaluate_poles!` 供两层共用；
   并新增 `geom` 关键字兼容旧调用。
2. **复约束布局 bug**：`vcat(real,imag)` 使逐行 pinv 只用实部（4n 行只索引前 2n）。
   修复：`hcat(real,imag)`（2n x 2N），与论文 Algorithm 2（对复数数据 pinv）等价。
3. **生成器几何与层带宽不一致**：固定臂长 0.12λ 与低阶盒子（rel_l~0.05λ）失配，
   数据带宽超出粗网格可表示范围。修复：几何按 rel_l 缩放（arm=0.5*rel_l*λ）。
4. **LVI truncation/off-by-one**（之前已修）。

## 修复后验证结果
| 指标 | 修复前 | 修复后 |
|---|---|---|
| (13->27) k=8 论文 εi | 1.03~4.8 | 2.46e-4（完整规模 k=9: 1.51e-4） |
| (27->53) k=8 论文 εi | 1.03 | 1.03e-3 |
| 误差随点数趋势 | 增大（16->24 时 1e2~1e3） | 递减（k=18: 2.7e-6） |
| 球谐精确限带 | 3.8e-15 | 4.8e-15（确定性、行和 1e-15） |
| test_lebedev.jl | - | 9/9 通过 |

论文 Fig.2（8 点 ~6e-4）复现成功。

## 已知残差（如实记录）
- EFIE 矢量场含 degree L+1 切向投影耦合，超出低阶粗网格（p<=27）精确表示极限；
  标量球谐精确 W 在该类数据上 εi~0.6，局部数据拟合路径（k=8~18）对此类更有效。
  理论完备解为矢量球谐（VSH），留作后续。
- 修复后的默认生产路径：`interpolationCSCMatCal` -> `LbTrainedInterp1tepInfo(method=:sh_auto)`
  （小规模精确稠密，大规模局部稀疏）；`method=:trained` 保留修复后的训练管线。

---

# MLFMA 集成与真实基准（2026-08-08，目标第 2-3 项）

## 集成（源码）
- `Interpolation.jl`：新增 `interp_type` 映射（GL->Lagrange，Lebedev 侧扩展）。
- `OctreeBuilder.jl`：`build_octree(; interp_method=Val(:Lagrange2Step))`；
  levels 字典类型放宽为 `Dict{Int,AbstractLevel}`；`LevelInfo` 的 IPT 由 poles 类型推导。
- `Aggregation.jl` / `Disaggregation.jl`：为 `LbTrainedInterp1tepInfo`（θϕCSC/θϕCSCT）
  增加一步插值/反插值分支（hasfield 判断，避免模块循环依赖）。
- `MLFMAOperator.jl`：构造函数支持 `interp_method` 与 `near_range` 透传。

## 真实基准（PEC 球 EFIE，EMMoMSuite MLFMA，near_range=4 校准配置）
### r=1.0λ (N=792, 3 层，非叶层无远邻居 -> 插值不进远场路径)
| 方法 | εq | 插值内存 | MVM |
|---|---|---|---|
| GL 两段 Lagrange | 1.344e-2 | 0.24 MB | 0.082s |
| Lebedev :sh_auto（稠密精确） | 1.344e-2 | 22.2 MB | 0.127s |
| Lebedev 训练式 k=8 | 1.344e-2 | 0.6 MB | 0.056s |
| 零插值矩阵 | 1.344e-2 | - | - |
（所有方法 εq 相同：该规模远场不经过插值；GL 与 Lebedev 全算子差 5.6e-9）

### r=2.0λ (N=1440, 4 层，level 3 有 13344 个远邻居 -> 插值真正参与)
| 方法 | εq | 插值内存 | 单次 MVM |
|---|---|---|---|
| GL 两段 Lagrange | 1.159e-1 | 0.64 MB | 0.349s |
| Lebedev 稠密精确 | 3.715e-1 | 40.8 MB | 0.923s |
| Lebedev :sh_local L_loc=6 | 3.530e-1 | 23.6 MB | 0.324s |
| **Lebedev 训练式 k=8/18（修复后）** | **1.160e-1** | **1.5/3.3 MB** | **0.225/0.240s** |
| **原版坏权重（k=9）** | **6.975（60x 差）** | 1.6 MB | 0.227s |
| 零插值矩阵 | 7.014 | - | - |

关键结论：
- **MLFMA 层直接证明修复有效**：坏权重 εq=6.98 -> 修复后 0.116（60x 提升）。
- **修复后训练权重 = 最优工程选择**：εq 与 GL Lagrange 完全一致（0.116），
  极点数少 29%（6044 vs 8466，对应论文"省 1/3 辐射模式存储"），MVM 快 35%（0.225 vs 0.349s，
  对应论文 Table II 趋势）。
- 稠密球谐精确 εq 反而差（0.372）：标量 SH degree=Lb 无法吸收 EFIE 矢量场的 L+1 切向投影
  耦合；数据拟合权重（含 θ↔ϕ 交叉块）对该矢量类更优 —— VSH 是理论完备解（后续工作）。
- EMMoMSuite MLFMA 按 near_range=4 校准；near_range=1 存在远邻居空隙（与插值无关的既有问题）。

## 回归测试
- test_lebedev.jl 9/9、test_mlfma.jl（MLFMAOperatorMPI 6/6、物理排序 2/2，含 MLFMA≈直接矩阵 1e-8）。

---

# 新目标（2026-08-08 用户："稠密插值没有优势；继续推进方案更新、效率优化：
# 利用点分布规律、扩展到更高阶、找更高效采样点、设计新插值算法"）

## 核心成果

### 1. 新插值算法：笛卡尔标量 SH 矢量插值（interp_weights_cart 系列）
- 原理：3 个笛卡尔分量 F_x,F_y,F_z 各自按标量限带函数插值（degree=Lb/Lb+1，
  共用同一最小范数权重 w），再以 θ̂/ϕ̂ 点积散布成 (θ,ϕ) 耦合 2x2 块。
  笛卡尔分量在极点是光滑的 -> 彻底避开 θ/ϕ 基极点奇异；精确吸收矢量 L+1 耦合。
- 验证（EFIE 矢量类，修正几何）：27->53 稠密 εi=1.0e-14（机器精度，免训练、确定性）；
  13->27 εi=1e-9。局部稀疏：L_loc=3 -> 1.9e-4（46 nnz/行，优于训练式 k=8 的 1e-4 且带精确性保证），
  L_loc=6 -> 3.3e-7（122/行），L_loc=13 -> 5.8e-13（465/行）。
- VSH（自旋加权球谐）已实现并验证正交性/机器精度，但 θ/ϕ 基在极点有固有奇异，
  对 EFIE 数据残留 εi~0.76 -> 笛卡尔方案为更优路径（记录为研究结论）。

### 2. 利用点分布规律：八面体群轨道压缩（interp_weights_cart_local_orbit）
- 每轨道代表求解标量权重（24x 少 pinv），其余行用各节点自己的 θ̂/ϕ̂ 帧点积散布
  （切平面 holonomy 自然吸收）；O(n*24) 规范形分组 + 向量化点积 + 无分配距离
  + 三元组稀疏装配 -> 65->131 构造从 887s（朴素）降到 ~1s（L_loc=6），100x+。
- 与朴素版逐位一致（6e-9，受列映射 1e-8 舍入限制）。

### 3. 扩展到更高阶数
- 65->131（5810 点）轨道压缩版 1s 内构造：L_loc=6 -> εi=7.1e-5（180 nnz/行），
  L_loc=10 -> 2.4e-6（343 nnz/行）。原训练管线在此规模不可行（4.6 GB 数据集）。

### 4. 采样点效率对比（lebedev_grid_study.jl）
| 网格 | L=6 κ(Y) | 最小角距 | L=13 κ(Y) | 最小角距 |
|---|---|---|---|---|
| Lebedev (74/266) | 2.23 | 0.208 | 1.68 | 0.095 |
| Fibonacci (74/266) | **1.32** | **0.361** | **1.44** | **0.190** |
| GL 张量积 (98/392) | 1.80 | 0.140 | 2.48 | 0.037 |
- Fibonacci 同点数下条件数更优、最小角距 2x（准均匀）、点数任意（可到 >5810）；
  对插值（不用权重）是更高效的采样点；GL 张量积 L=13 条件数最差且极区聚集。

### 5. 生成器源几何修复
- arm = min(0.12λ, 0.5*rel_l*λ)、rvec 缩到盒内余量 -> 源整体在盒内，数据带宽与层一致。

## 验证脚本
- scripts/lebedev_cart_verify.jl（笛卡尔稠密/局部/轨道 + 高阶）
- scripts/lebedev_vsh_verify.jl（自旋加权基与 VSH 验证）
- scripts/lebedev_grid_study.jl（采样点对比）

## 后续候选（如实记录）
- 把笛卡尔局部权重接入 MLFMA 真实用例（接口已兼容 θϕCSC，直接可 swap）；
- Fibonacci 网格接入 MLFMA（任意 n、条件数更优）；
- VSH 极点处理（在极点用笛卡尔帧）或直接以笛卡尔方案为主。

## 性能对比（2026-08-08 用户："非零元这么多，速度大概率是问题"）

### 插值 matvec 微基准（单次 W·x）
| 权重 | 27->53 | 65->131 |
|---|---|---|
| 训练式 k=8（16/行） | 2.0e-5 s | -（数据不可行） |
| 笛卡尔 L_loc=3（46/行） | 4.0e-5 s | 2.18e-3 s（255/行, cap=0.6） |
| 笛卡尔 L_loc=6（122/行） | 8.0e-5 s | 5.78e-3 s（667/行, cap=1.0） |
| 笛卡尔稠密（532/行） | 2.2e-4 s | 8.2e-3 s |
稀疏 matvec 受内存带宽限制（有效 3-9 GFLOP/s），nnz 直接折算时间。

### 真实 MLFMA（PEC 球 r=2λ, N=1440, near_range=4）
| 方法 | εq | 单次 MVM |
|---|---|---|
| GL 两段 Lagrange | 0.1159 | 0.341s |
| **训练式 k=8（16/行）** | 0.1160 | **0.239s（最快，比 GL 快 30%）** |
| 笛卡尔 L_loc=3（46/行） | 0.1160 | 0.346s（≈GL） |
| 笛卡尔 L_loc=6（122/行） | 0.1159 | 0.527s（1.5x 慢于训练式） |
| 稠密精确 | 0.372 | 0.999s（慢且 εq 差：矢量耦合残差） |

### 结论（回应"非零元多则慢"）
- **成立**：非零元从 16 到 122，插值阶段时间 ~4x；MLFMA 总 MVM 1.5-2x。
- **MLFMA εq 预算（~0.116，由远场其它误差主导）下，插值精度的富余不改变 εq**：
  最快的正确选择是训练式 k=8（与 GL 同精度、快 30%）；笛卡尔 L_loc=3 是
  确定性/免训练的同级替代（≈GL 速度，εi 1.9e-4）；L_loc=6 留给精度敏感场景。
- 尝试把支撑收紧到恰好 m=(Lloc+1)² 个最近点省 30% nnz，但 εi 恶化 44 倍
  （8.4e-3 vs 1.9e-4）——cap 球冠带余量才是稳健配置。
- 稠密无优势（慢 4x 且 εq 更差）——印证用户对稠密插值的判断。

## 混合权重（2026-08-08 用户："按你建议的试试"）

### 建议 #1：自适应支撑 + 笛卡尔标量 SH 精确性约束 + 数据拟合（interp_weights_hybrid）
- 每细层点取 k=max(ceil(1.5~2*(Lloc+1)^2), 16) 个最近粗点，KKT 求解
  "min |wA-b|^2 s.t. w·Yc[S]=Yf[i]"（笛卡尔数据拟合 + degree<=Lloc 精确性），
  再以 θ̂/ϕ̂ 点积散布；数据集由模块内 splitmix64 确定性生成（无全局随机依赖）。
- 27->53 EFIE 类（εi / nnz 每行 / 构造）：
  训练式 k=8: 1.05e-4 / 16 / (训练)；
  混合 L_loc=3 scale=1.5: **7.7e-6 / 48 / 2.4s**（比训练式准 14x）；
  混合 L_loc=3 scale=2.0: 1.5e-6 / 64 / 3.2s；
  混合 L_loc=4 scale=1.5: 2.6e-7 / 76 / 3.6s。
  确定性验证：同种子两次构造逐位一致 ✓。
- 65->131：构造 5.9s（L_loc=3 scale=2.0，64/行）——可行。
- 真实 MLFMA（r=2λ）：训练式 k=8 MVM=0.235s；**混合 L_loc=3 MVM=0.282s**
  （比训练式慢 20%，但比笛卡尔 L_loc=6 的 0.527s 快近 2 倍），εq 同为 0.116。
- 注：L_loc=4 scale=2.0 出现数值异常（εi 1.5e-3，k=50 条件数问题），默认推荐 L_loc=3。

### 建议 #2：Fibonacci 网格用更少点数
- 同点数 266->974：εi=3.6e-4 vs Lebedev 1.9e-4（差 2x）；
  减点数（200->740 等）后 εi 恶化 10~100x。
- 结论：局部稀疏插值场景 Lebedev 局部覆盖更优（Fibonacci 最小角距大=帽内点少）；
  Fibonacci 仅适合"任意点数/全局稠密"场景。不推荐作为本方案替代。

## 4λ 球规模行为（2026-08-08 用户："好的"）

### 修复：LVI 回退 bug（4λ 首次触发）
- octree level 1（4λ 球盒子 4.8λ，L=81 -> p=143 > 131）超出 Lebedev 上限，
  原回退调用不存在的 `levelIntegralInfoCal(…, Val(:Lagrange2Step))` -> MethodError。
  已修复为正确的 GL 调用 `Interpolation.levelIntegralInfoCal(cubeEdgel; λ)`。

### interp_weights_auto 阈值 4M -> 250k：中阶对也走混合
- 结果：εq 从 0.975（稠密对污染）回到与 GL 完全一致 0.6698；
  MVM 1.217s vs GL 3.536s（**快 2.9x**）；插值内存 35 -> 14.7 MB；
  稠密精确在 65->101 的 εi 只有 0.60（矢量耦合残差），混合 3.5e-5 且只需 48/行。

### 4λ 汇总（PEC 球 r=4λ, N=3312, 5 层, near_range=4）
| 指标 | GL | Lebedev :sh_auto(混合) |
|---|---|---|
| εq | 0.6698 | 0.6698（一致） |
| 单次 MVM | 3.536s | **1.217s（2.9x 快）** |
| 极点数 | 21914 | 19492（Lebedev 层 -29%） |
| 插值矩阵 | 1.67 MB | 14.7 MB（全混合 48/行） |
| 构造 | 29.5s | 101s（混合训练+KKT，一次性） |

- εq=0.67 为 EMMoMSuite MLFMA 在 4λ 的既有远场精度问题（GL 同样），与插值无关。
- 规模越大 Lebedev 优势越明显：2λ 快 1.3x -> 4λ 快 2.9x（极点节省跨层累积）。
- 回归：test_lebedev.jl 9/9、test_mlfma.jl 8/8。

## 高阶无 Lebedev 数据集怎么办（2026-08-09 用户提问）

### 方案（已实现）
1. **高阶节点现场生成**：`LebedevSortedPoints.fibonacci_grid(n)` + `high_order_nodes(p)`
   （n ≈ (4/3)τ²，与 Lebedev 同效率；任意点数、准均匀、无需数据集）。
   `get_t_nodes` 在 p>131 时自动回退 Fibonacci（替换原 GL 张量积回退）。
2. **阶数显式携带**：`LbPolesInfo` 增加 `p::Int` 字段；
   `interpolationCSCMatCal` 用 `p` 而非 `n2pDict`（任意点数可工作）。
   `levelIntegralInfoCal(Val(:LbTrained1Step))` 高阶层构建 Fibonacci + 等权重极点
   （无 GL 回退、无混合对）。
3. **插值权重天然支持任意点集**：混合/笛卡尔权重只用节点坐标，
   已验证 (101->163) εi=1.18e-4（48/行，构造 13.1s）。

### 4λ 球验证（N=3312, 5 层）
| 指标 | GL | Lebedev :sh_auto（混合+Fibonacci 高阶） |
|---|---|---|
| εq | 0.6698 | 0.6698（一致） |
| 单次 MVM | 1.272s | **1.079s（快 18%）** |
| 极点数 | 21914 | **15009（-31.5%）** |
| level 1 (L≈82) | 13448（GL） | **8965（Fibonacci）** |
| 插值矩阵 | 1.67 MB | 11.26 MB（全混合 48/行） |

### 已知限制
- Fibonacci 在 τ>=86（p>=173，盒子约 4.6λ+）时实球谐合成矩阵条件数退化（κ=Inf），
  属 SH 递归精度/网格分辨极限；p<=163（盒子 <=4λ 层）正常（κ 2-19）。
  更远高阶需要高精度 SH 或更均匀格点（如 HEALPix/等面积格）。
- dataset_generator 的高阶守卫已放宽（不再 throw）。

### 回归
- test_lebedev.jl 9/9、test_mlfma.jl 8/8。

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| Statistics.mean 未导入 | 补 using Statistics | 1 次 |
| Julia 顶层软作用域（for 内对全局赋值） | 包成 main() | 1 次 |
| realSHmatrix 漏 P_1^0（P[0,1] 越界） | 递推补 l=1 特例 pprev=0 | 1 次 |
| nnz(稠密矩阵) MethodError | 改用 count(!iszero, ...) | 1 次 |
