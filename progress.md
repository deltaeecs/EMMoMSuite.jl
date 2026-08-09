# Progress Log: MoM_Lebedev 调研

## 2026-08-08
- [x] 创建 goal（/goal 已生效）。
- [x] 通读 `Lebedev.jl` / `LebedevSortedPoints.jl` / `LVI.jl` / `pinv2interpW.jl` / `dataset_generator.jl` / `test_lebedev.jl`。
- [x] 通读 MLFMA `Interpolation.jl`（GL/Lagrange 两段式插值基线）。
- [x] 确认当前 EMMoMSuite 主 MLFMA 未启用 Lebedev（OctreeBuilder/Precomputations 用 GL+Lagrange）。
- [x] 数值验证：truncation 放大 2.45-3.11x、p=131 死代码、deps/InterpolationWeights 为空、
      现行 W 自身留出集误差 0.79-1.48、行和偏差 278、球谐 3.98e-15、只实部约束、不可复现、
      大阶数内存 4.65 GB。
- [x] 文献对照：He&Kong 2022（方法来源）、SE-MLFMA、FFT 插值、球谐 LS 插值。
- [x] 最终报告。

## 新目标（2026-08-08 续）
- [x] R1: 找到并克隆 github.com/deltaeecs/MoM_Lebedev.jl；确认原始算法与合并版一致；
      旧权重从未进 git；论文 PDF 被 IEEE paywall 挡住（待用户提供）。
- [x] R2: 复现原始管线：13->27 (EFIE 1.475 / 限带 3.6e2)、27->53 (0.558 / 1.1e2)、
      65->131 不可行（4.6GB）；发现训练数据带宽超粗网格可表示范围。
- [x] R3: 结构分析脚本 scripts/lebedev_structure_analysis.jl：
      24 元素八面体群轨道（6/8/12/24）、φ 非等间距、旋转不变性、正交性。
- [x] R4: 新增 src/FastAlgorithms/Lebedev/SHInterp.jl（exact/local/orbit + vectorize）；
      LVI.jl 修复 truncation（*2π bug）与 off-by-one；构造函数支持 method=:sh_*。
- [x] R5: 验证：test_lebedev.jl 9/9；sh_exact 3e-15；轨道压缩与朴素一致；
      demo 脚本 scripts/lebedev_optimized_demo.jl；基准脚本 scripts/lebedev_opt_benchmark.jl。

## 权重矩阵修复（2026-08-08 用户指令）
- [x] 定位根因：粗/细层数据来自不同随机源几何（generate_dataset_on_pkpt 两次调用
      generate_dataset_on_poles，各自内部随机）；vcat(real,imag) 只实部约束；臂长与层失配。
- [x] 修复 dataset_generator.jl（random_source_geometry + evaluate_poles! 共享几何、
      rel_l 缩放、geom 关键字）、pinv2interpW.jl（hcat 全复约束）、SHInterp.jl（degree=Lb、
      interp_weights_auto）、LVI.jl（interpolationCSCMatCal 默认 :sh_auto）。
- [x] 验证：论文 Fig.2 复现（k=8 εi 2.5e-4~1e-3，随点数递减）；
      球谐精确机器精度/确定性；轨道压缩一致；test_lebedev.jl 9/9。
- [x] 新脚本：lebedev_paper_metric_repro.jl / lebedev_paper_repro_v2.jl /
      lebedev_data_spectrum.jl / lebedev_arm_sweep.jl / lebedev_weights_fix_verify.jl。

## MLFMA 集成与基准（2026-08-08）
- [x] 集成：interp_type、OctreeBuilder(interp_method)、Aggregation/Disaggregation 一步分支、
      MLFMAOperator(interp_method, near_range)。
- [x] 真实基准 scripts/lebedev_mlfma_benchmark.jl + scripts/lebedev_mlfma_tradeoff.jl：
      r=1λ / r=2λ PEC 球 EFIE；εq、插值内存、MVM 时间、极点数。
- [x] 结论：修复后训练式 = 最优工程选择（εq==GL、极点-29%、MVM-35%）；
      坏权重 εq 60x 差（MLFMA 层证明修复有效）；稠密球谐在矢量类上受 L+1 耦合限制。
- [x] 回归：test_lebedev.jl 9/9、test_mlfma.jl 全过。

## 新目标（2026-08-08：方案更新 + 效率优化）
- [x] 实现自旋加权球谐（VSH）基与插值：正交性 1e-13、自旋限带机器精度；
      发现 θ/ϕ 基极点奇异残差（εi~0.76），转向笛卡尔方案。
- [x] 新算法：笛卡尔标量 SH 矢量插值（稠密/局部/轨道）：
      27->53 εi=1e-14（免训练确定性）；局部 L_loc=3 46nnz/行 1.9e-4 优于训练式；
      轨道压缩 65->131 构造 887s->1s（100x），与朴素版一致（6e-9）。
- [x] 高阶扩展：65->131 εi=7.1e-5（L_loc=6, 180nnz/行）/ 2.4e-6（L_loc=10）。
- [x] 采样点研究：Fibonacci 同点数 κ 1.3-1.4 vs Lebedev 1.7-2.2、最小角距 2x、
      任意 n；GL 张量积 L=13 κ=2.48 最差。
- [x] 生成器源几何修复（源在盒内）+ 效率优化（O(n*24) 轨道分组、向量化、三元组装配）。
- [x] 脚本：lebedev_cart_verify.jl / lebedev_vsh_verify.jl / lebedev_grid_study.jl；
      test_lebedev.jl 9/9。

## 性能对比（2026-08-08）
- [x] 微基准 + 真实 MLFMA（r=2λ）速度/εq：训练式 k=8 最快（0.239s，比 GL 快 30%，
      εq 0.116 与 GL 一致）；笛卡尔 L_loc=3 ≈GL（0.346s）；L_loc=6 1.5x 慢；稠密 4x 慢且 εq 差。
- [x] kNN 紧支撑实验：省 30% nnz 但 εi 恶化 44x -> cap 球冠为稳健选择。
- [x] 结论：MLFMA 默认用训练式 k=8；笛卡尔局部为确定性/高精度选项，速度成本如实报告。

## 混合权重与 Fibonacci（2026-08-08 用户："按你建议的试试"）
- [x] interp_weights_hybrid：自适应支撑 + 笛卡尔 SH 精确约束 + 数据拟合（KKT），
      模块内确定性 RNG；27->53 εi 7.7e-6（48/行，2.4s），65->131 5.9s；
      MLFMA r=2λ MVM 0.282s（εq 0.116 同训练式）；确定性逐位一致。
- [x] Fibonacci 更少点数测试：同点数 εi 差 2x、减点后差 10-100x -> 局部稀疏场景
      Lebedev 更优；脚本 lebedev_fibonacci_test.jl。
- [x] 脚本 lebedev_hybrid_verify.jl；回归 test_lebedev.jl。

## 4λ 球规模行为（2026-08-08 用户："好的"）
- [x] 修复 LVI 回退 bug（Val(:Lagrange2Step) 不存在 -> 改为 GL 调用）；
      阈值 4M->250k 使中阶对也走混合。
- [x] 4λ（N=3312, 5 层）：εq 与 GL 完全一致 0.6698，MVM 1.217s vs 3.536s（快 2.9x），
      极点 -29%（Lebedev 层），插值 14.7 MB 全混合；稠密对曾使 εq 恶化到 0.975。
- [x] 脚本 lebedev_4lambda_benchmark.jl；回归 test_lebedev 9/9、test_mlfma 8/8。

## 高阶无数据集方案（2026-08-09）
- [x] fibonacci_grid / high_order_nodes 接入 get_t_nodes（p>131 自动 Fibonacci）；
      LbPolesInfo 增加 p 字段，interpolationCSCMatCal 弃用 n2pDict；
      levelIntegralInfoCal 高阶层构建 Fibonacci+等权重极点；dataset_generator 守卫放宽。
- [x] 验证：(101->163) 混合 εi=1.18e-4（48/行，13.1s）；
      4λ 球 level1 8965 点（vs GL 13448，-33%），极点总数 -31.5%，
      εq=0.6698 与 GL 一致，MVM 1.079s 快 18%。
- [x] 已知限制：τ>=86 时 Fibonacci SH 条件数退化（κ=Inf）。
- [x] 回归 test_lebedev 9/9、test_mlfma 8/8。
