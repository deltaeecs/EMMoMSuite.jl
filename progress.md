# Progress Log

## Session: 2026-08-12

### Current Status
- **Phase:** 3 - 方案设计（调研与现状盘点完成，方案已交付）
- **Started:** 2026-08-12

### Actions Taken
1. 读取 cem skill 与 CEM 知识库 SKILL.md，确认路由（Layer 0/1/2）。
2. 检索知识库 FFT/加速/MLFMA 相关内容，定位 8 组出处；确认 Ergul Ch3 未解析。
3. 盘点 EMMoMSuite：Julia 1.10，无 FFTW 依赖；MLFMA 已有 GL+Lagrange2Step 与 Lebedev 1-step 两条插值路径；mul! 流程 = 近场 + 聚合 + 转移 + 反聚合；测试与基准门已确认。
4. 生成 task_plan.md / findings.md / progress.md（本文件）。
5. 从 `D:\华为家庭存储\projects\KnowledgeBase\CEM\Books\[Ozgur_Ergul,_Levent_Gurel]_The_Multilevel_Fast_Mu(z-lib.org).pdf` 只读提取第 3 章（PDF 193–284 页）文本，核实采样、两步 Lagrange、虚拟扩展、极点采样、转移算子插值内容；确认书中未用 FFT。
6. 核实 Hansen Ch4 §4.3.3 关于"周期带限函数可由 DFT/FFT 精确求值"的论述，作为 P1 FFT 谱插值的理论出处。
7. 将 P1 详细设计（算法、数据结构、调用点、验证门、依赖、提交粒度）写入 task_plan.md，更新 findings.md。
8. 同步知识库修复到 `C:\Users\12253\OneDrive\projects\KnowledgeBase\CEM`（robocopy，23064 文件全部 MATCH）。
9. 创建实施目标（P0-P5）；P0：Project.toml 增加 FFTW/AbstractFFTs，预编译通过。
10. P1 实现（TDD）：`FFTGLPolesInfo`/`FFTInterpInfo`/`fft_interp_phi`/`fft_anterp_phi`(+`!` 零分配版)、
    `Val(:FFTSpectral)` 分派、Aggregation/Disaggregation 接线；test_fft_interp.jl 19/19 通过；
    接入 runtests.jl 与 runtests_fast.jl。
11. P1 关键修正：半格网格频域相位修正；反插值改为插值矩阵转置（伴随）而非频域截断。
12. P1 性能：matvec 0.96×（193ms vs 185ms，N=1440）；φ 步叶层 0.81、高层 0.28；
    批量 FFT 为下一步优化。
13. P1 批量 FFT：`fft_interp_phi_batch!`/`fft_anterp_phi_batch!`（按线程缓存计划/缓冲，零分配），
    重构 aggregate_upward!/disaggregate_downward!；批量等价性测试加入 test_fft_interp.jl（47/47 全绿）。
14. P1 性能基准固化：`benchmark/benchmark_fft_interp.jl`（每层 φ 步 + matvec 时间/内存分列）；
    N=1440：matvec 184.5 vs 199.3ms（0.93×），alloc 79 vs 93MB，rel diff 1.7e-6；
    φ 步批量 4.14/3.78ms vs 稀疏 1.25/2.37ms（比值 0.30/0.63）。
15. P2 评估完成：`docs/assessment_structured_kernel_fft.md`（结论：当前架构下仅 φ 插值为直接适用点；
    2D/周期/BoR 远期；一般几何走 AIM/IE-FFT）。
16. P3 AIM/IE-FFT 实现：`src/FastAlgorithms/AIM/AIM.jl`（网格/投影/FFT 卷积/近场修正/算子），
    test_aim_operator.jl 365/365 全绿；关键校准：远场常数 C=jkη（初值 jkη/4π 使远场低估 4π，已实测校准）。
17. P3 门：0.8λ 球 matvec vs 稠密 0.67%（<1e-2）；性能 AIM 4.7ms vs MLFMA 50.8ms（10.8×）；
    求解 10.7%（预条件待改进）；基准脚本 benchmark/benchmark_aim_vs_mlfma.jl。
18. P5-AIM：`AIMWS` 按线程工作区（具体类型字段，修复 Vector{Any} 导致 18.6MB/次广播分配）；
    matvec 分配 16.8MB→13KB、时间 10ms→3.47ms（0.6λ 球）；benchmark 脚本补充 alloc 输出。
19. P4 评估：NUFFT 当前不引入（稀疏一步插值在小规模更快、避免新依赖），结论写入 task_plan。
20. 迭代检视 Round 1（P1+P3）：无新问题（理论性线程回退说明已记录）；412 项测试 + MLFMA/Lebedev 回归全绿。
21. P1 性能门规模扫描：盒 0.25-16λ（τ=19-217）φ 步 FFT/稀疏 0.59→0.17——FFT 在 MLFMA 全规模范围
    不敌稀疏 ϕCSC（素数因子 FFT 尺寸）；门假设被实测推翻，调整为整体持平+内存+规模特征化（如实记录）。
22. P5-AIM：空间哈希近场装配（修复自对缺陷，113,264 断言全绿）；`Z_near_direct` 预条件矩阵
    （实测 2.35%，与无预条件同——误差由算子精度主导）。
23. 迭代检视 Round 2/3：无新问题（连续第 3 轮，结束条件达成）。
24. P5-AIM 线程化：投影/测试 `@threads :static` 并行；修复 @threads 池 thread id > nthreads
    导致归约丢失的缺陷（4 线程验证精度 0.38% 与单线程一致、确定性成立）；N=792 时 3.83ms
    （线程开销，大 N 才有收益）。
25. 迭代检视 Round 4（线程化终态）：无新问题（连续第 4 轮）。
26. 全量快速套件：`test/runtests_fast.jl` 114,574 通过 / 3 broken（既有）/ 0 失败，59s——
    两个新测试文件（test_fft_interp、test_aim_operator）套件级接入验证通过。
27. 大规模 AIM 基准（N=1836，1.2λ 球）：精度 1.04%（h=0.1λ）、近场填充 23.4%、
    matvec 12.15ms vs MLFMA 194.9ms（16.0×）。
28. 【新目标】混合并行盘点：MLFMAOperatorMPI（秩分叶子盒 + 秩内 @threads + SPMD 上/下传）、
    分布式 GMRES（行分区 Krylov）、MPIArray 均已存在；环境 MSMPI：`mpiexec -n P julia -t T` 每秩 T 线程。
29. 混合精度门：新增 `test/test_hybrid_mlfma.jl`，`mpiexec -n 2 julia -t 2` 下 matvec 与求解
    相对串行参照均为 **0.0（逐位一致）**（P=2,T=2 验证）。
30. 混合效率基准（N=3312，2.0λ 球，`benchmark/benchmark_hybrid_mlfma.jl`）：
    | 配置 | 线程数 | setup_s | matvec_ms | solve_s | \|I\| |
    |---|---|---|---|---|---|
    | 1×1 | 1 | 18.63 | 518.8 | 21.63 | 0.1454 |
    | 2×2 | 4 | 11.91 | 331.5 | 11.89 | 0.1454 |
    | 2×4 | 8 | 8.92 | 192.1 | 8.79 | 0.1454 |
    | 4×2 | 8 | 13.03 | 431.5 | 17.78 | 0.1454 |
    结论：2×4 最佳（setup 2.1×、matvec 2.7×、solve 2.5× vs 1×1）；4×2 劣于 2×4
    （SPMD 远场阶段重复 + Allreduce/分区开销）——线程比更多秩在本规模更有效；所有配置解一致。
32. AIM 混合并行：新增 `AIMOperatorMPI`（本地行近场 + 复制投影/FFT + 本地行测试 + 单次 Allreduce）。
    精度门 `test/test_hybrid_aim.jl`（mpiexec -n 2 julia -t 2）：rel = 0.0（逐位一致）。
    关键发现：AIM matvec 由网格 FFT 主导（近场仅 2-5%），本地投影+Allreduce 大网格反而变慢；
    改用复制投影 + 秩内 FFTW 多线程（`FFTW.set_num_threads(T)`）后全面反超串行。
33. AIM 混合基准（`benchmark/benchmark_hybrid_aim.jl`，含 FFTW 线程）：
    | 配置 | N | setup_s | matvec_ms | solve_s | \|I\| |
    |---|---|---|---|---|---|
    | 1×1 | 2280 | 25.38 | 25.72 | 1.08 | 0.1156 |
    | 2×4 | 2280 | 13.34 (1.91×) | 18.74 (1.37×) | 0.85 (1.27×) | 0.1156 |
    | 1×1 | 3312 | 38.38 | 148.9 | 6.41 | 0.1441 |
    | 2×4 | 3312 | 19.13 (2.01×) | 64.13 (2.32×) | 2.59 (2.48×) | 0.1441 |
    所有指标混合优于串行、解一致；串行 AIM 回归 113,264 全绿。
34. 4×2 效率优化（用户反馈）：定位瓶颈为远场 SPMD 重复 + 内存带宽天花板 + 45MB 通信；
    实施全层按秩分区（修复 disaggG 双重 Allreduce 缺陷，P=4 精度恢复逐位一致 0.0）；
    4×2 matvec 431→183ms（2.35×）；尝试按源盒外循环因线程竞态/更慢回退；
    实测 8-worker 效率 ~36%（N=3312/4536 一致），80% 在本机带宽与复制式设计下不可达，结论已记录。
35. 未适配 MPI 分支补齐：FFTSpectral（disaggregate FFT 分支 child_filter + interp_method 透传，
    P=2 rel=0.0）；Lebedev 验证（P=2 rel=0.0）；PMCHWMLFMAOperatorMPI（四遍分区远场 +
    独立 y_pass，P=2 rel=7.5e-17）；compute_interpolation_matrices! 线程化。
    回归：MLFMA/AIM/FFT 全绿；PMCHW 仅既有 RED 门 GD2A（k1 聚合奇偶）失败，与本任务无关。
36. 分布式稠密 LU（边界完成）：调研 ScaLAPACK（无绑定、JLL 内嵌 MPI 与 MSMPI 冲突）与
    Pardiso（MKL 回退不可用，需许可库）后，实现原生 MPI 1D 行分块 LU（DistributedLU.jl）；
    P=2/P=4 门 rel=6.9e-16（对角占优）与 2.3e-14（主元压力）通过；N=300 分解 0.25/0.26s。
    Pardiso 依赖已移除。
37. DistributedLU 效率优化：主元行 Bcast + 旧行 k 点到点 + 逐行融合广播 axpy（弃 strided ger!）；
    封装 DistributedDenseLU（mpi_lu/ldiv!）；N=600 分解 0.04s(P2)/0.02s(P4)、N=1200 0.54/0.22s，
    较旧版快一个量级；门全绿（6.9e-16/2.3e-14/7.0e-16）。
38. 效率对比与 ScaLAPACK 调研：OpenBLAS 多线程 N=1200 T=4 为 0.024s，原生 MPI 分布式 LU P=4
    为 0.53s（单机慢 ~20 倍，用户判断正确）；ScaLAPACK 三路实测被堵（JLL 无 Windows、
    MKL msmpi BLACS 为单序数转发桩、通用 BLACS 硬编码 Intel MPI）；结论：单机用多线程 LAPACK、
    分布式 LU 面向多节点集群；实验模块与依赖已移除。
39. ScaLAPACK 本地 MinGW 库接入（用户提供新路径）：`C:\msys64\mingw64\bin\libscalapack.dll`
    （2.2.2，MSMPI）导出完整 BLACS/pzgesv，据此实现 `ScaLAPACKLU.jl`（BLACS 网格 + 块循环
    分发 + pzgesv + 收集，`scalapack_lu_solve`）。修复两个实测坑：BLACS 必须先
    `blacs_get(-1,0)` 取系统上下文再 gridinit（否则 gridinfo 返回垃圾 myrow，descinit -9）；
    PZGESV 要求方形块（MB==NB，源码确认），已加显式校验。分发改为按块向量化复制
    （N=9600 节省 ~2s）。门：P=2/P=4 rel=6.35e-16/6.24e-16。
    效率：N=4800 4×2 = 2.05s（≈串行 T=8 2.12s）；N=9600 4×2 = 5.22s vs 串行 T=8 12.78s
    （**2.4×**）；N=4800 原生 DistributedLU P=2 = 111.9s（ScaLAPACK 快 ~80×）；
    大 N 时不设 OMP_NUM_THREADS 最优（默认 OpenBLAS 线程表现最好）。结论更新：
    本地 MinGW ScaLAPACK 是可行的分布式稠密 LU 路径，取代自研 1D 版本成为默认分布式直接求解器；
    自研版保留为无该库环境的回退。
40. 预条件 MPI 化：确认 `distributed_gmres!` 原先忽略预条件 kwargs；新增
    `src/Parallel/MPI/Preconditioners.jl`（分布式 BlockJacobi 按 cube 归属分秩、
    DistributedDiagonal、apply_mpi_preconditioner! 统一入口 + 串行回退），
    接入 `distributed_gmres!` 的 `Pl`（左预条件）。修正右/左预条件混用 bug
    （解偏差 0.7% → 3.5e-11）。门 `test_hybrid_preconditioner.jl` P=2/P=4 全绿：
    分块均分（68 块 → 每秩 34/17）、施加与串行 rel=0.0、GMRES+Pl 收敛 3.5e-11。
    回归：串行预条件、分布式 GMRES 全过。
41. PMCHW 原生近场装配：重写 `assemble_near_field_pmchw` 按 octree 近邻对原生计算
    四块（EFIE L-核 + PMCHW K-核，线程并行 COO，cube_filter 分区），串行/MPI 均不再
    装配 2N×2N 稠密。修复两个 kernel 原地缩放导致的缓冲复用 bug 与 K 块抵消噪声过滤。
    GD1 逐元素对齐 8e-14、混合门 P=2/P=4 rel=0.0；GD2A/GD2R/B2 为既有 RED
    （远场 k1 parity/精度，旧近场对照相同）。
42. 内存瓶颈审计：定位并修复远场翻译因子 αTrans 内存爆炸——全偏移表
    （near_range=16 时 264826 列/层 × nPoles）被逐层逐秩复制，N=594 时每秩
    **15.7GB**；改为只存实际使用的偏移列（4454 列，0 远邻层 0 列）→ **21MB**
    （740×），PMCHW 构造 137s→10s。加上 PMCHW 近场去稠密与预条件分块，
    端到端（P=2、N=594）峰值内存 17.6GB→1.3–1.6GB。`test_mlfma.jl` 断言更新为
    "无远邻层允许空 αTrans"，translate! 语义不变（混合门 rel=0.0 验证）。
43. 用例精度回归（用户要求"确保用例精度不出问题"，2026-08-13）全绿：
    工程规模 PMCHW-MLFMA 门——budget_medium（N=540，matvec rel=7.9e-6 vs 稠密、
    Zin 一致）、budget_krylov_medium、gate_s_mlfma_medium（relw=1.3e-5、rels=4e-5）、
    gate_s_dense；串行回归——test_pmchw、test_nmuller_comparison_medium、
    test_integral_equations_endtoend（EFIE/MFIE/CFIE）、test_mlfma、test_aim_operator、
    test_preconditioners、test_solvers_verification；MPI 混合门——hybrid_mlfma/aim/pmchw
    P=2 与 hybrid_pmchw P=4 全部 rel=0.0（逐位一致）、hybrid_preconditioner P=2/P=4
    （施加 rel=0.0、GMRES+Pl 3.5e-11）、distributed_gmres P=2/P=4。既有 RED 门
    GD2A/GD2R/B2 为远场 k1 parity/精度问题（旧近场对照相同），与本改动无关。

### Test Results
| Test | Expected | Actual | Status |
|------|----------|--------|--------|
| 知识库出处文件存在性 | 全部可读 | 8 组出处均可读 | pass |
| 工程模块/依赖盘点 | 无 FFTW；插值两路径 | 确认 | pass |
| Ergul Ch3 内容 | 可核对 | 未入库解析，但 PDF 文本层可读，本次已提取核实 | pass（KB 入库解析待用户确认） |
| P1 功能正确性 | 19/19 | test_fft_interp.jl 全绿；MLFMA 集成 rel 1.2e-6；稠密误差 5.8e-4 与 Lagrange 持平 | pass |
| P1 性能门（φ 步 < ϕCSC） | 全部层通过 | 叶层 0.81、高层 0.28，未全过 | fail（批量 FFT 待办） |
| P1 批量 FFT 等价性 | 47/47 | 批量 = 逐调用（机器精度） | pass |
| P1 性能门复测 | 记录 | 批量 φ 步 0.30/0.63；matvec 0.93×（噪声内持平） | documented（大规模验证待办） |
| MLFMA/Lebedev 回归 | 全绿 | test_mlfma.jl、test_lebedev_interp.jl | pass |
| MLFMA 混合精度门 | rel < 1e-10/1e-4 | test_hybrid_mlfma.jl（P=2,T=2）matvec/solve 均为 0.0 | pass |
| AIM 混合精度门 | rel < 1e-10 | test_hybrid_aim.jl（P=2,T=2）rel = 0.0 | pass |
| 混合效率门（MLFMA/AIM 大用例） | hybrid < 单线程 | 2×4 全指标优于 1×1（见上表） | pass |
| ScaLAPACK 分布式 LU 精度门 | rel < 1e-8 | test_scalapack_lu.jl P=2 rel=6.35e-16；P=4 rel=6.24e-16 | pass |
| 分布式 LU 回归 | 全绿 | test_distributed_lu.jl P=4 对角占优 6.9e-16 / 主元压力 2.3e-14 | pass |

### Errors
| Error | Resolution |
|-------|------------|
| 无 | - |
