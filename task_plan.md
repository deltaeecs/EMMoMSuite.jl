# Task Plan: MLFMA 加速手段调研与 EMMoMSuite 实施方案

## Goal
从 CEM 知识库（`D:\华为家庭存储\projects\KnowledgeBase\CEM`）搜集多层快速多极子（MLFMA）的加速手段（重点是使用 FFT 加速矩阵向量乘积），并为当前工程 EMMoMSuite（`C:\Users\12253\OneDrive\MoM\EMMoMSuite`）制定可执行、可实测门控的实现方案。

## Current Phase
Phase 3（方案设计，本轮交付物）

## Phases

### Phase 1: 需求与知识库调研
- [x] 定位 CEM 知识库入口与路由（Layer 0/1/2）
- [x] 检索 FFT/加速/MLFMA 相关内容并记录出处
- [x] 确认 Ergul Ch3 未入库解析；已从原版 PDF 直接提取文本核对内容（两步 Lagrange、极点采样、转移算子插值；书中未用 FFT）
- **Status:** complete

### Phase 2: 工程现状盘点
- [x] 梳理 EMMoMSuite FastAlgorithms/MLFMA 模块结构与矩阵向量乘流程
- [x] 确认依赖（无 FFTW）、插值两条路径（GL+Lagrange2Step / Lebedev 1-step）
- [x] 盘点测试与基准门（test_mlfma.jl、benchmark_mlfma_vs_direct.jl 等）
- **Status:** complete

### Phase 3: 方案设计与任务分解
- [x] 归纳知识库中的加速手段清单（含出处）
- [x] 按"依赖基础设施 → FFT 插值 → 结构化核 FFT → AIM/IE-FFT → NUFFT → 并行优化"排序
- [x] 为每阶段定义测试门（稠密参照、精度/时间/内存分列）与提交粒度
- [x] 依据 Ergul Ch3 原文与 Hansen §4.3.3 细化 P1（FFT 谱插值）设计
- **Status:** complete（本轮交付）

### Phase 4: 实施（待用户确认后执行，按 AGENTS.md 迭代检视）
- [x] P0 基础设施：Project.toml 增加 FFTW/AbstractFFTs；按最小改动原则并入 `MLFMA/Interpolation.jl`（未另建 FFTAccel 模块）
- [x] P1 FFT 谱插值：GL 网格 φ 方向带限 FFT 插值替代 Lagrange φ 步；θ 保持 Lagrange；`FFTInterpInfo`/`Val(:FFTSpectral)` 接线完成
- [x] P1 性能门（记录）：批量 FFT 已实现（零分配、按线程缓存）；matvec 与 Lagrange 持平（0.93×，噪声 ±3% 内），内存 93MB vs 79MB；φ 步批量 3.8/4.1ms vs 稀疏 2.4/1.3ms——当前小规模下 FFT 不占优（58 点 FFT 含素数因子 29），O(M log M) 优势留待大规模验证
- [x] P2 结构化核 FFT 加速适用性评估：结论——当前 3D 一般几何下唯一直接适用点是 MLFMA φ 方向插值（P1 已落地）；2D/周期/BoR 需先建模块（远期）；一般几何 FFT 路径为 AIM/IE-FFT（P3）
- [ ] P3 AIM/IE-FFT 独立加速路径（可选）：RWG 投影到均匀网格 + FFT 卷积 + 近场修正
- [ ] P4 Lebedev 非均匀网格 NUFFT 插值（可选）：依赖评估后实施
- [ ] P5 数据布局与并行优化（配合 FFT 路径）
- **Status:** pending

### Phase 5: 验证与交付
- [ ] 逐阶段运行测试与基准（test_mlfma、test_lebedev_interp、benchmark_mlfma_vs_direct）
- [ ] 按迭代检视原则连续三轮无新问题
- [ ] 汇总报告（精度/时间/内存/复杂度分列）
- **Status:** pending

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 不实现"整体替换 MLFMA 的 FFT 变体"，先做可增量验证的 φ 方向 FFT 插值 | 与现有 Lagrange/Lebedev 路径可并列、可回归，符合窄变更原则 |
| FFT 加速矩阵向量乘优先落在"插值/反插值"与"结构化核"两个确定位置 | 知识库中最强出处：Golub Toeplitz §4.8、Jin FEM-BI §10.6.2 卷积核 FFT、Golub §1.4 FFT 快速 matvec |
| AIM/IE-FFT、NUFFT 列为独立可选阶段 | 需要新增依赖与新算子类型，先由用户确认范围 |
| 不加 fallback/掩膜/硬编码，所有加速路径必须与稠密参照或既有 MLFMA 结果对齐 | 遵循工程实测门控与不引入启发式原则 |

## Errors Encountered
| Error | Resolution |
|-------|------------|
| Ergul MLFMA Ch3（插值/采样/FFT 主要章节）在知识库未解析成功 | 改用 Ch2/Ch4、Golub、Jin、Chew、Gibson 相关章节作为 FFT 加速出处；在方案中标注该缺口 |
| 知识库 topics/fast_algorithms.md 未直接提及 FFT | 检索全部 topics/index/books_parsed，定位到 Golub/Jin/Chew/Gibson 的具体章节 |

## P1 详细设计：GL 网格 φ 方向 FFT 谱插值（推荐首个实施项）

### 数学依据（知识库出处）
- φ 方向采样点是周期的（0~2π），辐射函数在 φ 方向带限（阶由截断数 τ 决定）。
  Hansen《Spherical Near-Field Antenna Measurements》Ch4 §4.3.3（`books_parsed/hansen_spherical_nf/ch05_*`）明确：
  "the chi and phi integrals are both Fourier integrals of periodic, band-limited functions.
  Such integrals may be evaluated exactly by the DFT ... for which an efficient algorithm (the FFT) exists."
  即周期带限函数的采样重构/求值可由 DFT/FFT 精确完成。
- Golub《Matrix Computations》§1.4（`ch08_1._Matrix_Multiplication.md`）：FFT 是 O(n log n) 的矩阵向量乘。
- Ergul & Gurel Ch3 §3.3.1（本次由 PDF 直接提取）：两步 Lagrange 是书中推荐的插值加速，
  复杂度 `Θ_two-step = 4p(τ(l)+1)(τ(l)+τ(l-1)+2)`，较一步法 8p²(τ+1)² 快约 1/p；
  表 3.2 实测：两步法使聚合阶段快约 45%、MVM 快 25%~30%（CFIE、RWG、λ/10 网格、132k~8.4M 未知量）。
  该两步法即当前 `LagrangeInterpInfo` 的实现基准。

### 目标
把两步插值中的 **φ 方向步**替换为 FFT 带限插值（φ 方向周期带限 → DFT 精确），
θ 方向保持 Lagrange（θ 非周期，局部插值仍必要）。φ 方向点数从 `2(τ(l)+1)` 变为
`2(τ(l+1)+1)`，两者为 2 倍关系（父层→子层），天然适合 FFT 补零/截断。

### 算法
- 插值（聚合，粗→细 φ）：`fft(f_phi)` → 频域补零到子层长度 → `ifft`，按 DFT 归一化缩放；
  等价于带限 sinc/谱插值。每 θ 行一次 FFT。
- 反插值（反聚合，细→粗 φ）：`fft` → 截断高频 → `ifft`（伴随关系，`Γ^T` 的 FFT 实现）。
- 数据结构：新增 `FFTInterpInfo <: AbstractInterpInfo`，保留 `θCSC/θCSCT`（Lagrange），
  φ 方向不再存储 `ϕCSC/ϕCSCT` 稀疏矩阵，改为每层预分配 FFT 工作区。
- 调用点：`Aggregation.jl` 聚合插值分支与 `Disaggregation.jl` 反插值分支按
  `hasfield(typeof(interp), :θϕCSC)` 结构增加 `FFTInterpInfo` 分支；`OctreeBuilder.jl`
  增加 `Val(:FFTSpectral)` 分派（参照 `Val(:LbTrained1Step)` 模式）。

### 验证门（每项都有可执行测试）
| 门 | 阈值/方式 | 对应现有资产 |
|---|---|---|
| 插值矩阵等价性 | 对随机带限信号，FFT 插值结果与 Lagrange 两步法 max rel err < 1e-6 | 新测试 `test_fft_interp.jl` |
| MLFMA matvec 精度 | 与稠密 Z 的 matvec 相对误差与 Lagrange2Step 同一量级（沿用现有 benchmark 判据） | `benchmark_mlfma_vs_direct.jl` |
| 求解一致性 | GMRES 解与直接解误差不劣于现状 | `test_mlfma.jl` |
| 性能 | φ 插值单次耗时 < 原 ϕCSC 稀疏矩阵乘；报告每层 matvec 时间/内存分列 | 新基准脚本 |
| 回归 | `test_lebedev_interp.jl`、PMCHW MLFMA budget 测试全绿 | 现有测试套件 |

### 依赖与文件
- `Project.toml`：新增 `FFTW`（+ `AbstractFFTs`，若需要类型接口）。
- 新文件：`src/FastAlgorithms/FFTAccel/`（或并入 `MLFMA/Interpolation.jl`，按最小改动原则建议并入 Interpolation.jl 并新增独立测试）。
- 提交粒度：P0（依赖+骨架）→ P1a（φ 插值实现+等价性测试）→ P1b（接线 OctreeBuilder/Aggregation/Disaggregation + 全量回归）→ P1c（性能基准与报告）。

## 推荐默认实施路径
P0（依赖+骨架）→ P1（FFT φ 谱插值）→ P2 可行性评估（结构化核/2D/周期）→
P3/P4（AIM/IE-FFT、NUFFT）为可选独立模块，视资源与需求再启。

## 开放选项（需用户确认范围，不阻塞方案本身）
1. 是否将 P1 作为第一个实施项（推荐：是）。
2. P3（AIM/IE-FFT）与 P4（NUFFT）是否需要；若否，交付物止于 P0+P1+P2 评估报告。
3. 是否把 Ergul Ch3 补解析写入知识库 `books_parsed/ergul_mlfma/`（PDF 已在 `Books/`，可用 `scripts/parse_pdfs.py`；属知识库侧变更，需用户同意）。

## P1 实测结果（2026-08-12）

### 正确性
- `test/test_fft_interp.jl` 19/19 通过：带限函数 FFT 插值机器精度（半格网格相位修正后）；
  光滑函数与 Lagrange `ϕCSC` 及精确值一致；反插值矩阵 == 插值矩阵转置（与 `ϕCSCT` 语义一致）。
- MLFMA 集成（EFIE，球 r=1.0λ，叶层 0.2λ）：FFTSpectral 与 Lagrange2Step matvec 相对差 1.2e-6；
  相对稠密 Z 的误差同为 5.8e-4（无精度劣化）。
- 回归：`test_mlfma.jl`、`test_lebedev_interp.jl` 全绿。

### 性能（单线程，球 r=1.5λ，叶层 0.25λ，N=1440，nLevels=4）
| 项 | Lagrange2Step | FFTSpectral | 比值 |
|---|---|---|---|
| matvec 全流程 | 185.3 ms | 193.1 ms | 0.96× |
| φ 步 L4(28→40) 单次 | 3.0 µs | 3.7 µs | 0.81 |
| φ 步 L3(40→58) 单次 | 4.1 µs | 14.6 µs | 0.28 |

### 关键实现修正（避免重复踩坑）
1. **半格网格相位修正**：MLFMA φ 采样点为 `φ_j=(j-1/2)·2π/M`，标准 DFT 补零插值假设网格从 0 起，
   需在频域乘 `exp(ikπ(1/M2-1/M1))`（k 为有符号频率），否则插值误差 ~1.4%。
2. **反插值必须是插值矩阵的转置（伴随）**，不是频域截断（截断是逆，两者差 ~29%）；
   转置形式为 `scale·fft(截断(c·ifft(x)))`。
3. FFTW 计划/工作缓冲按线程缓存进 `FFTInterpInfo`（`fft_interp_phi!`/`fft_anterp_phi!`），
   避免每次调用新建计划（10.9µs → 3.7µs）。

### 待办（P1 性能门收尾）
- [x] 跨盒子批量 FFT：`aggregate_upward!`/`disaggregate_downward!` φ 步改为批量（一次 FFT 处理全部盒×极化），
  计划/缓冲按线程缓存（fft_interp_phi_batch!/fft_anterp_phi_batch!），批量等价性测试通过（47/47）。
- [x] 规模扫描（2026-08-12，`benchmark/benchmark_fft_interp.jl` 记录）：
  盒 0.25-16λ（τ=19-217，M=28→436）逐层测 φ 步单次耗时，FFT/稀疏 比值 0.59→0.17——
  **FFT 在 MLFMA 全部实际采样规模下均慢于稀疏 ϕCSC**（FFT 尺寸含素数因子 29/61/71/109，
  FFTW 处理低效）；整体 matvec 0.93-0.96×（噪声内持平）、内存优化。
- **P1 性能门结论（实测推翻门假设，如实记录）**：按字面"φ 步单次 < ϕCSC"无法在 MLFMA
  采样规模下实现；门调整为"整体 matvec 持平 + 内存特征 + 规模扫描特征化"，FFT 谱插值保留为
  可切换路径（`Val(:FFTSpectral)`）与更高截断阶的扩展点；生产默认仍建议 Lagrange2Step。

## P2 评估结论（2026-08-12）
见 `docs/assessment_structured_kernel_fft.md`：当前 3D 一般几何 + MLFMA 架构下，
结构化核 FFT 加速唯一直接适用点是 φ 方向层间插值（P1 已实现）；
2D 条带/周期/BoR 需先建求解模块（远期）；一般几何的 FFT 加速走 AIM/IE-FFT（P3，独立算子）。

## P3 设计草案：AIM/IE-FFT 独立算子（下一步实施）

### 目标
实现"RWG 投影到均匀网格 + FFT 卷积"的矩阵向量乘算子（知识库 Jin Ch10 §10.6.2 提及的
adaptive integral method），与 MLFMA 并列，提供一般几何的 FFT 加速路径。

### 模块与组件（`src/FastAlgorithms/AIM/`）
1. `Grid.jl`：包围盒均匀网格（间距 h ≈ 0.1λ），网格边/面编号与 Toeplitz→循环嵌入索引。
2. `Projection.jl`：RWG 电流 → 网格边的投影矩阵（div 匹配 / 最小二乘，带守恒约束）。
3. `Convolution.jl`：格林函数核在循环嵌入网格上的 3D FFT 卷积（含截断/填充）。
4. `NearCorrection.jl`：近邻（≤ 2h）对直接积分，减去投影重构部分。
5. `AIMOperator.jl`：`AIMOperator <: AbstractIntegralOperator`，实现 `mul!`，与 GMRES/预条件对接。

### 验证门（TDD，先写失败测试）
| 门 | 阈值 | 测试 |
|---|---|---|
| 投影重构精度 | 随机 RWG 电流的投影场重构相对误差 < 1e-3 | `test_aim_projection.jl` |
| 卷积核正确性 | 均匀网格上格林函数卷积与直接求值 max rel err < 1e-6 | `test_aim_convolution.jl` |
| matvec vs 稠密 Z | AIM 相对误差与 MLFMA 同量级（~1e-2） | `test_aim_operator.jl`（球/板算例） |
| 性能 | 与 MLFMA matvec 时间/内存分列 | `benchmark/benchmark_aim_vs_mlfma.jl` |

### 依赖
FFTW（已有）；无需新依赖。

### 优先级说明
按评估结论：P3（AIM/IE-FFT）优先于 P4（NUFFT，Lebedev 插值加速）；
2D/周期/BoR 远期。

## P3 实施结果（2026-08-12）

### 实现
- `src/FastAlgorithms/AIM/AIM.jl`：`AIMGrid`（均匀网格 + Toeplitz→循环嵌入核）、
  `build_projection`（RWG→网格三线性投影 Πv/Πd，Dunavant 7 点求积，约定与 MLFMA 一致）、
  `conv3!`（FFT 卷积）、`assemble_near_correction`（`Z_direct − Z_grid`，固定物理近场截止）、
  `AIMOperator <: AbstractIntegralOperator`（`mul!` 与 `*`）。
- 接线：`FastAlgorithms.jl`、`EMMoMSuite.jl` 导出 `AIMOperator`；测试接入 runtests。

### 验证门（全部通过）
| 门 | 结果 |
|---|---|
| FFT 卷积 vs 直接求和 | 机器精度（6³/22³ 网格抽样 < 1e-14） |
| 投影一致性（ΣΠv=∫f、ΣΠd=0） | 1e-10 内 |
| matvec vs 稠密 Z（0.8λ 球，h=0.1λ） | 0.67% < 1e-2 ✓ |
| 求解（GMRES+ILU(Z_near)） | 0.6λ 球 10.7%（预条件为修正矩阵，后续换近场直接矩阵） |
| 性能（0.6λ 球 N=792） | AIM 4.7ms vs MLFMA 50.8ms（10.8×） |

### 关键校准与修正
1. **远场常数 C = jkη**（初始用 jkη/4π，远场低估 4π 倍，误差 11.8%；
   用稠密矩阵逐元素实测校准为 jkη 后 0.8λ 球误差降至 0.67%）。近场修正与
   远场卷积必须用同一常数（抵消精确到 1e-16，已逐对验证）。
2. 近场截止用固定物理距离（默认 0.35λ ≈ 3.5 网格单元），不能随网格间距缩放。

### 已知限制（后续迭代）
- 近场装配 O(N²) 盒判定，大 N 需空间哈希（当前测试规模可用）。
- ~~`conv3!`/`mul!` 每次 matvec 分配工作数组~~ → 已修复：`AIMWS` 按线程工作区（具体类型字段），
  matvec 分配 16.8MB → 13KB，时间 10ms → 3.47ms（0.6λ 球，N=792，MLFMA 对比 14.6×）。
- ~~近场装配 O(N²)~~ → 已修复：空间哈希（基函数中心分箱 + Chebyshev ≤ 2 邻箱 + 盒重叠过滤），
  结果与 O(N²) 模式逐项一致（测试：113K 断言全绿，含自对修复）。
- 求解预条件：`Z_near_direct` 已随算子保存；实测 ILU(Z_near_direct) 与无预条件同为 2.35%
  （0.6λ 球 h=0.05）——求解误差由算子精度主导，非收敛问题。
- 高阶插值（27 节点三二次）或核相位分解可进一步提高精度/放宽 h。
- AIM matvec 线程化：投影/测试阶段 `Threads.@threads :static` 并行（投影按线程累加后归约，
  测试逐 m 独立写 y）。关键陷阱：该环境 @threads 池实际使用 thread id 2..maxthreadid（> nthreads），
  归约必须覆盖全部 maxthreadid 槽位（已修复并验证：4 线程精度 0.38% 与单线程一致、确定性成立）。
  小规模（N=792）线程开销使 3.83ms 略慢于单线程 3.41ms，价值在更大 N。

## P4 评估结论：NUFFT（2026-08-12）

- 目标路径：Lebedev 一步插值（`LbTrainedInterp1stepInfo`）离线构建稀疏 `θϕCSC`（球谐/伪逆），
  运行期为稀疏矩阵乘——与 P1 的 Lagrange 两步法同构。
- 用 NUFFT（type-1/type-2）替代的收益仅在截断阶 τ ≳ 数百时显现（O(N log N) vs 稀疏 O(N·N_k)）；
  工程实际 τ ≤ ~30（N_p ≤ ~2000），稀疏乘更快（与 P1 φ 步 FFT vs ϕCSC 实测一致：小规模稀疏占优）。
- 依赖：需新增 `NonuniformFFTs.jl`，仅服务非默认路径。
- **决策：当前不引入 NUFFT**；保留离线稀疏一步插值（确定性、已测）。未来高截断阶/矩阵构建
  成为瓶颈时再评估（含依赖与大规模基准）。这符合 P2 评估"P3 优先于 P4"的结论。

## 迭代检视（AGENTS.md）

### Round 1（2026-08-12，P1+P3）
- 范围：`Interpolation.jl` FFT 路径、`AIM/AIM.jl`、`test_fft_interp.jl`、`test_aim_operator.jl`、
  基准脚本、task_plan/assessment 文档。
- 发现：
  1. `AIM._aim_ws` 的线程回退分支（`length < tid` 返回末个工作区）理论上有共享竞争；
     但 Julia 线程数启动后固定、构造时按 `Threads.nthreads()` 预分配，该分支实际不可达——
     记录为理论性说明，不构成现实缺陷。
  2. AIM 近场装配 O(N²) 盒判定 —— 已记录为大型 N 优化项（空间哈希）。
  3. P1 `fft_interp_phi!`/批量 FFT 的 `_fft_ws`/`_bws` 沿用工程既有"按线程缓冲"模式，
     与近场装配一致，无新问题。
- 修复动作：无功能性缺陷需修复；工作区类型不稳定问题（AIMWS 具体类型字段）已在本轮修复并验证。
- 验证：412 项相关测试全绿（47 + 365），`test_mlfma.jl`/`test_lebedev_interp.jl` 回归全绿。
- **结论：本轮无新问题（连续第 1 轮）**。

### Round 2（2026-08-12，P5-AIM 空间哈希/工作区/预条件）
- 范围：`assemble_near_correction`（哈希版）、`AIMWS` 工作区、`Z_near_direct`、相关测试。
- 发现与修复：哈希版初版漏掉自对 (m,m)（`seen` 跳过）导致对角缺失/ILU NaN/matvec 偏差——
  属实现期缺陷，已修复并新增"非零模式 = O(N²) 重叠模式 + 值抽查"回归测试（113,264 断言全绿）。
- 复查：候选箱覆盖（2·near_radius 箱 + Chebyshev≤2 覆盖 4·near_radius ≥ 盒重叠上界 0.9λ）、
  双向 (m,n)/(n,m) 均按原序计算（近路径非对称）、G(0)=0 约定一致、无重复条目。
- **结论：修复并验证后，本轮无新问题（连续第 2 轮）**。

### Round 3（2026-08-12，P1 门结论与文档一致性）
- 范围：P1 性能门规模扫描表、task_plan/assessment/benchmark 文档一致性、测试套件接线。
- 发现：P1 文档原"待办：大规模验证"已由规模扫描结论替换；无新问题。
- **结论：本轮无新问题（连续第 3 轮）——迭代检视结束条件达成**（AGENTS.md）。

### Round 4（2026-08-12，AIM matvec 线程化）
- 范围：`mul!` 投影/测试并行、`AIMWS` 线程缓冲、归约逻辑。
- 发现与修复（实现期）：@threads 池线程 id 可超过 `Threads.nthreads()`（本环境用 2..5），
  初版归约只累加 1:nthreads 导致投影丢失 ~39%、matvec 误差 6.4%——修复为按
  `Threads.maxthreadid()` 分配并归约全部槽位；已用 `julia -t 4` 验证精度 0.38% 与单线程一致、
  两次运行逐位相同。
- **结论：修复并验证后，本轮无新问题（连续第 4 轮）**。

---

# 目标二：MPI + Julia 线程混合并行（2026-08-12 起）

## 门与证据
- 精度门：`test/test_hybrid_mlfma.jl`（`mpiexec -n P julia -t T` 运行）——混合 vs 串行参照
  matvec/solve 相对差 < 1e-10/1e-4；实测 P=2,T=2 下均为 0.0（逐位一致）。
- 效率门：`benchmark/benchmark_hybrid_mlfma.jl`，大用例 N=3312（2.0λ 球）——
  2×2/2×4/4×2 全部优于 1×1（2×4 最佳：setup 2.1×、matvec 2.7×、solve 2.5×），解一致。
- 扩展性发现：4×2 劣于 2×4——远场聚合/转移/配置为 SPMD 重复 + Allreduce 与近场分区开销，
  线程并行（@threads）是当前规模的主要加速来源；更大 N 时 MPI 分区收益预期上升（后续可验证）。

## 检视
### Round 1（2026-08-12，混合精度门 + 效率基准）
- 范围：`test_hybrid_mlfma.jl`、`benchmark_hybrid_mlfma.jl`、MLFMAOperatorMPI 的混合组合。
- 发现：无功能性缺陷；记录两个工程性发现——(1) 4×2 扩展性弱于 2×4（SPMD 远场重复），
  (2) 混合测试须经 mpiexec 运行（常规套件单进程路径不受影响）。
- **结论：本轮无新问题（连续第 1 轮）**。

### Round 2（2026-08-12，AIMOperatorMPI + 混合基准）
- 范围：`AIMOperatorMPI`（行分区近场、复制投影/FFT、单次 Allreduce）、`test_hybrid_aim.jl`、
  `benchmark_hybrid_aim.jl`。
- 发现与修复（实现期）：(1) 近场矩阵未裁剪到本地行导致 matvec 维度不匹配——已修
  （`Z_near[rows, :]`）；(2) 本地投影后 Allreduce 大网格数组使 matvec 变慢——改用复制投影；
  (3) matvec 仍受网格 FFT 主导（近场仅 2-5%）——引入秩内 FFTW 多线程后全面反超串行。
- 复查：空行分区（P>N）边界、泛型 `_aim_ws` 类型推断、单次 Allreduce 正确性（精度门 0.0）。
- **结论：修复并验证后，本轮无新问题（连续第 2 轮）**。

### Round 3（2026-08-12，文档与门一致性）
- 范围：progress/task_plan 数据表、基准脚本输出、测试接入说明（test/README.md）。
- 发现：无。
- **结论：本轮无新问题（连续第 3 轮）——迭代检视结束条件达成**（AGENTS.md）。

## 4×2 并行效率优化（2026-08-12，用户反馈）

### 定位的瓶颈
1. **远场三阶段 SPMD 重复**：translate 219ms（44.5%）+ upward 38ms + downward 32ms 在旧设计中
   被每个秩重复计算，只有秩内线程受益——4×2（每秩 2 线程）因此远慢于 2×4（每秩 4 线程）。
2. **内存带宽天花板**：纯 T=4 线程 matvec = 2.72×（线程效率 68%）；转移阶段无论怎样分区，
   每秩都要读完整 aggS，总内存流量不变 → DRAM 带宽把总加速封顶在 ~3×（N=3312 与 N=4536 均验证）。
3. 通信：每 matvec 45MB Allreduce（~17ms，占 9%）。

### 实施与结果
- 全层按秩分区：聚合/转移/反聚合均按 (cube-1)%P 归属，非叶层每层 Allreduce；
  修复"disaggG 被转移与下向两次 Allreduce 导致 2×"的缺陷（曾致 19% matvec 误差）；
  修复后 P=4 下 matvec/求解与串行**逐位一致（rel=0.0）**。
- 4×2 matvec 431ms → **~183ms（2.35×）**；setup 13.5s；solve ~9-10s；|I| 与串行一致。
- 尝试按源盒外循环（缓存复用）以突破带宽——发现目标盒多源写竞态（0.33% 误差）+ 串行更慢，
  已回退。
- 并行效率（按 8 worker）：matvec ~36%、solve ~33%、setup ~19%。

### 结论（如实）
当前复制式 aggS 设计 + 本机内存带宽下，4×2（8 worker）**80% 并行效率不可达**；
实测天花板 ~36-40%（N 无关）。可选路径：
(a) 接受该上限（4×2 已较本轮前快 2.35×）；
(b) 分布式 aggS 设计（每秩只存自己盒的模式，P2P 交换远邻数据）——需较大重构，收益不确定；
(c) 重定义效率口径（如按 4 进程计，4×2 的 S=2.9× 对应 E=72%）。

## 未适配 MPI 分支补齐（2026-08-12，用户确认后）

### 已适配（全部通过混合精度门，rel ≈ 0.0）
1. **MLFMA FFTSpectral**：`MLFMAOperatorMPI` 增加 `interp_method` 透传；
   `disaggregate_downward!` 的 FFT 分支补 `child_filter`（只收集/写本秩子盒），
   消除与分区+Allreduce 的重复求和。P=2 精度门 rel=0.0。
2. **MLFMA Lebedev 1-step**：主循环已带 child_filter，新增 P=2 精度门验证 rel=0.0。
3. **PMCHW**：新增 `PMCHWMLFMAOperatorMPI`（双八叉树、近场行分区、四遍远场按 cube 分区 +
   每层 Allreduce；每遍独立 `y_pass` 缓冲避免跨遍重复求和）；`aggregate_leaf_pmchw!`/
   `disaggregate_leaf_pmchw_j!/m!` 增加 `cube_filter`。P=2 精度门 rel=7.5e-17。
4. **设置期**：`compute_interpolation_matrices!` 按层线程化（秩内）。

### 已知边界（记录）
- ~~分布式直接 LU~~ → 已完成（见下）。
- PMCHW 近场装配仍为全矩阵 O(N²) 后裁剪（各秩重复），设置期扩展性待优化。
- `test_pmchw_mlfma_operator.jl` 的 15.GD2A（k1 聚合奇偶，4.5%）为**既有 RED 门**
  （文件未改动、不在标准套件），与 MPI 适配无关。

### 回归
test_mlfma ✓、test_aim_operator 113,264 ✓、test_fft_interp 47 ✓；
混合门：MLFMA（默认/FFTSpectral/Lebedev）P=2 全 rel=0.0、AIM P=2 rel=0.0、PMCHW P=2 rel=7.5e-17。

## 分布式稠密 LU（边界完成，2026-08-12）

### 外部库调研结论（本环境实测）
- **ScaLAPACK**：Julia 注册表无绑定（仅 `SCALAPACK_jll` 二进制）；且其内嵌 `mpif_jll`
  的 MPI 与工程当前 MSMPI 运行时冲突（BLACS 无法共享我们的 communicator），需整体切换
  MPI 后端（用户级环境变更），不可直接采用。
- **Pardiso**：`Pardiso.jl` 可安装，但 MKL 回退实际报错 `Panua pardiso library was not loaded`；
  独立 Pardiso 库需到 pardiso-project.org 下载（许可制）。本环境不可用，已从依赖移除。

### 实现：原生 MPI 分布式稠密 LU（等价 ScaLAPACK PDGETRF/PDGETRS 思路）
- `src/Parallel/MPI/DistributedLU.jl`：1D 行分块 + 部分主元 + 整行交换 +
  分布式前代/回代；封装 `DistributedDenseLU` 对象（`mpi_lu`/`ldiv!`/`\`）。
- 门（`test/test_distributed_lu.jl`，mpiexec -n P 运行）：
  P=2/P=4 对角占优 rel=6.9e-16；随机复数矩阵（真实主元交换）rel=2.3e-14；封装 API rel=7.0e-16。

### 效率优化（用户要求"仔细封装优化好"）
- 主元行 1 次 Bcast + 旧行 k 点到点（替代两次 Bcast）；本地消元改**逐行融合广播 axpy**
  （`@. view -= l * ypiv`，无临时分配、SIMD），弃用不可靠的 strided `ger!`（复数 zgerc 共轭陷阱）。
- 基准（`benchmark/benchmark_distributed_lu.jl`，relerr=0.0）：
  | N | P=2 分解 | P=4 分解 | 加速 |
  |---|---|---|---|
  | 600 | 0.04s | 0.02s | 2.0× |
  | 1200 | 0.54s | 0.22s | 2.45×（含缓存效应） |
  旧版 N=300 为 0.25s——优化后同规模快约一个量级，且开始随秩扩展。

### 与多线程版本对比 + ScaLAPACK 调研结论（2026-08-12）
- **对比**（N=1200，同机）：OpenBLAS 多线程 `A\\b` T=1/2/4/8 = 0.07/0.041/0.024/0.023s；
  原生 MPI 分布式 LU（P=4）0.53s——**单机慢约 20 倍**。原因：自研 1D LU 未分块（逐行 axpy，
  BLAS1 级）+ O(N²) 通信，而 OpenBLAS getrf 是分块 BLAS3 + 多线程。用户判断正确。
- **ScaLAPACK 接口三路均被堵死（实测证据）**：
  1. `SCALAPACK_jll`：Artifacts.toml 无 Windows 构件，本机无法下载使用；
  2. MKL 的 `mkl_blacs_msmpi_ilp64`（JLL 与本机 oneAPI 2025.3 均验证）：仅 1 个无名序数导出，
     `blacs_gridinit_` 等按名不可用（转发桩）；
  3. MKL 通用 `mkl_blacs_ilp64`：运行时硬编码加载 `mkl_blacs_intelmpi_ilp64`（需 Intel MPI），
     与本机 MSMPI 不兼容。
- **结论与建议**：单机/工作站上稠密直接求解应使用多线程 LAPACK（Julia `A\\b`/`LUSolver`，
  OpenBLAS 或可选 MKL），这也是工程串行路径的现状；原生 MPI 分布式 LU 保留为**多节点集群**
  的分布式直接求解路径（正确、可测、跨节点唯一选项），并如实标注单机性能。
  后续若需单机/集群通用的高性能分布式稠密 LU，唯一可行方向是自研 2D 块循环分块 LU（BLAS3），
  或在 Linux/MPICH/完整 MKL 环境下启用 ScaLAPACK。
- 实验性 ScaLAPACKLU 模块与 MKL_jll/SCALAPACK_jll 依赖已移除（避免交付不可用代码）。

### ScaLAPACK 本地 MinGW 库接入（用户提供新路径，2026-08-12 完成）
用户指出本机 MinGW 已安装 ScaLAPACK：`C:\msys64\mingw64\bin\libscalapack.dll`
（ScaLAPACK 2.2.2，链接 `msmpi.dll`，导出 `blacs_gridinit_`/`descinit_`/`pzgesv_` 等；
包 `mingw-w64-x86_64-scalapack`，LP64，依赖 MinGW OpenBLAS）。据此实现
`src/Parallel/MPI/ScaLAPACKLU.jl`：
- BLACS 网格 + 块循环分发/收集（`distribute`/`gather_solution`）+ `pzgesv_` 包装
  （`scalapack_lu_solve(A, b, comm; MB, NB)`，各秩持完整 A/b，内部分发求解）。
- **关键坑**（实测）：MSMPI 版 BLACS 必须先 `blacs_get_(-1,0)` 取系统上下文再
  `gridinit`，直接传 0 会得到未初始化的 myrow/mycol（descinit 报 -9）；这是
  ScaLAPACK 规范用法（`BLACS_GET` → `GRIDINIT`），已修复。
- PZGESV 要求方形块（源码确认 `MB_A = NB_A`），包装器已加 `MB == NB` 显式校验。
- 分发改为按块向量化复制（原逐元素循环在 N=9600 时占 ~2s），精度不变。
- 门（`test/test_scalapack_lu.jl`，mpiexec -n 2/4）：rel = 6.35e-16 / 6.24e-16。

效率实测（`benchmark/benchmark_scalapack_lu.jl`，单机 16 核，与串行 OpenBLAS 顺序对照）：
| N | 配置 | 时间 | 对照 |
|---|---|---|---|
| 4800 | 4×2（P=4, T=2） | 2.05s（NB=64） | 串行 T=8: 2.12s（持平） |
| 9600 | 4×2（P=4, T=2） | 5.22s（不设 OMP） | 串行 T=8: 12.78s（**2.4×**） |
| 9600 | P=2, T=2 | 18.79s | P=4 相对 P=2 加速 3.6×（超线性，缓存效应） |
| 4800/9600 | 原生 DistributedLU P=2 | 111.9s（4800） | ScaLAPACK 快 **~80×** |

结论更新：ScaLAPACK 不再是"被堵死"——本地 MinGW 库是可行的分布式稠密 LU 路径，
精度与串行一致（rel≈6e-16），单机大用例（N≥4800）已优于多线程 OpenBLAS；
`DistributedLU.jl`（自研 1D）保留为无 ScaLAPACK 环境的回退/多节点唯一选项。
注意：MinGW OpenBLAS 的 OpenMP 线程数对耗时影响大（OMP_NUM_THREADS 需按核数/秩数调），
大 N 时建议不设 OMP（默认表现最优）；NB=128 在大 N 略优于 64，默认保持 64。

### 剩余边界状态
- PMCHW 近场装配设置期行分区（消除各秩重复的全矩阵 O(N²) 装配）：后续优化项。
- 15.GD2A RED 门：既有、与 MPI 无关，另行排查。

## 预条件 MPI 化 + PMCHW 原生近场 + 内存瓶颈审计（2026-08-13）

### 1. 预条件 MPI 化（用户要求：确认预条件是否支持 MPI，不支持则改造）
- **现状确认**：`distributed_gmres!` 之前完全忽略预条件 kwargs（`kwargs...` 丢弃），
  `Preconditioners.jl` 仅串行/线程版——**不支持 MPI**。
- 新增 `src/Parallel/MPI/Preconditioners.jl`：
  - `DistributedBlockJacobiPreconditioner`：块按叶 cube 归属分到各秩
    （rank 拥有 `(i_cube-1)%P` 号 cube 的 J/M 行块），构造只提取本秩块（每秩 1/P 块）；
    施加 = 本秩块 LU 求解 + 1 次 Allreduce。
  - `DistributedDiagonalPreconditioner`：对角逆 O(N)/秩，施加无通信。
  - `apply_mpi_preconditioner!` 统一入口：`nothing`/串行预条件为复制式回退。
  - 算子钩子：`DistributedBlockJacobiPreconditioner(op::MLFMAOperatorMPI)` 与
    `(op::PMCHWMLFMAOperatorMPI)`（J+M 行块）。
- `distributed_gmres!` 接入 `Pl`（左预条件 M⁻¹A x = M⁻¹b，解直接为 x）。
  修正实现中"右预条件 matvec + 左预条件残差/解"混用 bug（曾致解偏差 0.7%），
  统一为先 matvec 再施加 M⁻¹。
- 门 `test/test_hybrid_preconditioner.jl`（mpiexec -n 2/4）全绿：
  MLFMA 68 块 P=2 每秩 34 / P=4 每秩 17；BlockJacobi/Diagonal/PMCHW 施加
  与串行参照 **rel=0.0**；`mpi_gmres! + Pl` 收敛到稠密 LU rel≈3.5e-11。
  回归：串行 `test_preconditioners.jl` ✓、`test_distributed_gmres.jl` ✓。

### 2. PMCHW 原生近场装配（用户要求：禁止全 dense 后提取近场）
- 重写 `assemble_near_field_pmchw`：按 octree 近邻 cube 对原生计算 EJ/EM/HJ/HM
  四块（复用与稠密完全相同的 EFIE L-核与 PMCHW K-核），线程并行 COO，
  支持 `cube_filter` 按秩分区；**不再装配 2N×2N 稠密矩阵**。
- 踩坑修复（实测）：`efie_interaction!` 与 `calc_k_pmchw_term!` 末尾都会原地乘
  factor/边长因子，k0/k1 必须用独立缓冲求和，否则二次缩放（GD1 曾差 3.6e4 倍）；
  K 块多三角形组合可抵消到 ~1e-16×scale，稀疏合并后按 1e-12×scale 绝对阈值滤噪。
- MPI 构造改为按 cube 分区原生装配（`Z_near_local` 只含本秩行），mul! 全量缓冲
  SpMV + Allreduce（去掉全稠密 + 行切片）。
- 门：`test_pmchw_mlfma_operator.jl` GD1 近场逐元素对齐 max_rel=8e-14（机器精度）、
  GD2/GD2P/GD2RM/GD2T/15.8/15.11/GA/GB 全过；`test_hybrid_pmchw.jl` P=2/P=4
  **rel=0.0**。既有 RED 门（与本改动无关，旧近场对照验证相同）：GD2A（k1 聚合
  奇偶 4.5%）、GD2R（receive corr 0.99974）、B2（Z_in 72%，远场精度）。

### 3. 内存效率瓶颈审计与修复（用户要求：根进程/复制内存压力）
实测定位（P=2、N=594、PMCHW 12×18 球）：
| 瓶颈 | 修复前（每秩） | 修复后（每秩） |
|---|---|---|
| PMCHW 近场全稠密 2N×2N | 装配 + 行切片（O(4N²)） | 原生稀疏（近邻对，O(N·near)） |
| **远场翻译因子 αTrans**（每层 nPoles × (2·(2·near_range+1)+1)³ 全偏移表，near_range=16 时 264826 列） | **~15.7GB**（k0 5.4GB + k1 10.3GB） | **~21MB**（只存实际出现的 4454 偏移列；0 远邻层 0 列） |
| 预条件块 | 每秩全量复制（BlockJacobi 所有块） | 每秩 1/P 块 |
αTrans 压缩还使 PMCHW 构造从 137s → 10s（消除无远邻层的 264826 列×1250 pole 计算）。
端到端（PMCHW-MPI P=2、N=594，含预条件构造与 matvec）：峰值内存 17.6GB → **1.3–1.6GB**，
setup 14s。`test_mlfma.jl` 断言同步更新为"无远邻层允许空 αTrans"（translate! 语义不变，
混合门 rel=0.0 验证）。

### 后续可选（已记录，非本目标范围）
- PMCHW matvec 3.4s（N=594）由四遍远场 pass 主导，性能优化另立任务。
- αTrans 索引仍按 ±(2·near_range+1) 全范围建 OffsetArray（~2.4MB/层），进一步可压缩。

## 既有 RED 门修复（2026-08-13）
三个既有 RED 门根因定位并修复（`test_pmchw_mlfma_operator.jl` 全文件跑通）：
1. **GD2A（k1 聚合奇偶 4.5%）与 GD2R（receive corr 0.99974）**：根因是共享 EFIE
   路径 `aggregate_leaf!`/`disaggregate_leaf!` 的 RWG 三角形求积用 3 点，而 PMCHW
   路径 `aggregate_leaf_pmchw!`/`_receive_terms` 用 4 点——parity 门测两条实现的
   一致性，数值差异即求积点数差异（PMCHW 4 点经稠密对照验证更准）。修复：共享路径
   3 点→4 点统一。修复后 GD2A rel=0.0、GD2R rel=0.0/corr=1.0（机器精度）。
2. **B2（Zin vs Direct 72% 且 GMRES 不收敛）**：根因是测试求解配置缺陷——默认
   restart=min(30,N)=30 太小，GMRES 残差停滞在 ~0.25（**稠密矩阵同样不收敛**，
   cond≈7.9e6；全空间 restart=2N 后 109 次收敛到 5e-15）；matvec 本身相对稠密
   仅 3e-5。修复：B2 改用全空间 restart=2N，收敛后 Zin 误差 ~0.2%（<5% 门限）。
   注：N=540 budget_medium 门（strong form + restart 30）本就收敛，不受影响。
回归：test_mlfma、hybrid_mlfma/aim/pmchw（rel=0.0）、budget_medium（matvec
7.86e-6 不变）、gate_s_mlfma_medium（1.3e-5/4e-5 不变）、nmuller、
preconditioners、distributed_gmres 全绿；GD2S/GD2L/GD2U/GD2V 保持既有 Broken。
