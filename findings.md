# Findings & Decisions

## Requirements
- 搜集 CEM 知识库中 MLFMA 的加速手段（示例：FFT 加速矩阵向量乘积）
- 为 EMMoMSuite 制定这些加速算法的实施方案
- 方案需可执行、可实测门控、可审计

## Research Findings

### CEM 知识库相关出处（均已核实文件存在）
1. `topics/fast_algorithms.md`：FMM/MLFMA/ACA/MLACA 总览；指出 MLFMA 通过插值/反插值与对角化转移实现 O(N log N) 矩阵向量乘。
2. `books_parsed/ergul_mlfma/ch08_*`（Ch2，已解析）：MLFMA 对 MoM 矩阵向量乘的加速；`ch10_*`（Ch4，已解析）：并行化（聚合/转移/配置、负载均衡）。
   - 注意：Ch3（插值/反插值、采样、转移算子插值、FFT 相关内容主章节）**未解析成功**（索引多处标注 "Ch3 未提取成功"）。
3. `books_parsed/golub_matrix_computations/ch08_1._Matrix_Multiplication.md` §1.4：FFT 作为结构化矩阵向量乘的 O(n log n) 方法；注释中给出 NUFFT（Dutt & Rokhlin 1993）、Greengard & Lee 2004 等参考文献。
4. `books_parsed/golub_matrix_computations/ch11_4._Special_Linear_Systems.md` §4.8：Circulant/Toeplitz 系统用 FFT 对角化并求解，O(n log n)。
5. `books_parsed/jin_fem_em_3rd/ch17_Chapter_10_Finite_Element–Boundary_Integral_Methods.md` §10.6.2：平面/圆等特殊边界上的边界积分是卷积形式，可用 FFT 直接求矩阵向量积而不显式存储矩阵；一般几何用"adaptive integral 与 fast multipole"。
6. `books_parsed/jin_fem_em_3rd/ch20_Chapter_13_Finite_Element_Analysis_of_Periodic_Structures.md`：周期结构用离散傅里叶变换/FFT 把分块循环矩阵对角化，K 个激励一次 FFT 处理。
7. `books_parsed/chew_waves_inhomogeneous/ch03_2.md` §2.7：Sommerfeld 积分在多个 ρ 处取值时用 FFT 求值 + 插值降低采样成本。
8. `books_parsed/gibson_mom/ch17_7._Bodies_of_Revolution.md`：BoR 用 FFT/模式分解；`ch16_6._Two-Dimensional_Problems.md`：2D 条带矩阵为对称 Toeplitz 结构。

### 归纳出的加速手段清单
| 加速手段 | 原理 | 知识库出处 | 对 EMMoMSuite 的相关性 |
|---|---|---|---|
| MLFMA 远场压缩（聚合-转移-配置，对角化转移） | 近远场拆分，远场 O(N log N) | fast_algorithms.md；ergul Ch2/Ch4 | 已实现（MLFMAOperator） |
| FFT 加速矩阵向量乘（结构化核） | 卷积/Toeplitz 核用 FFT 对角化，O(n log n) | Golub §1.4/§4.8；Jin Ch10 §10.6.2；Chew §2.7 | 待实施：2D/周期/近场块；3D 一般核不直接适用 |
| 层间插值/反插值加速 | 低层→高层采样率匹配；Lagrange/谱插值；反插值为伴随 | fast_algorithms.md；ergul Ch3（未解析）；工程 docs/theory/fast_algorithms.md | 已实现 Lagrange2Step 与 Lebedev 1-step；**FFT 谱插值是主要增量** |
| FFT/谱插值（φ 方向带限插值） | φ 方向周期带限，DFT/FFT 可精确重构/求值，O(K log K) | Golub §1.4；Hansen Ch4 §4.3.3（周期带限函数 DFT/FFT 精确）；Ergul Ch3 §3.3.1（两步法基准） | 待实施（P1，推荐首选） |
| 非均匀 FFT（NUFFT） | 不规则采样（如 Lebedev 点）快速插值 | Golub ch08 注释（Dutt & Rokhlin 1993；Greengard & Lee 2004） | 待评估（P4，可选） |
| AIM/IE-FFT（网格投影 + FFT 卷积） | RWG 投影到均匀网格，远场 FFT 卷积，近场修正 | Jin Ch10 §10.6.2（"adaptive integral method"） | 待评估（P3，可选独立算子） |
| 周期结构 FFT 对角化 | 分块循环矩阵 FFT 对角化 | Jin Ch20 §13.10 | 远期（工程暂无周期结构支持） |
| 并行化 MLFMA | 混合/分层并行、负载均衡 | ergul Ch4（已解析） | 已实现 MPI+Threads；配合 FFT 布局优化（P5） |
| ACA/MLACA（代数低秩压缩） | 远场块低秩近似 | fast_algorithms.md；gibson ch19/ch20 | 工程 FastAlgorithms 未含 ACA 模块（基准结果文件为历史实验，无源码） |

### Ergul Ch3 原文要点（2026-08-12 由原版 PDF 直接提取，`Books/[Ozgur_Ergul,_Levent_Gurel]_The_Multilevel_Fast_Mu(z-lib.org).pdf` 第 193–284 页）
知识库中 Ch3 未入库解析，但本次只读提取了 PDF 文本层，核实以下内容：
- §3.2.5 采样：采样率由截断数决定，`S_θ = τ(l)+1`、`S_φ = 2(τ(l)+1)`（与工程 `levelIntegralInfoCal` 一致）。
- §3.3.1 两步法（Two-Step Method）：一步 Lagrange 为 8p²(τ+1)² 次操作；
  两步法 θ→φ 分解后为 `Θ_two-step = 4p(τ(l)+1)(τ(l)+τ(l-1)+2)`，比值 < 1/p，永远更快；
  表 3.2（CFIE、RWG、λ/10 网格、132,003~8,447,808 未知量、p=3、6×6 模板）实测聚合阶段快约 45%、MVM 快 25%~30%。
  结论：**Ergul 书中插值加速用的是两步 Lagrange，未使用 FFT**。
- §3.3.2 虚拟扩展 θ-φ 空间：边界外插值点的处理（工程 `Interpolation.jl` 的 `pickθ/pickϕ` 即实现此逻辑）。
- §3.3.3 极点采样：在 θ=0,π 采样降低极区插值误差（工程 Lebedev 路径的极点处理与此对应）。
- §3.3.4 转移算子插值：翻译算子用插值在 O(N) 设置时间内生成，减少 setup 成本。

### Hansen Ch4 §4.3.3（`books_parsed/hansen_spherical_nf/ch05_PBEW026E_ch4.md`）
"the chi and phi integrals are both Fourier integrals of periodic, band-limited functions.
Such integrals may be evaluated exactly by the Discrete Fourier Transform (DFT) technique
for which an efficient algorithm (the Fast Fourier Transform, FFT) exists."
——这是知识库中"周期带限函数可用 FFT 精确求值"的直接论述，
构成 P1（MLFMA φ 方向 FFT 谱插值）的数学依据：φ 方向正是周期带限。

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| P1 先做 GL 网格 φ 方向 FFT 谱插值 | 改动局部（Interpolation.jl + 插值矩阵构造），与现有 2-step Lagrange 结构同构，便于等价性回归 |
| 新增 FFTW/AbstractFFTs 依赖 | Julia 标准 FFT 生态；MLFMA 现有代码无 FFT 依赖 |
| AIM/IE-FFT、NUFFT 独立成模块 | 避免污染既有 MLFMA 路径；各自配稠密参照门 |
| 每阶段先写失败测试再实现（TDD） | 工程 AGENTS.md 要求；可审计 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| Ergul Ch3 未入库解析 | 已从原版 PDF（Books/）只读提取文本核实：书中插值加速为两步 Lagrange 而非 FFT；FFT 出处改用 Golub/Hansen/Jin/Chew；KB 入库补解析列为可选后续（需用户同意写入 KB） |
| 工程无 FFTW 依赖、无 ACA/周期支持 | 方案按"可立即实施 / 条件性 / 远期"分层，避免过度承诺 |

## Resources
- 知识库：`D:\华为家庭存储\projects\KnowledgeBase\CEM`（topics/fast_algorithms.md、index/concepts.md、books_parsed/{ergul_mlfma,golub_matrix_computations,jin_fem_em_3rd,chew_waves_inhomogeneous,gibson_mom}）
- 工程：`C:\Users\12253\OneDrive\MoM\EMMoMSuite\src\FastAlgorithms\MLFMA\`（Interpolation.jl、Translation.jl、Aggregation.jl、Disaggregation.jl、MLFMAOperator.jl、OctreeBuilder.jl）、`src\FastAlgorithms\Lebedev\`（LVI.jl、SHInterp.jl）
- 文档：`docs\src\theory\fast_algorithms.md`（博士论文式 MLFMA 实现说明，含插值/反插值公式）
- 测试/基准：`test\test_mlfma.jl`、`test\test_lebedev_interp.jl`、`benchmark\benchmark_mlfma_vs_direct.jl`
