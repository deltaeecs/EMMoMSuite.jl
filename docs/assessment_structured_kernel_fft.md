# 结构化核 FFT 加速适用性评估（P2）

日期：2026-08-12

## 1. 背景与知识库出处

FFT 加速矩阵向量乘的本质：把**卷积/Toeplitz（循环）结构的核**对角化，
矩阵向量乘从 O(n²) 降到 O(n log n)。CEM 知识库中的相关出处：

- Golub & Van Loan《Matrix Computations》§4.8（`books_parsed/golub_matrix_computations/ch11_4._Special_Linear_Systems.md`）：
  circulant/Toeplitz 系统用 FFT 求解/矩阵向量乘，O(n log n)。
- Jin《FEM in Electromagnetics》Ch10 §10.6.2（`books_parsed/jin_fem_em_3rd/ch17_...md`）：
  平面/圆等特殊边界上的边界积分是**卷积形式**，可用 FFT 直接求矩阵向量积而不显式存储矩阵；
  一般几何则用 adaptive integral（AIM）与 fast multipole。
- Jin Ch13（`books_parsed/jin_fem_em_3rd/ch20_...md`）：周期结构分块循环矩阵 → FFT 对角化。
- Gibson《MoM in Electromagnetics》Ch6（`books_parsed/gibson_mom/ch16_6._Two-Dimensional_Problems.md`）：
  2D 条带 MoM 矩阵为对称 Toeplitz。
- Chew《Waves in Inhomogeneous Media》§2.7（`books_parsed/chew_waves_inhomogeneous/ch03_2.md`）：
  Sommerfeld 积分在多个 ρ 处取值用 FFT + 插值。

## 2. 适用条件

1. 核函数只依赖**位置差**（卷积核），采样在**均匀网格或周期域**上；
2. 或问题具有周期/旋转对称性（分块循环、模式分解）。

## 3. 对 EMMoMSuite 的逐项评估

| 候选场景 | 卷积结构 | 工程现状 | 结论 |
|---|---|---|---|
| 3D 一般曲面 MoM（EFIE/MFIE/CFIE，RWG） | 无全局 Toeplitz 结构（格林函数依赖 |r−r′|，曲面任意） | 已有 MLFMA | **不直接适用**；这正是 MLFMA 存在的意义 |
| 2D 条带 MoM（Toeplitz） | 有 | 工程无 2D MoM 模块（仅 3D 表面/体积积分） | 远期：若新增 2D 支持，按 Gibson Ch6 + Golub §4.8 实现 |
| 周期结构（单元胞/周期边界） | 分块循环 | 工程无周期边界支持 | 远期：按 Jin Ch13 FFT 对角化 |
| 旋转体 BoR（模式分解） | 方位角循环 | 工程无 BoR 支持 | 远期 |
| MLFMA φ 方向层间插值/反插值 | 周期带限卷积 | 已在 P1 实现 FFT 谱插值 | **已落地**（本项目内唯一直接适用的卷积 FFT 点） |
| 一般几何均匀网格投影（AIM/IE-FFT/P-FFT） | 投影后远场为卷积 | 未实现 | P3 独立算子路径（见下） |

## 4. 结论与建议

1. **当前架构下的直接适用面**：只有 MLFMA 的 φ 方向层间插值/反插值（周期带限），
   已在 P1 用 FFT 谱插值落地并通过等价性/精度门。
2. **2D/周期/BoR 的结构化核 FFT**：需要先建立相应求解模块，属工程外新功能，按远期处理；
   不应在现有 3D MLFMA 路径中硬编码或引入 fallback。
3. **一般几何的 FFT 加速**：可行路径是 AIM/IE-FFT（Jin Ch10 提及），
   作为与 MLFMA 并列的独立算子（P3），需要 RWG→均匀网格投影、近场修正与 FFT 卷积三部分。
4. **优先级建议**：
   - P1 收尾：在更大规模（N ≥ 1e4、更高截断阶）下验证 φ 步 FFT 的 O(M log M) 优势
     （当前小规模下稀疏 ϕCSC 更快，见 benchmark_fft_interp.jl）；
   - P3 或 P4：AIM/IE-FFT（新算子）或 NUFFT（Lebedev 插值加速），按资源取舍；
   - 2D/周期/BoR：按产品需求另行立项。

## 5. 证据
- 知识库出处：`D:\华为家庭存储\projects\KnowledgeBase\CEM\books_parsed\{golub_matrix_computations,jin_fem_em_3rd,gibson_mom,chew_waves_inhomogeneous}`
  （OneDrive 副本：`C:\Users\12253\OneDrive\KnowledgeBase\CEM\books_parsed\...`，已同步）。
- 工程现状：`src/` 仅有 3D 表面/体积积分（RWG/SWG）与 MLFMA；无 2D/周期/BoR 模块。
- P1 基准：`benchmark/benchmark_fft_interp.jl`。
