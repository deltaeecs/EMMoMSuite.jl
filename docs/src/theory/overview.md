# Theory Overview

本部分从理论层面总览 EMMoMSuite 使用的建模与求解技术栈。全部公式均以
《2023 年北京大学博士学位论文》（贺晓阳，下称"论文"）的 LaTeX 源码为准，
与 `docs/src/theory/` 下各篇文档一一对应。

## 1. Problem Classes

- PEC 面散射：S-EFIE、S-MFIE、CFIE（论文式 (2-24)~(2-26)）。
- 均匀/非均匀介质：V-EFIE、V-MFIE、PMCHWT、N-Muller、SCFIE。
- 金属-介质混合（VSIE）：VS-EFIE、VS-MFIE（论文式 (2-29)~(2-32)）。
- 端口激励下的辐射/散射工作流：平面波、delta-gap、波导端口。

## 2. Discretization

- 面基函数：RWG；体基函数：SWG、RBF（屋顶）、PWC（三种朝向）、PWCHex。
- 统一矢量基函数记号（论文式 (2-34)~(2-36)）：

  $$
  f_n(r) = \frac{a_n^{\pm}}{C_\Omega J^{\pm}} \rho, \qquad
  \rho = r - r_0^{\pm}, \qquad \nabla \cdot f_n = \frac{a_n^{\pm}}{J^{\pm}}
  $$

  其中 $C_\Omega$ 对三角形取 2、四面体取 3、其它网格取 1。
- 矩阵装配采用伽辽金测试（$t_m = f_m$），配合几何感知的三角形/四面体/六面体求积。

## 3. Operator Layer

- 场-源算子（论文式 (2-21)~(2-22)）：

  $$
  E(r) = \eta \mathcal{L}[J] + \mathcal{K}[J^m], \qquad
  H(r) = \frac{1}{\eta} \mathcal{L}[J^m] - \mathcal{K}[J]
  $$

- 稠密算子提供参考实现与基准真值；PMCHWT 组织为 $2N \times 2N$ 分块系统，
  显式区分 EJ/EM/HJ/HM 语义。

## 4. Fast Algorithms

- MLFMA 通过八叉树分层与球面波加法定理加速远场相互作用；
  近场以稀疏矩阵装配并与矩阵自由远场算子合并。
- 球面求积支持球面高斯求积（$\theta$ 方向 Gauss-Legendre + $\phi$ 方向均匀）与
  **Lebedev 求积**（效率接近 1，点数约为球面高斯的 $2/3$）。
- Lebedev 层间采用**球面非规则矢量插值**：插值矩阵按行伪逆求解或神经网络训练，
  反插值使用插值矩阵转置（论文第 4 章）。
- 截断项经验公式（论文式 (2-42)）：

  $$
  \tau(l) \approx 1.73\, k a_l + 2.16\, d_0^{2/3} (k a_l)^{1/3}, \qquad d_0 = 3
  $$

## 5. Krylov and Preconditioning

- 直接法：LU；迭代法：GMRES、BiCGSTAB 等，配合 MLFMA 矩阵向量乘。
- 预条件：块 Jacobi（叶层盒分块）、ILU、稀疏近似逆（SAI）。
- SAI 论文推荐用 LU 分解而非 QR 分解（论文式 (2-48)~(2-49)）：

  $$
  (P_l^i)^H = (LU)^{-1} Z_{near}^{(C_{in}, C_i)}
  $$

## 6. Validation Philosophy

- Legacy 实现是配平检查的基准（源码级 parity）。
- 稠密-快速探针与中等规模回归门用于分离求解器误差与算子保真度。
- 介质传输问题显式追踪公式级对比（PMCHWT vs N-Muller）。
- 奇异/近奇异积分以论文的"奇异值提取 + 递推公式"整体方案为准，
  不引入经验常数校准。
