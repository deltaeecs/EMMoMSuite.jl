# 求解器 (Solvers)

建立矩阵方程 ``\mathbf{Z} \mathbf{I} = \mathbf{V}`` 后，需采用高效的数值线性代数方法求解。

## 1. 直接求解器 (Direct Solvers)

适用于中小规模问题（``N < 10,000``）或需要多次求解不同右端项的情况。

### 1.1 LU 分解
将矩阵分解为下三角矩阵 ``\mathbf{L}`` 和上三角矩阵 ``\mathbf{U}``：
```math
\mathbf{P} \mathbf{Z} = \mathbf{L} \mathbf{U}
```
其中 ``\mathbf{P}`` 为置换矩阵（用于选主元）。
*   **计算复杂度**：``O(N^3)``
*   **内存复杂度**：``O(N^2)``
*   **求解步骤**：
    1.  前代 (Forward Substitution): ``\mathbf{L} \mathbf{y} = \mathbf{P} \mathbf{V}``
    2.  回代 (Backward Substitution): ``\mathbf{U} \mathbf{I} = \mathbf{y}``

### 1.2 Cholesky 分解
仅适用于对称正定矩阵（MoM 矩阵通常不是，但在某些变分形式下可能是）。
```math
\mathbf{Z} = \mathbf{L} \mathbf{L}^T
```
计算量约为 LU 分解的一半。

## 2. 迭代求解器 (Iterative Solvers)

适用于大规模问题，通常配合 MLFMA 使用。迭代法通过不断修正近似解，使残差 ``\mathbf{r}_k = \mathbf{V} - \mathbf{Z} \mathbf{I}_k`` 减小。

对 EMSuite 当前实现，更贴近代码的写法是左预条件形式
```math
\mathbf{M}^{-1} \mathbf{Z} \mathbf{I} = \mathbf{M}^{-1} \mathbf{V},
```
其中 `GMRESSolver` 会把左预条件子作为 `Pl` 传给底层 `IterativeSolvers.gmres`。因此，迭代历史中直接被最小化的是预条件残差
```math
	ilde{\mathbf{r}}_k = \mathbf{M}^{-1}(\mathbf{V} - \mathbf{Z} \mathbf{I}_k),
```
而最终数值验收时，仍建议同时检查物理残差
```math
\mathbf{r}_k = \mathbf{V} - \mathbf{Z} \mathbf{I}_k.
```

### 2.1 Krylov 子空间方法
在 Krylov 子空间 ``\mathcal{K}_m(\mathbf{Z}, \mathbf{r}_0) = \text{span}\{\mathbf{r}_0, \mathbf{Z}\mathbf{r}_0, \dots, \mathbf{Z}^{m-1}\mathbf{r}_0\}`` 中寻找近似解。

#### 2.1.1 GMRES (Generalized Minimal Residual)
*   **原理**：寻找 ``\mathbf{I}_m \in \mathbf{I}_0 + \mathcal{K}_m``，使得 ``\lVert \mathbf{V} - \mathbf{Z} \mathbf{I}_m \rVert_2`` 最小。
*   **特点**：
    *   适用于任意非奇异非对称矩阵。
    *   收敛平稳，单调下降。
    *   **缺点**：随着迭代步数 ``m`` 增加，需要存储所有正交基向量，内存消耗线性增长，计算量二次增长。通常使用重启 GMRES(m)。

在当前求解器里，`restart`、`maxiter`、`tol` 都直接映射到底层 GMRES 的参数。对 MLFMA 算子而言，一个实用习惯是同时记录：

*   预条件残差是否达到 `tol`；
*   物理残差 ``\lVert \mathbf{V} - \mathbf{Z} \mathbf{I} \rVert_2 / \lVert \mathbf{V} \rVert_2`` 是否同样收敛；
*   结果向量是否保持物理 basis 顺序，而不是误用 MLFMA 内部的 `sorted_ids` 顺序。

#### 2.1.2 BiCGSTAB (Biconjugate Gradient Stabilized)
*   **原理**：结合了 BiCG 和 GMRES 的思想，通过短递归公式更新解。
*   **特点**：
    *   内存消耗小且固定（仅需存储几个辅助向量）。
    *   每步迭代需要两次矩阵向量乘法 (MVP)。
    *   收敛曲线可能出现剧烈震荡，但在很多电磁问题中收敛速度快于 GMRES。

#### 2.1.3 CGS (Conjugate Gradient Squared)
*   **特点**：收敛速度通常比 BiCG 快，但震荡更剧烈，甚至可能发散。

#### 2.1.4 QMR (Quasi-Minimal Residual)
*   **特点**：旨在平滑 BiCG 的收敛曲线，适用于病态矩阵。

## 3. 预处理技术 (Preconditioning)

为了加速迭代收敛，求解预处理后的方程：
```math
\mathbf{M}^{-1} \mathbf{Z} \mathbf{I} = \mathbf{M}^{-1} \mathbf{V}
```
其中 ``\mathbf{M}`` 是 ``\mathbf{Z}`` 的近似，且易于求逆。理想的预处理子应使 ``\mathbf{M}^{-1} \mathbf{Z}`` 的特征值聚类于 1 附近。

### 3.1 块雅可比预处理 (Block-Jacobi)
利用 MLFMA 的近场矩阵（稀疏块对角）作为预处理子。
```math
\mathbf{M} = \text{diag}(\mathbf{Z}_{near})
```
*   **实现**：对每个对角块进行 LU 分解。
*   **优点**：计算简单，易于并行，对大多数散射问题效果显著。

在 EMSuite 当前 MLFMA 路径中，这些块通常直接来自叶层盒子的 basis 区间，也就是 `get_leaf_intervals` 或按 `sorted_ids[cube.bfInterval]` 抽取出来的近场子块。换言之，块雅可比不是任意切块，而是与 MLFMA 叶层几何分组绑定的预条件子。

### 3.2 稀疏近似逆 (SAI / SPAI)

稀疏近似逆 (Sparse Approximate Inverse, SAI) 旨在构造一个稀疏矩阵 ``\mathbf{M}``，使其直接逼近 ``\mathbf{Z}^{-1}``。预处理操作仅需进行一次稀疏矩阵向量乘法 (SpMV)，非常适合并行计算。

#### 3.2.1 数学表述
SAI 通过最小化 Frobenius 范数来构造 ``\mathbf{M}``：
```math
\min_{\mathbf{M} \in \mathcal{S}} \lVert \mathbf{M}\mathbf{Z} - \mathbf{I} \rVert_F^2 = \min_{\mathbf{M} \in \mathcal{S}} \sum_{k=1}^N \lVert \mathbf{M} \mathbf{z}_k - \mathbf{e}_k \rVert_2^2
```
其中 ``\mathcal{S}`` 是预先设定的稀疏模式，``\mathbf{z}_k`` 是 ``\mathbf{Z}`` 的第 ``k`` 列（实际上通常使用 ``\mathbf{Z}\mathbf{M} \approx \mathbf{I}`` 的形式，即最小化 ``\lVert \mathbf{Z}\mathbf{M} - \mathbf{I} \rVert_F``，此时问题分解为对 ``\mathbf{M}`` 的每一列 ``\mathbf{m}_k`` 的独立最小二乘问题）：
```math
\min_{\mathbf{m}_k} \lVert \mathbf{Z} \mathbf{m}_k - \mathbf{e}_k \rVert_2^2
```

#### 3.2.2 求解过程
由于 ``\mathbf{m}_k`` 是稀疏的，上述 ``N`` 维最小二乘问题可以简化为极小规模的子问题。
对于第 ``k`` 列 ``\mathbf{m}_k``：
1.  **确定稀疏模式**：设 ``\mathcal{J}_k = \{j \mid m_{jk} \neq 0\}`` 为 ``\mathbf{m}_k`` 中非零元素的索引集合。
2.  **确定相关行**：设 ``\mathcal{I}_k = \{i \mid \exists j \in \mathcal{J}_k, Z_{ij} \neq 0\}`` 为 ``\mathbf{Z}`` 中与 ``\mathcal{J}_k`` 相关的非零行索引。
3.  **构建子矩阵**：提取小矩阵 ``\hat{\mathbf{Z}}_k = \mathbf{Z}(\mathcal{I}_k, \mathcal{J}_k)``。
4.  **求解 LS 问题**：
    ```math
    \min_{\hat{\mathbf{m}}_k} \lVert \hat{\mathbf{Z}}_k \hat{\mathbf{m}}_k - \hat{\mathbf{e}}_k \rVert_2
    ```
    通常使用 QR 分解求解此小规模问题。

#### 3.2.3 稀疏模式选择
*   **静态模式 (Static Pattern)**：直接采用 ``\mathbf{Z}``（或其稀疏部分 ``\mathbf{Z}_{near}``）的稀疏模式，或者其幂次 ``\mathbf{Z}^p`` 的模式。
*   **动态模式 (SPAI 算法)**：自适应地添加非零元素以降低残差。如果当前残差 ``\lVert \mathbf{Z} \mathbf{m}_k - \mathbf{e}_k \rVert_2`` 过大，则根据梯度方向寻找新的索引加入 ``\mathcal{J}_k``，直到满足容差。

#### 3.2.4 优势与局限
*   **优势**：
    *   **完全并行**：每一列的计算完全独立。
    *   **应用快速**：预处理步骤仅为 SpMV。
    *   **鲁棒性**：对于不定矩阵或非对角占优矩阵通常比 ILU 更稳定。
*   **局限**：构造过程（尤其是动态模式）计算量较大。

### 3.3 不完全 LU 分解 (ILU)
对稀疏矩阵（如近场矩阵）进行不完全的 LU 分解，忽略小于阈值的填充元素。
*   **特点**：比 Block-Jacobi 更强，但构建和求解过程难以并行化。

## 4. 并行计算 (Parallel Computing)

### 4.1 MPI 分布式并行
适用于多节点集群。
*   **矩阵分布**：将阻抗矩阵按行块或二维块循环分布在不同进程中。
*   **MLFMA 并行**：
    *   八叉树的上层节点（共享部分）复制到所有进程。
    *   下层节点（独立部分）分布存储。
    *   通信主要发生在聚合/解聚合阶段的边界交换，以及转移阶段的辐射方向图交换。

### 4.2 OpenMP 线程并行
适用于单节点多核 CPU。
*   **循环并行化**：利用 `#pragma omp parallel for` 加速矩阵填充循环和数值积分循环。
*   **任务并行化**：在 MLFMA 中，不同盒子的计算任务可以动态分配给不同线程。

