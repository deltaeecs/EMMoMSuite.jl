# 快速算法 (Fast Algorithms)

为了突破传统矩量法 $O(N^2)$ 的内存和计算瓶颈，必须引入快速算法。

## 1. 多层快速多极子算法 (MLFMA)

MLFMA 适用于电大尺寸问题，将复杂度降低至 $O(N \log N)$。

> 实现状态说明：本页描述的是 MLFMA 的理论流程。当前仓库中的通用 surface MLFMA 实现，已确认 `nLevels >= 3` 时的主要偏差集中在 upward/downward pass 的插值与相移链路，属于仍在排查的实现问题，而不是文档未刷新导致的假象。

### 1.1 加法定理
利用格林函数的加法定理，将源点和场点分离：
$$
\frac{e^{-jk|\mathbf{r}_{obs} - \mathbf{r}_{src}|}}{|\mathbf{r}_{obs} - \mathbf{r}_{src}|} \approx \int_{S^2} e^{-j\mathbf{k} \cdot (\mathbf{r}_{obs} - \mathbf{r}_{c,obs})} T_L(\mathbf{k}, \mathbf{D}) e^{j\mathbf{k} \cdot (\mathbf{r}_{src} - \mathbf{r}_{c,src})} d^2\hat{k}
$$

### 1.2 算法流程
1.  **八叉树分组**：建立层级结构。
2.  **聚合 (Aggregation)**：计算最底层的辐射方向图，并向上传递。
3.  **转移 (Translation)**：在同层间转移辐射方向图。
4.  **解聚合 (Disaggregation)**：向下传递接收场。

### 1.3 误差控制
截断项数 $L \approx kd + 1.8 (d_0)^{2/3} (kd)^{1/3}$，其中 $d$ 为盒子尺寸。

## 2. 低频 MLFMA (Low-Frequency MLFMA)

标准 MLFMA 在低频（盒子尺寸远小于波长）时会出现数值崩溃（Low-frequency breakdown），因为汉克尔函数发散。
**解决方案**：采用基于多极子展开（Multipole Expansion）而非平面波展开的形式，或者使用归一化的平面波展开 (N-MLFMA)。

## 3. 自适应交叉近似 (ACA)

ACA (Adaptive Cross Approximation) 是一种纯代数压缩方法，不依赖于格林函数的解析形式。

### 3.1 原理
对于远场相互作用矩阵块 $\mathbf{Z}_{block}$，其数值秩远小于维数。
ACA 通过自适应地选取行和列，将矩阵分解为低秩形式：
$$
\mathbf{Z}_{block} \approx \mathbf{U} \mathbf{V}^T
$$
其中 $\mathbf{U} \in \mathbb{C}^{M \times r}, \mathbf{V} \in \mathbb{C}^{N \times r}$，且 $r \ll M, N$。

### 3.2 优势
*   **核无关 (Kernel-Independent)**：适用于多层介质格林函数或其他复杂核。
*   **易于实现**：无需复杂的解析展开。
