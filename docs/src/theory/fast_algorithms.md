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

更贴近 EMSuite 当前实现的写法，可以把远场 MVP 链路记成
$$
\mathbf{y}_{far} = \mathcal{D}_{leaf}
\left(
\prod_{\ell=2}^{L-1} \mathcal{D}_{\ell \to \ell+1}
\right)
\left(
\sum_{\ell=2}^{L} \mathcal{T}_{\ell}
\right)
\left(
\prod_{\ell=L-1}^{2} \mathcal{A}_{\ell+1 \to \ell}
\right)
\mathcal{A}_{leaf}(\mathbf{x}).
$$
其中：

*   $\mathcal{A}_{leaf}$：叶层把物理基函数电流投影成辐射方向图；
*   $\mathcal{A}_{\ell+1 \to \ell}$：child-to-parent upward pass；
*   $\mathcal{T}_{\ell}$：同层 far-neighbor translation；
*   $\mathcal{D}_{\ell \to \ell+1}$：parent-to-child downward pass；
*   $\mathcal{D}_{leaf}$：叶层把接收场再测试回基函数空间。

### 1.3 叶层聚合与测试的离散形式

对 RWG 基函数，EMSuite 当前叶层聚合使用的是“对每个支撑三角形做固定阶高斯积分，再在球面极化基上投影”的形式。若继续使用上一章的统一符号记号，则叶层辐射方向图可写成
$$
S_{n,p}(\hat{k})=
\sum_{i=1}^2
\int_{T_{n,i}}
\hat{e}_p(\hat{k}) \cdot
\left[
\frac{\tilde l_{n,i}}{2}
\left( \mathbf{r} - \mathbf{v}_{n,i}^{\mathrm{opp}} \right)
e^{jk\hat{k}\cdot(\mathbf{r} - \mathbf{r}_c)}
\right]
\, dS,
\qquad p \in \{\theta, \phi\}.
$$

叶层测试则是与之对偶的回投影：把叶层接收到的 $\theta/\phi$ 极化场在每个 RWG 支撑三角形上做积分，累加到对应基函数系数。当前实现中这一步与叶层聚合共享同一组支撑三角形几何与高斯点结构。

### 1.4 upward / downward pass 的实现约定

在当前仓库中：

*   upward pass 先在叶层得到 `aggS`，然后逐层做 child-to-parent 插值与相移；
*   translation 在每层对 far-neighbor 盒子使用预计算的 $\alpha_{trans}(\hat{k}, \mathbf{D})$ 逐极化相乘；
*   downward pass 先做相移，再做 anterpolation，把父层接收场累加回子层 `disaggG`。

理论上可记为
$$
\mathbf{S}_{\ell}^{parent}(\hat{k}) = \sum_{c \in \mathrm{kids}(parent)}
e^{jk\hat{k}\cdot(\mathbf{r}_c - \mathbf{r}_{parent})}
\, \mathcal{I}_{c \to parent}
\mathbf{S}_{\ell+1}^{c}(\hat{k}),
$$
$$
\mathbf{G}_{\ell}^{obs}(\hat{k}) = \sum_{src \in \mathrm{far}(obs)}
\alpha_{trans}(\hat{k}, \mathbf{D}_{obs,src})
\, \mathbf{S}_{\ell}^{src}(\hat{k}),
$$
$$
\mathbf{G}_{\ell+1}^{child}(\hat{k}) = \mathcal{I}^{-1}_{parent \to child}
\left[
e^{jk\hat{k}\cdot(\mathbf{r}_{child} - \mathbf{r}_{parent})}
\mathbf{G}_{\ell}^{parent}(\hat{k})
\right].
$$

这里 $\mathcal{I}$ / $\mathcal{I}^{-1}$ 分别表示 interpolation / anterpolation。`nLevels >= 3` 的当前开放问题，正是集中在这条插值与相移链是否与 Legacy 完全对齐。

### 1.5 一个容易误用的实现细节：向量顺序

MLFMA 八叉树内部会维护 `sorted_ids`，以便按盒子内连续区间遍历基函数；但这只是内部数据局部性优化。EMSuite 当前 `MLFMAOperator` 的外部接口约定仍是：

*   输入向量 $\mathbf{x}$ 采用物理 basis 的原始顺序；
*   输出向量 $\mathbf{y}$ 也回到物理 basis 的原始顺序；
*   `sorted_ids` 只在 leaf aggregation / leaf disaggregation 内部把“盒子区间索引”映射回原始 basis id。

因此，benchmark 或回归测试不能在调用 `mul!` 前自行把 RHS 改成 `sorted_ids` 顺序，也不能在结果返回后再做一轮反排。此前某些 MLFMA 偏差曾被这类顺序误用放大为伪失败。

### 1.6 误差控制
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
