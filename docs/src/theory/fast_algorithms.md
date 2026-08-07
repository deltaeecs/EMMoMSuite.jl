# 快速算法

传统矩量法在矩阵显式装配和矩阵向量乘阶段通常需要 $O(N^2)$ 的时间与存储成本。对电大尺寸问题，必须引入快速算法来降低复杂度。

## 1. 多层快速多极子算法

MLFMA 通过八叉树分层、方向图聚合、远场转移和解聚过程，把远场相互作用的复杂度降低到接近 $O(N \log N)$。

### 1.1 基本思想

核心出发点是利用格林函数的加法定理，把源点与观测点依赖拆开，改写为“源盒辐射方向图 + 盒间转移 + 观测盒接收方向图”的组合。

一个典型的远场链路可抽象写成

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

- $\mathcal{A}_{leaf}$ 表示叶层聚合，把物理基函数系数投影成方向图。
- $\mathcal{A}_{\ell+1 \to \ell}$ 表示 upward pass 中的子盒到父盒传递。
- $\mathcal{T}_{\ell}$ 表示同层远邻盒间的转移。
- $\mathcal{D}_{\ell \to \ell+1}$ 表示 downward pass 中的父盒到子盒传递。
- $\mathcal{D}_{leaf}$ 表示叶层测试，把接收方向图投回物理基函数空间。

### 1.2 算法流程

1. 八叉树划分：按层级把几何对象组织进空间盒结构。
2. 聚合：在叶层计算方向图，并逐层向上传递。
3. 转移：在同层远邻盒之间应用平移算子。
4. 解聚：把远场贡献逐层向下传回叶盒。
5. 叶层测试：把接收场重新积分回离散基函数系数。

### 1.3 与 EMMoMSuite 实现相关的叶层形式

对 RWG 基函数，EMMoMSuite 当前的叶层聚合采用“支撑三角形高斯积分 + 球面极化投影”的实现路线。若沿用统一符号记号，则一个极化分量的叶层方向图可写为

$$
S_{n,p}(\hat{k}) =
\sum_{i=1}^{2}
\int_{T_{n,i}}
\hat{e}_p(\hat{k}) \cdot
\left[
\frac{\tilde l_{n,i}}{2}
\left(\mathbf{r} - \mathbf{v}_{n,i}^{\mathrm{opp}}\right)
e^{j k \hat{k} \cdot (\mathbf{r} - \mathbf{r}_c)}
\right] dS,
\qquad p \in \{\theta, \phi\}.
$$

叶层测试则是上述过程的对偶回投影，即把接收的极化场在同一组支撑三角形上重新积分并累加回离散自由度。

### 1.4 upward / downward pass

在实现层面，upward pass 和 downward pass 可以理解为插值、相移与反插值的串联。

$$
\mathbf{S}_{\ell}^{parent}(\hat{k}) = \sum_{c \in \mathrm{kids}(parent)}
e^{j k \hat{k} \cdot (\mathbf{r}_c - \mathbf{r}_{parent})}
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
e^{j k \hat{k} \cdot (\mathbf{r}_{child} - \mathbf{r}_{parent})}
\mathbf{G}_{\ell}^{parent}(\hat{k})
\right].
$$

当前仓库中，`nLevels >= 3` 时已知主要风险集中在 upward/downward pass 的插值与相移链路，而不是整体接口层面的命名或调用方式。

### 1.5 一个容易误用的接口细节

MLFMA 八叉树内部通常会维护按盒排序的 `sorted_ids` 以优化局部访问，但 EMMoMSuite 对外暴露的 `MLFMAOperator` 仍然使用物理基函数的原始自由度顺序作为输入和输出。

因此：

- 调用 `mul!` 或 Krylov 求解器前，不应自行把向量改写成 `sorted_ids` 顺序。
- 结果返回后，也不应再做一次“反排序”处理。

否则会把内部数据局部性实现误当成外部接口约定，导致伪回归差异。

### 1.6 误差控制

MLFMA 的精度通常受截断项数、插值阶数、盒尺寸和层级深度共同控制。工程上需要通过 benchmark 和 Legacy 对齐来校验这些参数组合，而不是依靠经验常数补偿偏差。

## 2. 低频 MLFMA

标准平面波形式的 MLFMA 在低频下会出现 low-frequency breakdown，因为展开形式随盒尺寸和波数缩小而变得病态。

常见解决路线包括：

- 使用基于多极子展开的低频稳定形式。
- 使用归一化平面波展开或其他低频重标定技术。
- 在极低频段切换到其他更适合的压缩策略。

## 3. 自适应交叉近似

ACA 是一种核无关的代数压缩方法，适用于远场块具有低数值秩的场景。

### 3.1 基本形式

对远场块矩阵 $\mathbf{Z}_{block}$，ACA 尝试构造

$$
\mathbf{Z}_{block} \approx \mathbf{U} \mathbf{V}^{T},
$$

其中 $\mathbf{U} \in \mathbb{C}^{M \times r}$，$\mathbf{V} \in \mathbb{C}^{N \times r}$，并且 $r \ll M, N$。

### 3.2 特点

- 不依赖核函数的解析加法定理。
- 对复杂介质核或非标准格林函数更灵活。
- 通常更容易接入已有稀疏或分块框架，但常数因子和稳定性需要单独评估。
