# 激励与端口 (Excitations and Ports)

本章介绍如何在矩量法中施加不同类型的激励源。

## 1. 平面波激励 (Plane Wave)

用于计算雷达散射截面 (RCS)。
$$
\mathbf{E}^{inc}(\mathbf{r}) = \mathbf{E}_0 e^{-j\mathbf{k} \cdot \mathbf{r}}
$$
激励向量元素：
$$
V_m = \int_{S_m} \mathbf{f}_m(\mathbf{r}) \cdot \mathbf{E}_0 e^{-j\mathbf{k} \cdot \mathbf{r}} dS
$$

## 2. Delta-Gap 电压源

用于线天线或简单缝隙天线的馈电。
假设端口位于第 $k$ 条边，施加电压 $V_g$。
$$
V_m = \begin{cases}
V_g, & m = k \\
0, & m \neq k
\end{cases}
$$
输入阻抗 $Z_{in} = V_g / I_k$。

## 3. 偶极子源 (Dipole Source)

用于模拟附近的点源激励。
电偶极子位于 $\mathbf{r}_0$，矩为 $\mathbf{p}$。其产生的场 $\mathbf{E}^{inc}$ 为解析公式（格林函数的导数）。
$$
\mathbf{E}^{inc}(\mathbf{r}) = \frac{1}{\epsilon} \nabla \times \nabla \times (\mathbf{p} G(\mathbf{r}, \mathbf{r}_0))
$$
直接将此场代入 $V_m$ 积分公式计算。

## 4. 波导端口 (Waveguide Port)

用于微带线、波导等传输线结构的馈电。
**原理**：
在端口面上，假设场分布为波导模式的组合（如 TEM, TE10）。
$$
\mathbf{E}_{port} = \sum_{i} a_i \mathbf{e}_i(\mathbf{r})
$$
通过模式匹配法将端口场与内部矩量法区域耦合。这通常需要引入额外的未知量（模式系数）。
