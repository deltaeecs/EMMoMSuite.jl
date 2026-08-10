# 激励与端口 (Excitations and Ports)

本章介绍矩量法中的激励向量 $V_m = \langle \bm{f}_m, \bm{g} \rangle$ 的构造，
其中 $\bm{g}$ 为入射场或端口等效源。

## 1. 平面波激励 (Plane Wave)

用于计算雷达散射截面 (RCS) 与平面波照射下的散射问题：

$$
\bm{E}^{inc}(\bm{r}) = \bm{E}_0\, e^{-{\rm j}\bm{k} \cdot \bm{r}}
$$

激励向量元素为测试基函数与入射场的积分：

$$
V_m = \int_{S_m} \bm{f}_m(\bm{r}) \cdot \bm{E}_0\, e^{-{\rm j}\bm{k} \cdot \bm{r}}\, dS
$$

对体问题（VEFIE）相应地使用体积分 $\int_{V_m} \bm{f}_m \cdot \bm{E}^{inc} dV$。

## 2. Delta-Gap 电压源

用于线天线或简单缝隙天线的馈电。假设端口位于第 $k$ 个未知量处，施加电压
$V_g$：

$$
V_m = \begin{cases} V_g, & m = k \\ 0, & m \neq k \end{cases}
$$

输入阻抗为 $Z_{in} = V_g / I_k$。

## 3. 偶极子源 (Dipole Source)

用于模拟附近的点源激励。电偶极子位于 $\bm{r}_0$、偶极矩为 $\bm{p}$，
其产生的入射场由格林函数解析表示：

$$
\bm{E}^{inc}(\bm{r}) = \frac{1}{\varepsilon} \nabla \times \nabla \times
\left(\bm{p}\, G(\bm{r}, \bm{r}_0)\right)
$$

直接代入激励向量积分即可。

## 4. 波导端口 (Waveguide Port)

用于微带线、波导等传输线结构的馈电。在端口面上，把场分布展开为波导模式的
组合（如 TEM、TE10）：

$$
\bm{E}_{port} = \sum_i a_i \bm{e}_i(\bm{r})
$$

通过模式匹配将端口场与矩量法区域耦合，需要引入额外的模式系数未知量。
EMMoMSuite 的端口实现（`Ports` 模块）将端口面离散为等效源并装配到激励向量
或扩展阻抗矩阵中。
