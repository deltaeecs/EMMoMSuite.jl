# 矩量法实现 (Method of Moments Implementation)

本章详细介绍矩量法（MoM）的数值实现细节，包括矩阵填充、数值积分和奇异性处理。

## 1. 伽辽金法 (Galerkin Method)

将积分方程写为算子形式 ``L(\mathbf{J}) = \mathbf{E}^{inc}``。
将电流展开为 ``\mathbf{J} = \sum I_n \mathbf{f}_n``。
选用测试函数 ``\mathbf{t}_m = \mathbf{f}_m`` 对方程两边做内积：
```math
\sum_{n=1}^N I_n \langle \mathbf{f}_m, L(\mathbf{f}_n) \rangle = \langle \mathbf{f}_m, \mathbf{E}^{inc} \rangle
```
得到线性方程组 ``\mathbf{Z} \mathbf{I} = \mathbf{V}``。

## 2. 矩阵元素计算

### 2.1 EFIE 阻抗矩阵
```math
Z_{mn}^{EFIE} = j\omega\mu \int_{S_m} \int_{S_n} \mathbf{f}_m(\mathbf{r}) \cdot \mathbf{f}_n(\mathbf{r}') G(R) dS' dS + \frac{1}{j\omega\epsilon} \int_{S_m} \int_{S_n} (\nabla \cdot \mathbf{f}_m) (\nabla' \cdot \mathbf{f}_n) G(R) dS' dS
```
计算涉及 4 个三角形对 ``(T_m^\pm, T_n^\pm)`` 的积分累加。

若按上一章的统一记号，把 RWG 的方向折叠到带符号边长 ``\tilde l_{m,i}``、``\tilde l_{n,j}`` 中，则矩阵元更贴近实现地写成
```math
Z_{mn}^{EFIE} = \sum_{i=1}^2 \sum_{j=1}^2 \left( Z_{mn,L}^{(i,j)} + Z_{mn,\phi}^{(i,j)} \right),
```
其中
```math
Z_{mn,L}^{(i,j)} = j\omega\mu \frac{\tilde l_{m,i} \tilde l_{n,j}}{4 A_{m,i} A_{n,j}}
\int_{T_{m,i}} \int_{T_{n,j}}
\left( \mathbf{r} - \mathbf{v}_{m,i}^{\mathrm{opp}} \right)
\cdot
\left( \mathbf{r}' - \mathbf{v}_{n,j}^{\mathrm{opp}} \right)
G(R)
\, dS' dS,
```
```math
Z_{mn,\phi}^{(i,j)} = \frac{1}{j\omega\epsilon}
\frac{\tilde l_{m,i}}{A_{m,i}}
\frac{\tilde l_{n,j}}{A_{n,j}}
\int_{T_{m,i}} \int_{T_{n,j}} G(R)\, dS' dS.
```

这样四个子三角形配对使用同一套几何模板，差别只体现在 ``A``、``\mathbf{v}^{\mathrm{opp}}`` 和带符号边长 ``\tilde l`` 上。EMSuite 当前 RWG 路径正是按这种“局部几何 + 局部符号”分解来实现装配的。

### 2.2 MFIE 阻抗矩阵
```math
Z_{mn}^{MFIE} = \frac{1}{2} \int_{S_m} \mathbf{f}_m \cdot \mathbf{f}_n dS - \int_{S_m} \mathbf{f}_m \cdot \left( \hat{n} \times \int_{S_n} \mathbf{f}_n \times \nabla' G dS' \right) dS
```

在 RWG-RWG 离散下，MFIE 的奇异主值项与 EFIE 一样会拆成四个支撑子三角形配对。若继续采用带符号边长记号，则局部法向、叉乘方向和 RWG 正负号可以统一折叠到局部几何因子中，而不需要在公式层面单独维护“正半/负半”两套写法。

### 2.3 体积分方程 (VIE) 阻抗矩阵

对于非磁性介质（``\mu=\mu_0``），VIE 方程为：
```math
\mathbf{D}(\mathbf{r}) - \epsilon_0 (\epsilon_r(\mathbf{r}) - 1) \int_V \overline{\mathbf{G}}_e(\mathbf{r}, \mathbf{r}') \cdot \frac{\mathbf{D}(\mathbf{r}')}{\epsilon_0 \epsilon_r(\mathbf{r}')} dV' = \epsilon_0 \mathbf{E}^{inc}(\mathbf{r})
```
使用 SWG 基函数 ``\mathbf{f}_n`` 展开电通量密度 ``\mathbf{D} = \sum D_n \mathbf{f}_n``。
阻抗矩阵元素 ``Z_{mn}`` 涉及两个四面体 ``V_m, V_n`` 的体积分：
```math
Z_{mn}^{VIE} = \int_{V_m} \frac{\mathbf{f}_m \cdot \mathbf{f}_n}{\epsilon_0 \epsilon_r} dV - \int_{V_m} \mathbf{f}_m \cdot \left( k_0^2 \int_{V_n} (\epsilon_r^{-1}-1) \mathbf{f}_n G dV' + \nabla \int_{V_n} (\epsilon_r^{-1}-1) (\nabla' \cdot \mathbf{f}_n) G dV' \right) dV
```
注意第二项中的 ``\nabla\nabla`` 操作通常通过分部积分转移到测试函数上，类似于 EFIE。

### 2.4 面-体耦合积分方程 (Surface-Volume Coupled IE)

对于涂覆介质的金属目标，需联立求解表面电流 ``\mathbf{J}_S`` 和体极化电流 ``\mathbf{J}_V``。
方程组结构：
```math
\begin{bmatrix}
Z_{SS} & Z_{SV} \\
Z_{VS} & Z_{VV}
\end{bmatrix}
\begin{bmatrix}
I_S \\
I_V
\end{bmatrix} = \begin{bmatrix}
V_S \\
V_V
\end{bmatrix}
```
*   ``Z_{SS}``: 传统的 EFIE 矩阵（RWG-RWG）。
*   ``Z_{VV}``: 传统的 VIE 矩阵（SWG-SWG）。
*   ``Z_{SV}``: 体电流产生的场在表面上的测试（RWG-SWG）。
*   ``Z_{VS}``: 表面电流产生的场在体内的测试（SWG-RWG）。

计算 ``Z_{SV}`` 时，源点在四面体 ``V_n`` 内，场点在三角形 ``S_m`` 上。需处理混合维度的积分。

## 3. 数值积分 (Numerical Integration)

### 3.1 高斯求积 (Gaussian Quadrature)
对于非奇异积分，采用三角形高斯求积公式。
```math
\int_T g(\mathbf{r}) dS \approx A_T \sum_{i=1}^{N_{quad}} w_i g(\mathbf{r}_i)
```
通常外层积分采用低阶（如 3 点），内层积分采用高阶（如 7 点或更高）。

对于 EMSuite 当前 RWG surface 路径，理论文档中常见的双面积分通常会在实现里进一步拆成：

*   外层测试三角形上的固定阶高斯点；
*   内层源三角形上的固定阶高斯点；
*   对 self / near / far 三类配对分别走不同的核函数或奇异性处理路径。

因此，推导时最好先写出“支撑配对求和”，再写每个配对子项的数值积分结构，这样更容易与代码一一对应。

## 4. 奇异性处理 (Singularity Treatment)

当源三角形 ``T^{src}`` 和场三角形 ``T^{obs}`` 重合或共边时，格林函数 ``G \sim 1/R`` 出现奇异性。

### 4.1 奇异性减法 (Singularity Subtraction)
```math
\int_T F(\mathbf{r}') \frac{1}{R} dS' = \int_T \left( \frac{F(\mathbf{r}')}{R} - \frac{F(\mathbf{r})}{R} \right) dS' + F(\mathbf{r}) \int_T \frac{1}{R} dS'
```
第一项数值积分，第二项解析积分。

### 4.2 Wilton 解析积分
对于静态奇异项 ``\int_T \frac{1}{R} dS'``，Wilton (1984) 给出了基于几何参数的解析公式，涉及对数项和反正切项。

### 4.3 Duffy 变换 (Duffy Transformation)
另一种处理奇异性的方法是将三角形域变换为正方形域，通过变量代换消除雅可比行列式中的奇异性。
```math
\int_0^1 \int_0^{1-x} \frac{f(x,y)}{\sqrt{x^2+y^2}} dy dx \xrightarrow{x=u, y=u v} \int_0^1 \int_0^1 \frac{f(u, uv)}{u\sqrt{1+v^2}} u dv du
```
变换后被积函数平滑，可直接使用高斯积分。

### 4.4 近奇异性处理 (Near-Singularity)
当源与场非常接近但不重合时（如细导线、薄介质层），积分呈现剧烈变化。

#### 4.4.1 投影变换法 (Projection Method)
将源三角形 ``T'`` 投影到观测点 ``\mathbf{r}`` 在 ``T'`` 所在平面的投影点 ``\mathbf{r}_0`` 为中心的极坐标系 ``(\rho, \phi)`` 中。
```math
\frac{1}{R} = \frac{1}{\sqrt{\rho^2 + d^2}}
```
其中 ``d`` 是观测点到源平面的垂直距离。
积分转化为：
```math
\int_{T'} \frac{1}{R} dS' = \sum_{i=1}^3 \int_{\phi_i}^{\phi_{i+1}} \int_{0}^{R_i(\phi)} \frac{\rho}{\sqrt{\rho^2 + d^2}} d\rho d\phi
```
内层 ``\rho`` 积分有解析解 ``\sqrt{R_i(\phi)^2 + d^2} - d``。外层 ``\phi`` 积分通常采用数值积分。

#### 4.4.2 自适应细分 (Adaptive Subdivision)
根据距离 ``d`` 和源单元尺寸 ``h`` 的比值，递归地将源三角形划分为 4 个子三角形。
准则：如果 ``d < C \cdot h``，则细分。
对每个子三角形，如果满足准则则直接高斯积分，否则继续细分或使用投影法。
这种方法精度高但计算量大。

#### 4.4.3 变量代换法 (Sinh Transformation)
对于一维线积分的近奇异性，常用 Sinh 变换消除 ``1/\sqrt{x^2+d^2}`` 类型的奇异性。
```math
x = d \sinh(u) \implies dx = d \cosh(u) du
```
```math
\sqrt{x^2+d^2} = d \cosh(u)
```
被积函数中的奇异项被完全抵消。
