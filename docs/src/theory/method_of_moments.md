# 矩量法实现

本章聚焦矩量法在代码中的离散实现，包括伽辽金离散、矩阵元素计算、数值积分以及奇异和近奇异积分处理。

## 1. 伽辽金离散

设连续积分方程写成

$$
L(\mathbf{J}) = \mathbf{E}^{inc}.
$$

将未知电流展开为

$$
\mathbf{J} = \sum_{n=1}^{N} I_n \mathbf{f}_n,
$$

并取测试函数 $\mathbf{t}_m = \mathbf{f}_m$，则得到

$$
\sum_{n=1}^{N} I_n \langle \mathbf{f}_m, L(\mathbf{f}_n) \rangle
=
\langle \mathbf{f}_m, \mathbf{E}^{inc} \rangle.
$$

这对应离散线性系统

$$
\mathbf{Z} \mathbf{I} = \mathbf{V}.
$$

## 2. 矩阵元素计算

### 2.1 EFIE 矩阵

对表面 EFIE，典型阻抗矩阵元为

$$
Z_{mn}^{EFIE} =
j \omega \mu
\int_{S_m} \int_{S_n}
\mathbf{f}_m(\mathbf{r}) \cdot \mathbf{f}_n(\mathbf{r}') G(R)
\, dS' dS
+
\frac{1}{j \omega \epsilon}
\int_{S_m} \int_{S_n}
(\nabla \cdot \mathbf{f}_m)
(\nabla' \cdot \mathbf{f}_n)
G(R)
\, dS' dS.
$$

对 RWG-RWG 配对，实现中通常会先拆成四个支撑子三角形配对：

$$
Z_{mn}^{EFIE} = \sum_{i=1}^{2} \sum_{j=1}^{2}
\left(Z_{mn,L}^{(i,j)} + Z_{mn,\phi}^{(i,j)}\right).
$$

若沿用带符号边长记号，则每个局部子项可统一写成

$$
Z_{mn,L}^{(i,j)} =
j \omega \mu
\frac{\tilde l_{m,i} \tilde l_{n,j}}{4 A_{m,i} A_{n,j}}
\int_{T_{m,i}} \int_{T_{n,j}}
\left(\mathbf{r} - \mathbf{v}_{m,i}^{\mathrm{opp}}\right)
\cdot
\left(\mathbf{r}' - \mathbf{v}_{n,j}^{\mathrm{opp}}\right)
G(R)
\, dS' dS,
$$

$$
Z_{mn,\phi}^{(i,j)} =
\frac{1}{j \omega \epsilon}
\frac{\tilde l_{m,i}}{A_{m,i}}
\frac{\tilde l_{n,j}}{A_{n,j}}
\int_{T_{m,i}} \int_{T_{n,j}} G(R) \, dS' dS.
$$

这种写法与 EMMoMSuite 当前的局部几何加局部符号分解最接近。

### 2.2 MFIE 矩阵

MFIE 的离散形式包含主值积分与恒等项：

$$
Z_{mn}^{MFIE} =
\frac{1}{2} \int_{S_m} \mathbf{f}_m \cdot \mathbf{f}_n \, dS
-
\int_{S_m}
\mathbf{f}_m \cdot
\left(
\hat{\mathbf{n}} \times \int_{S_n} \mathbf{f}_n \times \nabla' G \, dS'
\right) dS.
$$

在 RWG-RWG 离散下，同样需要拆到局部支撑配对层面处理。局部法向、叉乘方向和正负支撑符号都应在统一的几何模板里表达，而不是再额外维护两套公式。

### 2.3 VIE 矩阵

对非磁介质，体积分方程常写为

$$
\mathbf{D}(\mathbf{r})
- \epsilon_0 (\epsilon_r(\mathbf{r}) - 1)
\int_V
\overline{\mathbf{G}}_e(\mathbf{r}, \mathbf{r}')
\cdot
\frac{\mathbf{D}(\mathbf{r}')}{\epsilon_0 \epsilon_r(\mathbf{r}')}
\, dV'
=
\epsilon_0 \mathbf{E}^{inc}(\mathbf{r}).
$$

若使用 SWG 基函数展开 $\mathbf{D} = \sum D_n \mathbf{f}_n$，则矩阵元涉及两个四面体上的体积分。实现时重点在于体格林函数核、散度项以及近邻四面体的积分精度控制。

### 2.4 面体耦合系统

对涂覆目标或导体和介质混合问题，常需联立面电流与体极化电流：

$$
\begin{bmatrix}
Z_{SS} & Z_{SV} \\
Z_{VS} & Z_{VV}
\end{bmatrix}
\begin{bmatrix}
I_S \\
I_V
\end{bmatrix}
=
\begin{bmatrix}
V_S \\
V_V
\end{bmatrix}.
$$

这里：

- $Z_{SS}$ 对应表面未知量之间的耦合。
- $Z_{VV}$ 对应体未知量之间的耦合。
- $Z_{SV}$ 与 $Z_{VS}$ 对应面体交叉耦合项。

混合维度积分往往比纯面或纯体积分更敏感，需要专门的近奇异处理策略。

## 3. 数值积分

### 3.1 高斯求积

对非奇异积分，最常用的是三角形或四面体上的高斯求积：

$$
\int_T g(\mathbf{r}) \, dS
\approx
A_T \sum_{i=1}^{N_{quad}} w_i g(\mathbf{r}_i).
$$

在实现中，常按 self、near、far 三类配对分别采用不同的积分核与积分阶数。对 RWG 表面路径，一个双面积分通常会展开为：

- 外层测试三角形高斯点循环。
- 内层源三角形高斯点循环。
- 针对近邻或奇异配对切换到专门的局部核或解析修正。

### 3.2 工程实现原则

为了与代码一一对应，推导文档最好先写支撑配对求和，再写每个局部配对的积分表达式。这样更容易与装配函数中的循环结构和数据访问方式对齐。

## 4. 奇异性处理

当源面元和观测面元重合或共享边顶点时，格林函数中的 $1 / R$ 项会导致奇异或近奇异行为。

### 4.1 奇异性减除

常见策略是把奇异核拆成“可解析项 + 光滑余项”：

$$
\int_T F(\mathbf{r}') \frac{1}{R} dS'
=
\int_T
\left(
\frac{F(\mathbf{r}')}{R} - \frac{F(\mathbf{r})}{R}
\right) dS'
+
F(\mathbf{r}) \int_T \frac{1}{R} dS'.
$$

第一项可数值积分，第二项交给解析公式。

### 4.2 Wilton 解析积分

对三角形上的静态奇异积分，Wilton 类公式能把结果写成对数项和反正切项的组合，是经典的 self 项处理工具。

### 4.3 Duffy 变换

Duffy 变换通过变量代换把三角形奇异积分映射到正方形参考域，并显式消去雅可比中的奇异因子。例如

$$
\int_0^1 \int_0^{1-x} \frac{f(x,y)}{\sqrt{x^2+y^2}} \, dy \, dx
\xrightarrow{x = u,\, y = u v}
\int_0^1 \int_0^1
\frac{f(u, u v)}{u \sqrt{1+v^2}} u \, dv \, du.
$$

变换后被积函数更平滑，可继续使用常规求积规则。

### 4.4 近奇异积分

当源和场很近但不重合时，积分虽然不真正发散，但会出现剧烈梯度变化。常见处理办法包括：

- 投影法：把源单元投影到观测点附近的局部坐标系中处理径向积分。
- 自适应细分：按距离和单元尺度递归细分源单元。
- 变量代换：对一维或准一维强变化核使用例如 sinh 变换来拉平积分核。

这些技术的目标都是把高梯度积分变回平滑可积的数值问题。
