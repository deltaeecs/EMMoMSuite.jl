# 后处理 (Post-Processing)

求解得到电流系数 $\bm{I}$ 后，可计算各类物理量。本章公式来自论文第 2 章
"后处理"一节（远场、雷达散射截面）。

## 1. 远场计算 (Far-Field Calculation)

远场条件下 $r \to \infty$，格林函数近似为（论文式 (2-69)）：

$$
G(R) \xrightarrow{r\to\infty}
\frac{e^{-{\rm j}k\left(r - \hat{\bm{r}} \cdot \bm{r}'\right)}}{4\pi r}
= \frac{e^{-{\rm j}kr}}{4\pi r}\, e^{{\rm j}k\hat{\bm{r}} \cdot \bm{r}'}
$$

对辐射远场做替换 $\nabla \to -{\rm j}k\hat{\bm{r}}$，代入
$\bm{E} = -{\rm j}\omega\bm{A} + \nabla\nabla\cdot\bm{A}/({\rm j}\omega\mu\varepsilon)$，
得到电场、磁场的远场形式（论文式 (2-70)~(2-71)）：

$$
\bm{E} = -{\rm j}k\eta\, \frac{e^{-{\rm j}kr}}{4\pi r}
\left(\hat{\bm{\theta}}\hat{\bm{\theta}} + \hat{\bm{\phi}}\hat{\bm{\phi}}\right) \cdot
\int_\Omega e^{{\rm j}k\hat{\bm{r}} \cdot \bm{r}'} \bm{J}(\bm{r}')\, d\Omega'
$$

$$
\bm{H} = -{\rm j}k\, \frac{e^{-{\rm j}kr}}{4\pi r}
\left(\hat{\bm{\phi}}\hat{\bm{\theta}} - \hat{\bm{\theta}}\hat{\bm{\phi}}\right) \cdot
\int_\Omega e^{{\rm j}k\hat{\bm{r}} \cdot \bm{r}'} \bm{J}(\bm{r}')\, d\Omega'
$$

投影算子 $\hat{\bm{\theta}}\hat{\bm{\theta}} + \hat{\bm{\phi}}\hat{\bm{\phi}}$ 等价于
$\overline{\bm{I}} - \hat{\bm{r}}\hat{\bm{r}}$，即只保留 $\theta$、$\phi$ 分量。
实际计算远场时一般忽略与距离有关的 $e^{-{\rm j}kr}/(4\pi r)$ 项。
`FarField.farField` 即按此公式实现，返回 `[2, Nθ, Nφ]` 的
$(E_\theta, E_\phi)$ 分量。

## 2. 雷达散射截面 (RCS)

RCS 表征目标在雷达波照射下产生回波的能力，定义为（论文式 (2-72)）：

$$
\sigma = \lim_{r\to\infty} \left[ 4\pi r^2
\frac{|\bm{E}^s|^2}{|\bm{E}^i|^2} \right]
$$

利用远场公式，单位幅度入射波下的单站/双站 RCS 可由辐射矢量
$\bm{N}(\theta,\phi) = \int_S \bm{J}\, e^{{\rm j}k\hat{\bm{r}}\cdot\bm{r}'} dS'$
计算：

$$
\sigma(\theta, \phi) = \frac{k^2\eta^2}{4\pi}
\left(|N_\theta|^2 + |N_\phi|^2\right)
$$

$$
\sigma_{dBsm} = 10\log_{10}\sigma
$$

`RCS.radarCrossSection` 返回分量 RCS、总 RCS（线性，m²）与 dBsm 三种形式。

## 3. 近场计算 (Near-Field Calculation)

近区必须使用精确的格林函数积分：

$$
\bm{E}(\bm{r}) = -{\rm j}\omega \bm{A}(\bm{r}) - \nabla \phi(\bm{r}), \qquad
\bm{H}(\bm{r}) = \frac{1}{\mu}\nabla \times \bm{A}(\bm{r})
$$

计算公式与阻抗矩阵元素类似，但观测点 $\bm{r}$ 不在源面上；观测点非常靠近
源面时需使用近奇异性处理（如论文第 2 章的奇异值提取方案或自适应细分）。
`NearField.jl` 与 `NearFieldAdvanced.jl` 提供相应实现。

## 4. 天线参数

- 方向性系数：$D(\theta, \phi) = 4\pi U(\theta, \phi) / P_{rad}$，
  $U$ 为辐射强度、$P_{rad}$ 为总辐射功率；
- 增益：$G(\theta, \phi) = \eta_{eff} D(\theta, \phi)$，计入损耗效率；
- 端口参数（S 参数、输入阻抗）见 `Ports` 模块与 `SParameterExtraction`。
