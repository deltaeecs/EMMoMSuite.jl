# 后处理 (Post-Processing)

求解得到电流系数 \$\mathbf{I}\$ 后，可计算各类感兴趣的物理量。

## 1. 远场计算 (Far-Field Calculation)

当观测距离 \$r \to \infty\$ 时，格林函数近似为：
\$\$
G(\mathbf{r}, \mathbf{r}') \approx \frac{e^{-jkr}}{4\pi r} e^{jk \hat{r} \cdot \mathbf{r}'}
\$\$
散射电场远场近似为：
\$\$
\mathbf{E}^{far}(\mathbf{r}) = \frac{e^{-jkr}}{r} \mathbf{F}(\theta, \phi)
\$\$
其中辐射矢量 \$\mathbf{F}(\theta, \phi)\$ 为：
\$\$
\mathbf{F}(\theta, \phi) = -j\omega\mu \frac{1}{4\pi} (\overline{\mathbf{I}} - \hat{r}\hat{r}) \cdot \int_S \mathbf{J}(\mathbf{r}') e^{jk \hat{r} \cdot \mathbf{r}'} dS'
\$\$
对于 RWG 基函数，积分转化为各个三角形上的傅里叶变换形式。

## 2. 雷达散射截面 (RCS)

RCS 定义为：
\$\$
\sigma(\theta, \phi) = \lim_{r \to \infty} 4\pi r^2 \frac{|\mathbf{E}^{scat}|^2}{|\mathbf{E}^{inc}|^2}
\$\$
对于单位幅度的入射平面波：
\$\$
\sigma_{dBsm} = 10 \log_{10} (4\pi |\mathbf{F}(\theta, \phi)|^2)
\$\$
通常计算双站 RCS (Bistatic RCS) 或单站 RCS (Monostatic RCS)。

## 3. 近场计算 (Near-Field Calculation)

在近区，必须使用精确的格林函数积分。
\$\$
\mathbf{E}(\mathbf{r}) = -j\omega \mathbf{A}(\mathbf{r}) - \nabla \Phi(\mathbf{r})
\$\$
\$\$
\mathbf{H}(\mathbf{r}) = \frac{1}{\mu} \nabla \times \mathbf{A}(\mathbf{r})
\$\$
计算公式与阻抗矩阵元素计算类似，但观测点 \$\mathbf{r}\$ 不在源面上。
**注意**：当观测点非常靠近源面时，需使用近奇异性处理技术（如自适应细分）。

## 4. 天线参数

### 4.1 方向性系数 (Directivity)
\$\$
D(\theta, \phi) = \frac{4\pi U(\theta, \phi)}{P_{rad}}
\$\$
其中 \$U\$ 为辐射强度，\$P_{rad}\$ 为总辐射功率。

### 4.2 增益 (Gain)
\$\$
G(\theta, \phi) = \eta_{eff} D(\theta, \phi)
\$\$
考虑了损耗效率。
