"""
    PMCHWModule — PMCHWT 均匀介质体表面积分方程

实现 Poggio-Miller-Chang-Harrington-Wu-Tsai (PMCHWT) 方法，用于均匀介质体的散射/辐射分析。

# 数学背景 (Gibson §3.6.3.1)

对于平均介质体（外部 R₀, k₀, η₀；内部 R₁, k₁, η₁），联立内外 EFIE+MFIE
并利用边界条件 J₁=-J₂, M₁=-M₂，化简得 2N×2N 线性系统：

```
[Z^EJ   Z^EM ] [J]   [V_E]
[Z^HJ   Z^HM ] [M] = [V_H]
```

其中：
- Z^EJ  = L(k₀, η₀) + L(k₁, η₁)    [因子 jkη/(16π)，等同于两个EFIE之和]
- Z^EM  = K(k₀) + K(k₁)             [纯PV K算子，因子 1/(16π)，无质量矩阵项]
- Z^HJ  = -(K(k₀) + K(k₁)) = -Z^EM  [K的负版本]
- Z^HM  = Lₑ(k₀, η₀) + Lₑ(k₁, η₁) [因子 jk/(η·16π)，"inverted-η" EFIE]

结构不变量：Z^EM + Z^HJ = 0（精确成立）。
"""
module PMCHWModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..MFIEModule
using ..Impedance
using StaticArrays
using LinearAlgebra

import ..CoreModule: assemble_impedance_matrix

export PMCHW, assemble_impedance_matrix

# ─────────────────────────────────────────────────────────────────────────────
# 常数
# ─────────────────────────────────────────────────────────────────────────────
const _C0    = 299792458.0   # 真空光速 (m/s)
const _MU0   = 4π * 1e-7     # 真空磁导率 (H/m)
const _EPS0  = 1.0 / (_C0^2 * _MU0)  # 真空介电常数 (F/m)
const _ETA0  = sqrt(_MU0 / _EPS0)    # 自由空间波阻抗 (≈ 377 Ω)

# ─────────────────────────────────────────────────────────────────────────────
# 结构体
# ─────────────────────────────────────────────────────────────────────────────

"""
    PMCHW{FT, CT} <: AbstractIntegralOperator

Poggio-Miller-Chang-Harrington-Wu-Tsai (PMCHWT) 表面积分方程算子。

用于求解均匀介质体（相对介电常数 `eps_r`，相对磁导率 `mu_r`）散射/辐射问题。

# 字段
- `freq`:   工作频率 (Hz)
- `k0`:     自由空间波数 k₀ = 2πf/c₀
- `eta0`:   自由空间波阻抗 η₀  ≈ 377 Ω
- `k1`:     内部介质波数 k₁ = k₀ √(εᵣ μᵣ)（复数，支持损耗）
- `eta1`:   内部介质波阻抗 η₁ = η₀ √(μᵣ/εᵣ)（复数，支持损耗）
- `eps_r`:  相对介电常数（可为复数，支持有损介质）
- `mu_r`:   相对磁导率（可为复数）
"""
struct PMCHW{FT<:AbstractFloat,CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    k0::FT
    eta0::FT
    k1::CT       # 内部介质波数（复数：实部传播，虚部衰减）
    eta1::CT      # 内部介质波阻抗
    eps_r::CT
    mu_r::CT
end

"""
    PMCHW(freq, eps_r, mu_r=1.0)

构造 PMCHWT 算子。

# 参数
- `freq`:  工作频率 (Hz)
- `eps_r`: 相对介电常数（实数=无损，复数=有损）
- `mu_r`:  相对磁导率（默认 1.0）
"""
function PMCHW(freq::FT, eps_r_in, mu_r_in = 1.0) where {FT<:AbstractFloat}
    CT = Complex{FT}
    k0   = FT(2π * freq / _C0)
    eta0 = FT(_ETA0)

    eps_r = CT(eps_r_in)
    mu_r  = CT(mu_r_in)

    k1   = CT(k0) * sqrt(eps_r * mu_r)
    eta1 = CT(eta0) * sqrt(mu_r / eps_r)

    return PMCHW{FT,CT}(freq, k0, eta0, k1, eta1, eps_r, mu_r)
end

# ─────────────────────────────────────────────────────────────────────────────
# 内部辅助：构造各子块算子
# ─────────────────────────────────────────────────────────────────────────────

"""
    _l_block_operator(k_real, eta_real, k_c, eta_c, mode) → EFIE

构造局部 L 算子（EFIE-like），供 PMCHW 内部使用。

- `k_real`:   用于 Green 函数的实数波数（lossless 积分核参数）
- `eta_real`: 仅用于 EFIE 记录字段
- `k_c`:      完整复数波数（用于计算 factor）
- `eta_c`:    完整复数波阻抗
- `mode`:     :EJ（factor = jk η / 16π）或 :HM（factor = jk / (η·16π)）
"""
function _l_block_operator(
    k_real::FT,
    eta_real::FT,
    k_c::CT,
    eta_c::CT,
    mode::Symbol,
) where {FT,CT<:Complex}
    π16 = FT(16π)
    factor = if mode === :EJ
        CT(im * k_c * eta_c / π16)   # jk η / (16π) → Z^EJ 因子 = jωμ 相关
    else
        CT(im * k_c / (eta_c * π16)) # jk / (η·16π) → Z^HM 因子 = jωε 相关
    end
    return efie_from_keta(k_real, eta_real, factor)
end

"""
    _k_real(k_complex)

从复数波数提取用于 Green 函数的实数部分。
- 对实数（无损）：直接返回
- 对复数（有损）：返回实部的绝对值（忽略损耗对相位特性的影响，仅为近似）
"""
function _k_real(k_c::CT) where {CT<:Complex}
    r = real(k_c)
    abs(r) < 1e-40 * abs(k_c) && (r = abs(k_c))
    return r
end

# ─────────────────────────────────────────────────────────────────────────────
# 矩阵装配
# ─────────────────────────────────────────────────────────────────────────────

"""
    assemble_impedance_matrix(pmchw::PMCHW, basis::RWGBasis) → Matrix{Complex}

装配 PMCHWT 全阻抗矩阵 Z（大小 2N×2N）。

```
Z = [Z^EJ   Z^EM ]
    [Z^HJ   Z^HM ]
```

- Z^EJ [1:N,   1:N  ]:  jωμ₀ L(k₀) + jωμ₁ L(k₁)
- Z^EM [1:N,   N+1:2N]: K(k₀) + K(k₁)      [纯 PV K，无质量矩阵]
- Z^HJ [N+1:2N, 1:N ]: -K(k₀) − K(k₁) = −Z^EM
- Z^HM [N+1:2N, N+1:2N]: jωε₀ L(k₀) + jωε₁ L(k₁)

# 精度说明
对于**有损介质**（eps_r 有虚部），Z^EJ 和 Z^HM 中的 Green 函数使用 Re(k₁)
作为近似（忽略指数衰减项），因子仍使用完整复数 k₁, η₁。
无损介质（eps_r 实数）结果精确。
"""
function assemble_impedance_matrix(pmchw::PMCHW{FT,CT}, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    N = num_basis(basis)

    # ── 波数与波阻抗 ──
    k0     = pmchw.k0
    eta0   = pmchw.eta0
    k0_c   = CT(k0)
    eta0_c = CT(eta0)

    k1_c   = pmchw.k1
    eta1_c = pmchw.eta1
    k1_r   = FT(_k_real(k1_c))     # 实数近似，用于 Green 函数
    eta1_r = FT(abs(real(eta1_c))) # 实数近似，用于 EFIE 记录字段

    Z = zeros(CT, 2N, 2N)

    # ── Z^EJ = L(k₀) + L(k₁),  factor = jk η / (16π) ──────────────────────
    efie_ej0 = _l_block_operator(k0, eta0, k0_c, eta0_c, :EJ)
    efie_ej1 = _l_block_operator(k1_r, eta1_r, k1_c, eta1_c, :EJ)
    Z[1:N, 1:N]       .+= assemble_impedance_matrix(efie_ej0, basis)
    Z[1:N, 1:N]       .+= assemble_impedance_matrix(efie_ej1, basis)

    # ── Z^HM = Lₑ(k₀) + Lₑ(k₁),  factor = jk / (η·16π) ───────────────────
    efie_hm0 = _l_block_operator(k0, eta0, k0_c, eta0_c, :HM)
    efie_hm1 = _l_block_operator(k1_r, eta1_r, k1_c, eta1_c, :HM)
    Z[N+1:2N, N+1:2N] .+= assemble_impedance_matrix(efie_hm0, basis)
    Z[N+1:2N, N+1:2N] .+= assemble_impedance_matrix(efie_hm1, basis)

    # ── Z^EM = K(k₀) + K(k₁)，Z^HJ = -Z^EM ─────────────────────────────────
    K0 = assemble_K_offdiag(basis, k0)
    K1 = assemble_K_offdiag(basis, k1_r)  # lossless approx for Green kernel
    K_total = K0 + K1

    Z[1:N,   N+1:2N] .=  K_total
    Z[N+1:2N, 1:N]   .= -K_total

    return Z
end

end # module PMCHWModule
