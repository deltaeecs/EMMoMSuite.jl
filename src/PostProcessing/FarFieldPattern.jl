"""
    FarFieldPatternModule — 远场方向图统一容器

提供 `FarFieldPattern` 结构体及其全部查询接口：
- `gain` / `gain_db`     — 增益（dBi）
- `hpbw`                 — 半功率波束宽度（度）
- `side_lobe_level`      — 副瓣电平（dB）
- `axial_ratio`          — 轴比（dB）
- `co_cross_decompose`   — Ludwig III 主/交叉极化分解
- `xpd`                  — 主/交叉极化比（dB）
"""
module FarFieldPatternModule

using LinearAlgebra
using StaticArrays

export FarFieldPattern
export gain, gain_db, hpbw, side_lobe_level
export axial_ratio, co_cross_decompose, xpd

# ─────────────────────────────────────────────────────────────────────────────
# 常数
# ─────────────────────────────────────────────────────────────────────────────
const ETA0_FP = 376.730313461  # 自由空间波阻抗 Ω

# ─────────────────────────────────────────────────────────────────────────────
# FarFieldPattern 结构体
# ─────────────────────────────────────────────────────────────────────────────

"""
    FarFieldPattern

多频率远场方向图容器。

# 字段
- `freqs`       : 频率向量 (Hz)，长度 Nf
- `theta`       : 极角向量 (rad)，长度 Nθ
- `phi`         : 方位角向量 (rad)，长度 Nφ
- `E_theta`     : θ 分量，形状 `(Nf, Nθ, Nφ)`
- `E_phi`       : φ 分量，形状 `(Nf, Nθ, Nφ)`
- `_gain`       : 懒惰缓存（增益 dBi）；初始 `nothing`
- `_axial_ratio`: 懒惰缓存（轴比 dB）；初始 `nothing`

# 构造
```julia
ff = FarFieldPattern(freqs, theta, phi, E_theta, E_phi)
```
"""
mutable struct FarFieldPattern
    freqs        :: Vector{Float64}
    theta        :: Vector{Float64}
    phi          :: Vector{Float64}
    E_theta      :: Array{ComplexF64, 3}
    E_phi        :: Array{ComplexF64, 3}
    _gain        :: Union{Nothing, Array{Float64, 3}}
    _axial_ratio :: Union{Nothing, Array{Float64, 3}}

    function FarFieldPattern(
        freqs  :: AbstractVector{<:Real},
        theta  :: AbstractVector{<:Real},
        phi    :: AbstractVector{<:Real},
        E_theta:: AbstractArray{<:Complex, 3},
        E_phi  :: AbstractArray{<:Complex, 3},
    )
        Nf = length(freqs)
        Nθ = length(theta)
        Nφ = length(phi)
        size(E_theta) == (Nf, Nθ, Nφ) ||
            throw(ArgumentError(
                "E_theta size $(size(E_theta)) ≠ (Nf=$Nf, Nθ=$Nθ, Nφ=$Nφ)"))
        size(E_phi) == (Nf, Nθ, Nφ) ||
            throw(ArgumentError(
                "E_phi size $(size(E_phi)) ≠ (Nf=$Nf, Nθ=$Nθ, Nφ=$Nφ)"))
        return new(
            collect(Float64, freqs),
            collect(Float64, theta),
            collect(Float64, phi),
            ComplexF64.(E_theta),
            ComplexF64.(E_phi),
            nothing,
            nothing,
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 内部：梯形积分辅助
# ─────────────────────────────────────────────────────────────────────────────

"""
    _trapz(y, x)

沿 1D 梯形数值积分 ∫ y dx。
"""
function _trapz(y::AbstractVector, x::AbstractVector)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("x and y must have same length"))
    s = zero(promote_type(eltype(y), eltype(x)))
    @inbounds for i in 1:(n-1)
        s += (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return s / 2
end

"""
    _trapz2d(U, theta, phi)

二维梯形积分 ∫∫ U(θ,φ) sin(θ) dθ dφ。
"""
function _trapz2d(U::AbstractMatrix, theta::AbstractVector, phi::AbstractVector)
    Nθ = length(theta)
    Nφ = length(phi)
    size(U) == (Nθ, Nφ) || throw(DimensionMismatch())

    # 沿 θ 积分（含 sinθ）
    inner = zeros(Float64, Nφ)
    @inbounds for j in 1:Nφ
        col_sinθ = [U[i,j] * sin(theta[i]) for i in 1:Nθ]
        inner[j] = _trapz(col_sinθ, theta)
    end
    return _trapz(inner, phi)
end

# ─────────────────────────────────────────────────────────────────────────────
# 增益
# ─────────────────────────────────────────────────────────────────────────────

"""
    gain(ff::FarFieldPattern; freq_idx=1) → Matrix{Float64}

计算增益方向图（dBi），返回 `(Nθ, Nφ)` 矩阵。

**公式**:
- 辐射强度 U(θ,φ) = (|E_θ|² + |E_φ|²) / (2η₀)
- P_rad = ∫∫ U sinθ dθ dφ（梯形法）
- D(θ,φ) = 4π U / P_rad（方向性系数，线性量）
- gain_dBi = 10·log10(D)
"""
function gain(ff::FarFieldPattern; freq_idx::Int=1)
    1 <= freq_idx <= length(ff.freqs) ||
        throw(ArgumentError("freq_idx=$freq_idx out of range [1,$(length(ff.freqs))]"))

    Eθ = ff.E_theta[freq_idx, :, :]   # Nθ × Nφ
    Eφ = ff.E_phi[freq_idx, :, :]

    U = (abs2.(Eθ) .+ abs2.(Eφ)) ./ (2 * ETA0_FP)

    P_rad = _trapz2d(U, ff.theta, ff.phi)

    if P_rad <= 0
        return fill(-Inf, size(U))
    end

    D = (4π * U) ./ P_rad
    return @. 10 * log10(max(D, eps(Float64)))
end

"""
    gain_db(ff::FarFieldPattern, θ, φ; freq_idx=1) → Float64

在 (θ, φ) 方向插值增益（最近邻）。
"""
function gain_db(ff::FarFieldPattern, θ::Real, φ::Real; freq_idx::Int=1)
    θ_arr = ff.theta
    φ_arr = ff.phi
    iθ = argmin(abs.(θ_arr .- θ))
    iφ = argmin(abs.(φ_arr .- φ))
    g = gain(ff; freq_idx)
    return g[iθ, iφ]
end

# ─────────────────────────────────────────────────────────────────────────────
# HPBW / SLL
# ─────────────────────────────────────────────────────────────────────────────

"""
    hpbw(ff::FarFieldPattern; plane::Symbol=:E, freq_idx=1) → Float64

计算半功率波束宽度（HPBW，度）。

- `plane=:E` — 在 φ=0（或最近）切面上计算
- `plane=:H` — 在 θ=π/2（或最近）切面，沿 φ 方向扫描
"""
function hpbw(ff::FarFieldPattern; plane::Symbol=:E, freq_idx::Int=1)
    g_dBi = gain(ff; freq_idx)  # Nθ × Nφ

    if plane === :E || plane === :H
        if plane === :E
            # φ=0 切面 → 沿 θ 方向
            iφ = argmin(abs.(ff.phi .- 0.0))
            pattern = g_dBi[:, iφ]
            angles  = ff.theta .* (180 / π)
        else   # :H — θ=π/2 切面，沿 φ 方向
            iθ = argmin(abs.(ff.theta .- π/2))
            pattern = g_dBi[iθ, :]
            angles  = ff.phi .* (180 / π)
        end

        peak_idx = argmax(pattern)
        peak_val = pattern[peak_idx]
        half_pow = peak_val - 3.0

        # 从峰值向两侧找 3dB 点
        left_idx  = peak_idx
        right_idx = peak_idx
        for i in (peak_idx-1):-1:1
            pattern[i] < half_pow && break
            left_idx = i
        end
        for i in (peak_idx+1):length(pattern)
            pattern[i] < half_pow && break
            right_idx = i
        end

        return Float64(angles[right_idx] - angles[left_idx])
    else
        throw(ArgumentError("plane must be :E or :H, got $plane"))
    end
end

"""
    side_lobe_level(ff::FarFieldPattern; plane::Symbol=:E, freq_idx=1) → Float64

计算副瓣电平（SLL，dB），相对于主瓣峰值的最高旁瓣。
"""
function side_lobe_level(ff::FarFieldPattern; plane::Symbol=:E, freq_idx::Int=1)
    g_dBi = gain(ff; freq_idx)

    if plane === :E
        iφ = argmin(abs.(ff.phi .- 0.0))
        pattern = g_dBi[:, iφ]
    elseif plane === :H
        iθ = argmin(abs.(ff.theta .- π/2))
        pattern = g_dBi[iθ, :]
    else
        throw(ArgumentError("plane must be :E or :H"))
    end

    peak_val = maximum(pattern)
    peak_idx = argmax(pattern)

    # 找主瓣范围：从峰值两边找第一个局部最小值作为边界
    left_bound  = _find_local_min_left(pattern, peak_idx)
    right_bound = _find_local_min_right(pattern, peak_idx)

    # 旁瓣区域
    sidelobe_vals = vcat(pattern[1:left_bound], pattern[right_bound:end])
    isempty(sidelobe_vals) && return -Inf
    sll = maximum(sidelobe_vals) - peak_val
    return sll
end

function _find_local_min_left(pat, peak)
    i = peak - 1
    while i > 1
        pat[i] <= pat[i-1] && return i
        i -= 1
    end
    return 1
end

function _find_local_min_right(pat, peak)
    n = length(pat)
    i = peak + 1
    while i < n
        pat[i] <= pat[i+1] && return i
        i += 1
    end
    return n
end

# ─────────────────────────────────────────────────────────────────────────────
# 轴比（Axial Ratio）
# ─────────────────────────────────────────────────────────────────────────────

"""
    axial_ratio(ff::FarFieldPattern; freq_idx=1) → Matrix{Float64}

计算极化椭圆轴比（dB），返回 `(Nθ, Nφ)` 矩阵。

**公式**（IEEE Std 149-1979, Eq. A-5）:
- E_R = (E_θ - j·E_φ) / √2   （右旋圆极化分量）
- E_L = (E_θ + j·E_φ) / √2   （左旋圆极化分量）
- AR  = (|E_R| + |E_L|) / ||E_R| - |E_L||
- AR_dB = 20 · log10(AR)

当 |E_R| = |E_L| 时（线极化） AR→∞；当两者之一为零时（圆极化）AR = 0 dB。
"""
function axial_ratio(ff::FarFieldPattern; freq_idx::Int=1)
    Eθ = ff.E_theta[freq_idx, :, :]
    Eφ = ff.E_phi[freq_idx, :, :]

    E_R = (Eθ .- im .* Eφ) ./ sqrt(2)   # RHCP component
    E_L = (Eθ .+ im .* Eφ) ./ sqrt(2)   # LHCP component

    aR = abs.(E_R)
    aL = abs.(E_L)

    AR = similar(aR, Float64)
    @inbounds for idx in eachindex(aR)
        num = aR[idx] + aL[idx]
        den = abs(aR[idx] - aL[idx])
        if den < 1e-12 * num
            AR[idx] = Inf   # 线极化
        else
            AR[idx] = num / den
        end
    end

    return @. 20 * log10(max(AR, eps(Float64)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Ludwig III 极化分解
# ─────────────────────────────────────────────────────────────────────────────

"""
    co_cross_decompose(ff::FarFieldPattern; freq_idx=1)
        → NamedTuple{:co, :cross}

Ludwig 第三定义的主极化/交叉极化分解，返回 `(Nθ, Nφ)` 复数矩阵对。

**Ludwig III**（IEEE Std 1720-2012 定义）:
对每个观测方向 (θ, φ)：
- E_co    = E_θ · cos(φ) - E_φ · sin(φ)   [垂直/V 极化主分量]
- E_cross = E_θ · sin(φ) + E_φ · cos(φ)   [水平/H 极化交叉分量]

**说明**: 在 φ=0 切面，E_co = E_θ（主极化即 θ 分量），E_cross = E_φ（交叉分量）。

**性质**: co + cross 功率之和 = 总辐射功率（变换为酉变换）。
"""
function co_cross_decompose(ff::FarFieldPattern; freq_idx::Int=1)
    Eθ = ff.E_theta[freq_idx, :, :]
    Eφ = ff.E_phi[freq_idx,   :, :]

    Nθ, Nφ = size(Eθ)
    E_co    = similar(Eθ, ComplexF64)
    E_cross = similar(Eφ, ComplexF64)

    @inbounds for j in 1:Nφ
        c = cos(ff.phi[j])
        s = sin(ff.phi[j])
        for i in 1:Nθ
            E_co[i,j]    = Eθ[i,j] * c - Eφ[i,j] * s
            E_cross[i,j] = Eθ[i,j] * s + Eφ[i,j] * c
        end
    end

    return (; co=E_co, cross=E_cross)
end

"""
    xpd(ff::FarFieldPattern; freq_idx=1) → Float64

主极化与交叉极化功率之比（XPD，dB）。

XPD = 10 · log10(P_co / P_cross)

当 P_cross ≈ 0 时（纯主极化），XPD → +Inf；返回 300.0 dB 作为上限。
"""
function xpd(ff::FarFieldPattern; freq_idx::Int=1)
    decomp = co_cross_decompose(ff; freq_idx)
    P_co    = sum(abs2, decomp.co)
    P_cross = sum(abs2, decomp.cross)

    if P_cross < 1e-30 * P_co
        return 300.0   # 实际 Inf 用 300 dB 表示
    end
    return 10 * log10(P_co / P_cross)
end

end  # module FarFieldPatternModule
