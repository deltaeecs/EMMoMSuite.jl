module Interpolation

using StaticArrays
using SparseArrays
using FastGaussQuadrature
using FFTW
using LinearAlgebra
using ....Geometry

export AbstractPolesInfo, AbstractInterpInfo
export GLPolesInfo, LagrangeInterpInfo, FFTGLPolesInfo, FFTInterpInfo
export levelIntegralInfoCal, interpolationCSCMatCal, truncation_kernel
export interpolate, anterpolate, interpolate!, anterpolate!
export fft_interp_phi, fft_anterp_phi
export fft_interp_phi!, fft_anterp_phi!
export fft_interp_phi_batch!, fft_anterp_phi_batch!
export interp_type

abstract type AbstractPolesInfo{FT<:AbstractFloat} end
abstract type AbstractInterpInfo{IT<:Integer,FT<:Real} end

function interpolate end
function anterpolate end
function interpolate! end
function anterpolate! end

const NBDIGITS = 9.0 # Default value

"""
    truncation_kernel(rel_l) -> L

MLFMA 转移函数截断项数经验公式（论文式 (2-42)）：

```math
\\tau(l) \\approx 1.73\\, k a_l + 2.16\\, d_0^{2/3} (k a_l)^{1/3}, \\qquad d_0 = 3
```

实现输入 `rel_l = a_l / λ`（盒子边长以波长为单位），利用
`k a_l = 2π a_l / λ = 2π rel_l` 与 `1.73·2π ≈ 2π√3` 改写为：

```math
\\tau = 2\\pi \\sqrt{3}\\, rel_l + 2.16\\, d_0^{2/3} (2\\pi\\, rel_l)^{1/3}
```

注意：实现中的精度参数 `NBDIGITS = 9.0` 对应公式中的 `d_0`（论文推荐
`d_0 = 3`），取 9 时第二项约为推荐值的 `(9/3)^{2/3} ≈ 2.08` 倍，截断更保守。
`levelIntegralInfoCal` 使用 `ceil` 取整为整数截断项。
"""
function truncation_kernel(rel_l)
    return 2π * rel_l * sqrt(3) + 2.16 * NBDIGITS^(2.0 / 3.0) * (2π * rel_l)^(1 / 3)
end

"""
    GLPolesInfo{FT}

Gauss-Legendre Poles Information.
"""
struct GLPolesInfo{FT<:Real} <: AbstractPolesInfo{FT}
    Xθs::Vector{FT}
    Xϕs::Vector{FT}
    Wθϕs::Vector{FT}
    r̂sθsϕs::Vector{r̂θϕInfo{FT}}
end

"""
    interp_type(::AbstractPolesInfo) -> Type

采样信息类型 -> 对应插值矩阵类型。GL 网格用两段 Lagrange；Lebedev（LbPolesInfo）
在 Lebedev/LVI.jl 中扩展为一步训练/球谐插值（LbTrainedInterp1tepInfo）。
"""
interp_type(::GLPolesInfo{FT}) where {FT} = LagrangeInterpInfo{Int,FT}

"""
    LagrangeInterpInfo{IT, FT}

Lagrange Interpolation Information.
"""
mutable struct LagrangeInterpInfo{IT<:Integer,FT<:Real} <: AbstractInterpInfo{IT,FT}
    θCSC::SparseMatrixCSC{FT,IT}
    ϕCSC::SparseMatrixCSC{FT,IT}
    θCSCT::SparseMatrixCSC{FT,IT}
    ϕCSCT::SparseMatrixCSC{FT,IT}

    function LagrangeInterpInfo{IT,FT}() where {IT<:Integer,FT<:Real}
        new{IT,FT}()
    end
end

"""
    FFTGLPolesInfo{FT} <: AbstractPolesInfo{FT}

GL 采样信息的标记包装：表示该层的层间插值使用 φ 方向 FFT 谱插值（方案 P1），
θ 方向仍为 Lagrange。字段访问自动转发到内部 `GLPolesInfo`。
"""
struct FFTGLPolesInfo{FT<:Real} <: AbstractPolesInfo{FT}
    inner::GLPolesInfo{FT}
end

Base.getproperty(p::FFTGLPolesInfo, s::Symbol) =
    s === :inner ? getfield(p, :inner) : getproperty(getfield(p, :inner), s)
Base.propertynames(p::FFTGLPolesInfo) = propertynames(p.inner)

"""
    FFTInterpInfo{IT,FT} <: AbstractInterpInfo{IT,FT}

φ 方向 FFT 谱插值 + θ 方向 Lagrange 的层间插值信息。
`θCSC`/`θCSCT` 为 θ 方向 Lagrange 插值矩阵及其转置（构造方式与 `LagrangeInterpInfo` 相同）；
φ 方向的插值/反插值由 `fft_interp_phi`/`fft_anterp_phi` 在调用处按相邻层极点数完成。
`nθ`/`M1`/`M2` 为子层 θ 行数与子/父层 φ 采样数；`corr` 为半格网格相位修正；
`workspaces` 缓存每线程的 FFTW 计划与缓冲（`fft_interp_phi!`/`fft_anterp_phi!` 使用）。
"""
mutable struct FFTInterpInfo{IT<:Integer,FT<:Real} <: AbstractInterpInfo{IT,FT}
    θCSC::SparseMatrixCSC{FT,IT}
    θCSCT::SparseMatrixCSC{FT,IT}
    nθ::IT
    M1::IT
    M2::IT
    nCubes::IT
    corr::Vector{ComplexF64}
    corrL::Vector{ComplexF64}
    corrH::Vector{ComplexF64}
    workspaces::Vector{Any}
    bworkspaces::Vector{Any}
end

function FFTInterpInfo(
    θCSC::SparseMatrixCSC{FT,IT},
    θCSCT::SparseMatrixCSC{FT,IT},
    nθ::IT,
    M1::IT,
    M2::IT,
    nCubes::IT,
) where {IT<:Integer,FT<:Real}
    half = M1 ÷ 2
    corr = zeros(ComplexF64, M2)
    corr[1:half] = exp.(im .* (0:(half-1)) .* (π / M2 - π / M1))
    corr[(M2-half+1):M2] = exp.(im .* ((half:(M1-1)) .- M1) .* (π / M2 - π / M1))
    corrL = corr[1:half]
    corrH = corr[(M2-half+1):M2]
    workspaces = [nothing for _ in 1:Threads.nthreads()]
    bworkspaces = [nothing for _ in 1:Threads.nthreads()]
    return FFTInterpInfo{IT,FT}(
        θCSC,
        θCSCT,
        nθ,
        M1,
        M2,
        nCubes,
        corr,
        corrL,
        corrH,
        workspaces,
        bworkspaces,
    )
end

"""
    interp_type(::FFTGLPolesInfo) -> Type

FFT 谱插值路径对应 `FFTInterpInfo`。
"""
interp_type(::FFTGLPolesInfo{FT}) where {FT<:Real} = FFTInterpInfo{Int,FT}

"""
    levelIntegralInfoCal(levelCubeEdgel, Val(:FFTSpectral); λ=1.0, L_min=0)

与 GL 采样相同的截断与采样点计算，但返回 `FFTGLPolesInfo` 包装，
标记层间插值走 φ 方向 FFT 谱插值（方案 P1）。
"""
function levelIntegralInfoCal(
    levelCubeEdgel::FT,
    ::Val{:FFTSpectral};
    λ = 1.0,
    L_min::Int = 0,
) where {FT<:Real}
    L, poles = levelIntegralInfoCal(levelCubeEdgel; λ = λ, L_min = L_min)
    return L, FFTGLPolesInfo(poles)
end

"""
    fft_interp_phi(x, nθ, M1, M2) -> y

φ 方向带限 FFT 谱插值（上采样 M1 → M2，比例因子 M2/M1）。
`x` 长度 `nθ*M1`，按 (θ 内、φ 外) 列主序展平；对每个 θ 行执行
DFT → 频域补零 → IDFT。MLFMA 的 φ 采样网格为半格偏置（`φ_j=(j-1/2)·2π/M`），
因此补零前需乘相位修正 `exp(ikπ(1/M2-1/M1))`（k 为有符号频率），
带限信号（带宽 < M1/2）下达到机器精度。
"""
function fft_interp_phi(
    x::AbstractVector,
    nθ::Integer,
    M1::Integer,
    M2::Integer,
) 
    length(x) == nθ * M1 || throw(ArgumentError("length(x) != nθ*M1"))
    A = reshape(ComplexF64.(x), nθ, M1)
    F = fft(A, 2)
    scale = M2 / M1
    half = M1 ÷ 2
    corr = zeros(ComplexF64, M2)
    corr[1:half] = exp.(im .* (0:(half-1)) .* (π / M2 - π / M1))
    corr[(M2-half+1):M2] = exp.(im .* ((half:(M1-1)) .- M1) .* (π / M2 - π / M1))
    Fp = hcat(
        view(F, :, 1:half),
        zeros(ComplexF64, nθ, M2 - M1),
        view(F, :, (half+1):M1),
    )
    B = scale .* Fp .* reshape(corr, 1, M2)
    return vec(ifft(B, 2))
end

"""
    fft_anterp_phi(x, nθ, M2, M1) -> y

φ 方向 FFT 反插值（下采样 M2 → M1）——即 `fft_interp_phi` 矩阵的**转置（伴随）**，
与 Lagrange 路径 `ϕCSCT` 的语义一致（Ergul §3.3 反插值 = 转置插值）。
对每个 θ 行执行 `scale·fft( 截断( c·ifft(x) ) )`，比例因子 M2/M1；
`x` 长度 `nθ*M2`，按 (θ 内、φ 外) 列主序展平。
"""
function fft_anterp_phi(
    x::AbstractVector,
    nθ::Integer,
    M2::Integer,
    M1::Integer,
) 
    length(x) == nθ * M2 || throw(ArgumentError("length(x) != nθ*M2"))
    A = reshape(ComplexF64.(x), nθ, M2)
    u = ifft(A, 2)
    scale = M2 / M1
    half = M1 ÷ 2
    corr = zeros(ComplexF64, M2)
    corr[1:half] = exp.(im .* (0:(half-1)) .* (π / M2 - π / M1))
    corr[(M2-half+1):M2] = exp.(im .* ((half:(M1-1)) .- M1) .* (π / M2 - π / M1))
    v = u .* reshape(corr, 1, M2)
    w = hcat(view(v, :, 1:half), view(v, :, (M2-half+1):M2))
    return vec(scale .* fft(w, 2))
end

function _fft_ws(info::FFTInterpInfo)
    ws = info.workspaces
    tid = Threads.threadid()
    if length(ws) < tid || ws[tid] === nothing
        buf1 = zeros(ComplexF64, info.nθ, info.M1)
        buf2 = zeros(ComplexF64, info.nθ, info.M2)
        p1f = FFTW.plan_fft!(buf1, 2)
        p1i = FFTW.plan_ifft!(buf1, 2)
        p2f = FFTW.plan_fft!(buf2, 2)
        p2i = FFTW.plan_ifft!(buf2, 2)
        if length(ws) < tid
            resize!(ws, tid)
        end
        ws[tid] = (buf1, buf2, p1f, p2i, p2f, p1i)
    end
    return ws[tid]
end

"""
    fft_interp_phi!(out, x, info::FFTInterpInfo)

零分配版本：φ 方向 FFT 谱插值（子层 M1 → 父层 M2），使用 `info` 缓存的
FFTW 计划、相位修正与每线程缓冲。`x` 长度 `nθ*M1`，`out` 长度 `nθ*M2`。
"""
function fft_interp_phi!(out::AbstractVector{CT}, x::AbstractVector, info::FFTInterpInfo) where {CT}
    nθ, M1, M2 = info.nθ, info.M1, info.M2
    length(out) == nθ * M2 || throw(ArgumentError("length(out) != nθ*M2"))
    buf1, buf2, p1f, p2i, p2f, p1i = _fft_ws(info)
    copyto!(buf1, 1, x, 1, nθ * M1)
    p1f * buf1
    half = M1 ÷ 2
    scale = M2 / M1
    corrL, corrH = info.corrL, info.corrH
    fill!(buf2, 0)
    @views buf2[:, 1:half] .= scale .* buf1[:, 1:half] .* transpose(corrL)
    @views buf2[:, (M2-half+1):M2] .=
        scale .* buf1[:, (half+1):M1] .* transpose(corrH)
    p2i * buf2
    copyto!(out, 1, vec(buf2), 1, nθ * M2)
    return out
end

"""
    fft_anterp_phi!(out, x, info::FFTInterpInfo)

零分配版本：φ 方向 FFT 反插值（父层 M2 → 子层 M1，插值矩阵的转置），
使用 `info` 缓存的 FFTW 计划、相位修正与每线程缓冲。
`x` 长度 `nθ*M2`，`out` 长度 `nθ*M1`。
"""
function fft_anterp_phi!(out::AbstractVector{CT}, x::AbstractVector, info::FFTInterpInfo) where {CT}
    nθ, M1, M2 = info.nθ, info.M1, info.M2
    length(out) == nθ * M1 || throw(ArgumentError("length(out) != nθ*M1"))
    buf1, buf2, p1f, p2i, p2f, p1i = _fft_ws(info)
    copyto!(buf2, 1, x, 1, nθ * M2)
    p2i * buf2
    half = M1 ÷ 2
    scale = M2 / M1
    buf2 .*= transpose(info.corr)
    @views buf1[:, 1:half] .= scale .* buf2[:, 1:half]
    @views buf1[:, (half+1):M1] .= scale .* buf2[:, (M2-half+1):M2]
    p1f * buf1
    copyto!(out, 1, vec(buf1), 1, nθ * M1)
    return out
end

"""
    fft_interp_phi_batch!(out, agg, info::FFTInterpInfo)

批量 φ 方向 FFT 谱插值：所有子盒、两个极化的模式一次批量完成
（子层 M1 → 父层 M2）。`agg` 形状 `(nθ*M1, 2, nCubes)`，`out` 形状 `(nθ*M2, 2, nCubes)`；
结果与逐调用 `fft_interp_phi` 一致（机器精度）。
"""
function fft_interp_phi_batch!(out, agg, info::FFTInterpInfo)
    nθ, M1, M2 = info.nθ, info.M1, info.M2
    nCubes = size(agg, 3)
    Ws, Wb, pWsf, pWbi = _bws(info)
    for c in 1:nCubes, p in 1:2
        copyto!(view(Ws, :, :, p, c), 1, view(agg, :, p, c), 1, nθ * M1)
    end
    pWsf * Ws
    half = M1 ÷ 2
    scale = M2 / M1
    fill!(Wb, 0)
    @views Wb[:, 1:half, :, :] .=
        scale .* Ws[:, 1:half, :, :] .* reshape(info.corrL, 1, half, 1, 1)
    @views Wb[:, (M2-half+1):M2, :, :] .=
        scale .* Ws[:, (half+1):M1, :, :] .* reshape(info.corrH, 1, half, 1, 1)
    pWbi * Wb
    for c in 1:nCubes, p in 1:2
        copyto!(view(out, :, p, c), 1, view(Wb, :, :, p, c), 1, nθ * M2)
    end
    return out
end

"""
    fft_anterp_phi_batch!(out, temp, info::FFTInterpInfo)

批量 φ 方向 FFT 反插值（父层 M2 → 子层 M1，插值矩阵转置）：
所有子实例一次批量完成。`temp` 形状 `(nθ*M2, 2, nInst)`，`out` 形状 `(nθ*M1, 2, nInst)`。
"""
function fft_anterp_phi_batch!(out, temp, info::FFTInterpInfo)
    nθ, M1, M2 = info.nθ, info.M1, info.M2
    nInst = size(temp, 3)
    Ws, Wb, pWsf, pWbi = _bws(info)
    for c in 1:nInst, p in 1:2
        copyto!(view(Wb, :, :, p, c), 1, view(temp, :, p, c), 1, nθ * M2)
    end
    pWbi * Wb
    Wb .*= reshape(info.corr, 1, M2, 1, 1)
    half = M1 ÷ 2
    scale = M2 / M1
    @views Ws[:, 1:half, :, :] .= scale .* Wb[:, 1:half, :, :]
    @views Ws[:, (half+1):M1, :, :] .= scale .* Wb[:, (M2-half+1):M2, :, :]
    pWsf * Ws
    for c in 1:nInst, p in 1:2
        copyto!(view(out, :, p, c), 1, view(Ws, :, :, p, c), 1, nθ * M1)
    end
    return out
end

function _bws(info::FFTInterpInfo)
    ws = info.bworkspaces
    tid = Threads.threadid()
    if length(ws) < tid || ws[tid] === nothing
        nθ, M1, M2, nCubes = info.nθ, info.M1, info.M2, info.nCubes
        Ws = zeros(ComplexF64, nθ, M1, 2, nCubes)
        Wb = zeros(ComplexF64, nθ, M2, 2, nCubes)
        pWsf = FFTW.plan_fft!(Ws, 2)
        pWbi = FFTW.plan_ifft!(Wb, 2)
        if length(ws) < tid
            resize!(ws, tid)
        end
        ws[tid] = (Ws, Wb, pWsf, pWbi)
    end
    return ws[tid]
end

function integral1DXW(lb::FT, hb::FT, Nsample::IT, mod::Symbol) where {IT<:Integer,FT<:Real}
    X, W = zeros(FT, Nsample), zeros(FT, Nsample)

    if mod == :uni # Uniform integration
        dl = (hb - lb) / Nsample
        for j = 1:Nsample
            X[j] = lb + (j - 1) * dl + dl / 2
            W[j] = dl
        end
    elseif mod == :glq # Gauss-Legendre integration
        Dx = 0.5 * (hb - lb)
        center = 0.5 * (hb + lb)
        XGL, WGL = gausslegendre(Nsample)
        X .= center .+ Dx .* XGL
        W .= abs(Dx) .* WGL
    else
        error("Only :uni and :glq modes are supported")
    end
    return X, W
end

function octreeXWNCal(lb::FT, hb::FT, L::IT, mod::Symbol) where {IT<:Integer,FT<:Real}
    N = if mod == :uni
        2 * (L + 1)
    elseif mod == :glq
        L + 1
    else
        error("Only :uni and :glq modes are supported")
    end
    Xs, Ws = integral1DXW(lb, hb, N, mod)
    return Xs, Ws
end

function truncationLCal(cubel::FT; λ = 1.0) where {FT<:Real}
    rel_l = cubel / λ
    L = floor(Int, truncation_kernel(rel_l))
    return L
end

"""
    levelIntegralInfoCal(levelCubeEdgel; λ=1.0, L_min=0) -> (L, GLPolesInfo)

计算球面高斯求积（GL）层的截断项与角谱采样信息（论文 4.1.2 节）。

对截断项为 `L` 的层，球面多项式最高阶为 `p = 2L + 1`，需要 `L+1` 个
Gauss-Legendre 点覆盖 `θ ∈ [0, π]`（实际对 `cosθ` 求积）与 `2(L+1)` 个
均匀点覆盖 `φ ∈ [0, 2π]`，总采样点数 `N_p = 2(L+1)²`，权重为两方向权重之积
（论文式 (4-1)）：

```math
\\int f(\\hat{k})\\, d^2\\hat{k} \\approx \\sum_{p=1}^{N_p} W_p f(\\hat{k}_p), \\qquad
\\sum_p W_p = 4\\pi
```

# Arguments
- `levelCubeEdgel`: 层盒子边长（米）。
- `λ`: 局部波长（默认 1，配合以 λ 为单位输入的 `levelCubeEdgel`）。
- `L_min`: 截断项下限（用于小盒子/低频保护）。

# Returns
- `L`: 整数截断项（`ceil` 取整后与 `L_min` 取最大）。
- `GLPolesInfo`: `θ`/`φ` 采样点、权重与球面方向信息。
"""
function levelIntegralInfoCal(levelCubeEdgel::FT; λ = 1.0, L_min::Int = 0) where {FT<:Real}
    ## Calculate truncation number
    L = max(truncationLCal(levelCubeEdgel; λ = λ), L_min)

    ## Integration points and weights
    # Theta direction (Gauss-Legendre)
    Xcosθs, Wθs = octreeXWNCal(one(FT), -one(FT), L, :glq)
    Xθs = acos.(Xcosθs)

    # Phi direction (Uniform)
    Xϕs, Wϕs = octreeXWNCal(zero(FT), convert(FT, 2π), L, :uni)

    # Calculate all poles info
    r̂sθsϕs = [r̂θϕInfo(θ, ϕ) for ϕ in Xϕs for θ in Xθs]
    Wθϕs = [Wθ * Wϕ for Wϕ in Wϕs for Wθ in Wθs]

    Poles = GLPolesInfo{FT}(Xθs, Xϕs, Wθϕs, r̂sθsϕs)

    return L, Poles
end

"""
    interpolationCSCMatCal(pLevelPoles::GLPolesInfo{FT}, tLevelPoles::GLPolesInfo{FT}, nlocalInterp::IT)

Calculate sparse interpolation matrices from parent level to current level (or vice versa).
"""
function interpolationCSCMatCal(
    pLevelPoles::GLPolesInfo{FT},
    tLevelPoles::GLPolesInfo{FT},
    nlocalInterp::IT,
) where {IT<:Integer,FT<:Real}

    # Parent and Current level theta coordinates
    pXθs = pLevelPoles.Xθs
    tXθs = tLevelPoles.Xθs
    npXθs = length(pXθs)
    ntXθs = length(tXθs)

    # Parent and Current level phi coordinates
    pXϕs = pLevelPoles.Xϕs
    tXϕs = tLevelPoles.Xϕs
    npXϕs = length(pXϕs)
    ntXϕs = length(tXϕs)

    ################################################################
    # Theta direction

    nlocalInterpTheta = nlocalInterp
    if nlocalInterpTheta > npXθs
        nlocalInterpTheta = npXθs
    end

    pθsIntθs = cooraInCoorb(pXθs, tXθs)

    interWθs = ones(FT, (nlocalInterpTheta, npXθs))
    interIDθs = zeros(IT, (nlocalInterpTheta, npXθs))
    RelativeOffsets = (1:nlocalInterpTheta) .- nlocalInterpTheta ÷ 2

    @inbounds for ipXθs = 1:npXθs
        pθIntθ = pθsIntθs[ipXθs]
        iInterIDθs = pθIntθ .+ RelativeOffsets
        interIDθs[:, ipXθs] .= iInterIDθs
    end

    @inbounds for i = 1:nlocalInterpTheta
        θsInterp = [pickθ(idx, tXθs) for idx in interIDθs[i, :]]
        sinHalfDiffθpLevel = sin.((pXθs .- θsInterp) ./ 2)
        for j = 1:nlocalInterpTheta
            if i != j
                θsInterpLocal = [pickθ(idx, tXθs) for idx in interIDθs[j, :]]
                sinHalfDiffθtLevel = sin.((θsInterpLocal .- θsInterp) ./ 2)
                interWθs[j, :] .*= sinHalfDiffθpLevel ./ sinHalfDiffθtLevel
            end
        end
    end

    @inbounds for i = 1:npXθs
        interWθs[:, i] ./= sum(interWθs[:, i])
    end

    ################################################################
    # Phi direction

    nlocalInterpPhi = nlocalInterp
    if nlocalInterpPhi > npXϕs
        nlocalInterpPhi = npXϕs
    end

    pϕsIntϕs = cooraInCoorb(pXϕs, tXϕs)

    interWϕs = ones(FT, (nlocalInterpPhi, npXϕs))
    interIDϕs = zeros(IT, (nlocalInterpPhi, npXϕs))
    RelativeOffsets = (1:nlocalInterpPhi) .- nlocalInterpPhi ÷ 2

    @inbounds for ipXϕs = 1:npXϕs
        pϕIntϕ = pϕsIntϕs[ipXϕs]
        iInterIDϕs = collect(pϕIntϕ .+ RelativeOffsets)
        interIDϕs[:, ipXϕs] .= iInterIDϕs
    end

    @inbounds for i = 1:nlocalInterpPhi
        ϕsInterp = [pickϕ(idx, tXϕs) for idx in interIDϕs[i, :]]
        sinHalfDiffϕpLevel = sin.((pXϕs .- ϕsInterp) ./ 2)
        for j = 1:nlocalInterpPhi
            if i != j
                ϕsInterpLocal = [pickϕ(idx, tXϕs) for idx in interIDϕs[j, :]]
                sinHalfDiffϕtLevel = sin.((ϕsInterpLocal .- ϕsInterp) ./ 2)
                interWϕs[j, :] .*= sinHalfDiffϕpLevel ./ sinHalfDiffϕtLevel
            end
        end
    end

    @inbounds for i = 1:npXϕs
        interWϕs[:, i] ./= sum(view(interWϕs, :, i))
    end

    ################################################################
    # Construct Sparse Matrices (2-step interpolation)

    npSample = npXθs * npXϕs
    ntSample = ntXθs * ntXϕs
    ntempSample = ntXθs * npXϕs

    ### Step 1: Phi direction
    tSampleIndexes = reshape(collect(1:ntSample), ntXθs, ntXϕs)
    interIDGlobalϕs = repeat(interIDϕs, inner = (1, ntXθs))

    for itθ = 1:ntXθs
        for ipϕ = 1:npXϕs
            interIDGlobalϕs[:, (ipϕ-1)*ntXθs+itθ] .=
                [pickCycleVec(interIDϕ, tSampleIndexes[itθ, :]) for interIDϕ in interIDϕs[:, ipϕ]]
        end
    end

    rawIDϕs = repeat(collect(IT, 1:ntempSample); inner = nlocalInterpPhi)
    interWGlobalϕs = repeat(interWϕs, inner = (1, ntXθs))
    interpϕCSC = sparse(rawIDϕs, view(interIDGlobalϕs, :), view(interWGlobalϕs, :))


    #### Step 2: Theta direction
    for i in eachindex(interIDθs)
        ((interIDθs[i] < 1) | (interIDθs[i] > ntXθs)) && (interWθs[i] *= -1)
    end

    tempSampleIndexes = reshape(collect(1:ntempSample), ntXθs, npXϕs)
    interIDGlobalθs = repeat(interIDθs, outer = (1, npXϕs))

    halfnpϕ = npXϕs ÷ 2
    for ipϕ = 1:npXϕs
        for ipθ = 1:npXθs
            inGlobalIDs = zeros(IT, nlocalInterpTheta)
            for jInter = 1:nlocalInterpTheta
                interIDθ = interIDθs[jInter, ipθ]
                targetIdxInTempSampleIndexes = [interIDθ, ipϕ]

                if (interIDθ < 1) | (interIDθ > ntXθs)
                    if interIDθ < 1
                        targetIdxInTempSampleIndexes[1] = -interIDθ + 1
                    elseif interIDθ > ntXθs
                        targetIdxInTempSampleIndexes[1] = 2ntXθs + 1 - interIDθ
                    end

                    if ipϕ <= halfnpϕ
                        targetIdxInTempSampleIndexes[2] = ipϕ + halfnpϕ
                    else
                        targetIdxInTempSampleIndexes[2] = ipϕ - halfnpϕ
                    end
                end
                inGlobalIDs[jInter] = tempSampleIndexes[targetIdxInTempSampleIndexes...]
            end
            interIDGlobalθs[:, (ipϕ-1)*npXθs+ipθ] .= inGlobalIDs
        end
    end

    rawIDθs = repeat(collect(IT, 1:npSample); inner = nlocalInterpTheta)
    interWGlobalθs = repeat(interWθs; outer = (1, npXϕs))
    interpθCSC = sparse(rawIDθs, view(interIDGlobalθs, :), view(interWGlobalθs, :))

    dropzeros!(interpθCSC)
    dropzeros!(interpϕCSC)

    # Create struct
    info = LagrangeInterpInfo{IT,FT}()
    info.θCSC = interpθCSC
    info.ϕCSC = interpϕCSC
    info.θCSCT = sparse(transpose(interpθCSC))
    info.ϕCSCT = sparse(transpose(interpϕCSC))

    return info
end

function cooraInCoorb(coora::Vector{T}, coorb::Vector{T}) where {T<:Number}
    targetIDs = fill!(similar(coora, Int), 0)
    @inbounds for i in eachindex(coora)
        a = coora[i]
        for j in eachindex(coorb)
            b = coorb[j]
            if a >= b
                targetIDs[i] = j
                continue
            end
        end
    end
    return targetIDs
end

function pickϕ(index::Integer, ϕs::Vector{TT}) where {TT<:Real}
    re = zero(TT)
    ln = length(ϕs)
    if index < 1
        re = (ϕs[index+ln] - 2π)
    elseif index > ln
        re = ϕs[index-ln] + 2π
    else
        re = ϕs[index]
    end
    return re
end

function pickCycleVec(index::Integer, cycleVec::Vector{T}) where {T<:Real}
    ln = length(cycleVec)
    if index < 1
        re = cycleVec[index+ln]
    elseif index > ln
        re = cycleVec[index-ln]
    else
        re = cycleVec[index]
    end
    return re
end

function pickθ(index::Integer, θs::Vector{T}) where {T<:Real}
    re = zero(T)
    ln = length(θs)
    if index < 1
        re = -θs[-index+1]
    elseif index > ln
        re = 2pi - θs[2ln+1-index]
    else
        re = θs[index]
    end
    return re
end


end
