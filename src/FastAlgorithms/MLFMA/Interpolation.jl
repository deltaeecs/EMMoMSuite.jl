module Interpolation

using StaticArrays
using SparseArrays
using FastGaussQuadrature
using LinearAlgebra
using ....Geometry

export AbstractPolesInfo, AbstractInterpInfo
export GLPolesInfo, LagrangeInterpInfo
export levelIntegralInfoCal, interpolationCSCMatCal, truncation_kernel
export interpolate, anterpolate, interpolate!, anterpolate!
export interp_type

abstract type AbstractPolesInfo{FT<:AbstractFloat} end
abstract type AbstractInterpInfo{IT<:Integer,FT<:Real} end

function interpolate end
function anterpolate end
function interpolate! end
function anterpolate! end

const NBDIGITS = 9.0 # Default value (d0 in the truncation formula)

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
function truncation_kernel(rel_l, nbdigits::Real = NBDIGITS)
    return 2π * rel_l * sqrt(3) + 2.16 * nbdigits^(2.0 / 3.0) * (2π * rel_l)^(1 / 3)
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

function truncationLCal(cubel::FT; λ = 1.0, precision_digits::Real = NBDIGITS) where {FT<:Real}
    rel_l = cubel / λ
    L = floor(Int, truncation_kernel(rel_l, precision_digits))
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
function levelIntegralInfoCal(levelCubeEdgel::FT; λ = 1.0, L_min::Int = 0,
                              precision_digits::Real = NBDIGITS) where {FT<:Real}
    ## Calculate truncation number
    L = max(truncationLCal(levelCubeEdgel; λ = λ, precision_digits = precision_digits), L_min)

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

    rawIDϕs = repeat(collect(IT, 1:ntempSample); inner = nlocalInterp)
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
            inGlobalIDs = zeros(IT, nlocalInterp)
            for jInter = 1:nlocalInterp
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

    rawIDθs = repeat(collect(IT, 1:npSample); inner = nlocalInterp)
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
