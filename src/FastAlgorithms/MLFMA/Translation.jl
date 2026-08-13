module Translation

using LinearAlgebra
using StaticArrays
using OffsetArrays
using SpecialFunctions
using ....CoreModule
using ..Level

export compute_translation_factors!, translate!

"""
    compute_translation_factors!(nLevels::Int, levels::Dict{Int, LevelInfo}, k::Real)

预计算各层远亲盒子间的转移函数（论文式 (2-41)）。

```math
T_\\tau(k, \\hat{k}, R_{ba}) = \\frac{-{\\rm j}k}{(4\\pi)^2}
\\sum_{l=0}^{\\tau} (-{\\rm j})^l (2l+1)\\, h_l^{(2)}(kR_{ba})\\, P_l(\\hat{k} \\cdot \\hat{R}_{ba})
```

实现中对每个远亲相对坐标 `Δ = (Δx, Δy, Δz)`（`|Δ|∞ > near_range`）计算
`R_{ba} = cubeEdgel * Δ`，再按采样点累加球汉克尔函数与勒让德多项式，
最后乘 `-jk/(4π)` 与求积权重 `W_p`（`αTrans` 即含权重的 `W_p T_τ`，
与论文式 (2-43) 的球面积分加权形式一致）。

# Arguments
- `nLevels`: 八叉树层数（从第 2 层到叶层计算）。
- `levels`: 层信息字典，每个 `LevelInfo` 写入 `αTrans` 与 `αTransIndex`。
- `k`: 波数。
- `near_range`: 近邻判定半径（默认 4，远邻条件为任一维度 `|Δ| > near_range`）。
"""
function compute_translation_factors!(
    nLevels::Int,
    levels::Dict{Int,LV},
    k::Real;
    near_range::Int = 4,
) where {LV<:AbstractLevel}
    # Compute for each level (from 2 to nLevels)
    for iLevel = 2:nLevels
        level = levels[iLevel]
        cal_alpha_trans_on_level!(level, k, near_range)
    end
end

function cal_alpha_trans_on_level!(
    level::LevelInfo,
    k::Real,
    near_range::Int,
)
    FT = eltype(level.cubeEdgel)
    CT = Complex{FT}

    truncL = level.L
    poles = level.poles
    nPoles = length(poles.r̂sθsϕs)
    Wθϕs = poles.Wθϕs

    # Constants
    # Standard code uses -im * k / (4 * π)
    # Verification shows we need this standard factor to match Direct Solver.
    const_factor = -im * k / (4 * FT(π))

    # 只计算/存储本层实际出现的远场相对偏移列：
    # 全量偏移表为 (2·(2·near_range+1)+1)³ − (2·near_range+1)³ 列，稀疏八叉树中
    # 绝大多数列从未被任何 cube 的 farneighbors 引用（实测 N=594、near_range=16 时
    # 仅 1.7% 被使用，低层远邻集为空却仍分配数 GB）。改为按使用集压缩。
    max_range = 2 * near_range + 1
    seen = Set{Tuple{Int,Int,Int}}()
    used = Tuple{Int,Int,Int}[]
    for cube in level.cubes
        for iFarNei in cube.farneighbors
            far = level.cubes[iFarNei]
            off = (
                cube.ID3D[1] - far.ID3D[1],
                cube.ID3D[2] - far.ID3D[2],
                cube.ID3D[3] - far.ID3D[3],
            )
            if !(off in seen)
                push!(seen, off)
                push!(used, off)
            end
        end
    end
    n_used = length(used)

    αTransIndex = OffsetArray(
        zeros(Int, 2 * max_range + 1, 2 * max_range + 1, 2 * max_range + 1),
        -max_range:max_range,
        -max_range:max_range,
        -max_range:max_range,
    )
    for (i, off) in enumerate(used)
        αTransIndex[off[1], off[2], off[3]] = i
    end
    level.αTransIndex = αTransIndex

    αTrans = zeros(CT, nPoles, n_used)

    # Loop over used far directions
    Threads.@threads for iFarNei = 1:n_used
        off = used[iFarNei]
        # Relative vector * cube edge length
        RabVec = SVector{3,FT}(
            off[1] * level.cubeEdgel,
            off[2] * level.cubeEdgel,
            off[3] * level.cubeEdgel,
        )
        Rab = norm(RabVec)

        R̂ab = RabVec / Rab

        # Spherical Hankel H2 (0 to truncL)
        # argument x = k * Rab
        x = k * Rab
        h2lxs = spherical_h2l_array(truncL, x)

        for iPole = 1:nPoles
            r̂ = poles.r̂sθsϕs[iPole].r̂
            cosϕ = clamp(dot(r̂, R̂ab), -one(FT), one(FT))

            # Legendre Polynomials P_l(cosϕ)
            legendrePls = collectPl(truncL, cosϕ)

            # Summation
            val = zero(CT)
            j_term = im
            for l = 0:truncL
                j_term *= -im
                val += j_term * (2 * l + 1) * h2lxs[l+1] * legendrePls[l+1]
            end

            αTrans[iPole, iFarNei] = val * const_factor * Wθϕs[iPole]
        end
    end

    level.αTrans = αTrans
    return nothing
end

function spherical_h2l_array(lmax::Int, x::T) where {T<:Real}
    # h_l^(2)(x) = j_l(x) - i y_l(x)
    res = Vector{Complex{T}}(undef, lmax + 1)
    for l = 0:lmax
        res[l+1] = sphericalbesselj(l, x) - im * sphericalbessely(l, x)
    end
    return res
end

function collectPl(lmax::Int, x::T) where {T<:Real}
    # Compute Legendre polynomials P_l(x) for l=0 to lmax
    res = Vector{T}(undef, lmax + 1)
    res[1] = one(T) # P_0
    if lmax >= 1
        res[2] = x # P_1
    end

    for l = 1:(lmax-1)
        # P_{l+1} = ((2l+1)x P_l - l P_{l-1}) / (l+1)
        res[l+2] = ((2 * l + 1) * x * res[l+1] - l * res[l]) / (l + 1)
    end
    return res
end

"""
    translate!(level::LevelInfo)

Perform translation (Horizontal Pass) for a single level.
"""
function translate!(level::LevelInfo; cube_filter = nothing)
    cubes = level.cubes
    aggS = level.aggS

    # Ensure disaggG is allocated and zeroed
    if !isdefined(level, :disaggG)
        level.disaggG = zeros(Complex{eltype(level.cubeEdgel)}, size(aggS))
    else
        fill!(level.disaggG, 0)
    end
    disaggG = level.disaggG

    αTrans = level.αTrans
    αTransIndex = level.αTransIndex

    # Loop over cubes
    count_trans = 0
    max_factor = 0.0

    Threads.@threads for iCube = 1:length(cubes)
        cube_filter !== nothing && !cube_filter(iCube) && continue
        cube = cubes[iCube]
        farNeighborIDs = cube.farneighbors

        for iFarNei in farNeighborIDs
            farNeiCube = cubes[iFarNei]

            # Relative 3D ID: Target - Source
            relative3DID = cube.ID3D .- farNeiCube.ID3D

            # Get index in αTrans array
            idx = αTransIndex[relative3DID[1], relative3DID[2], relative3DID[3]]

            # Apply translation
            factor = view(αTrans, :, idx)

            src = view(aggS, :, :, iFarNei)
            dest = view(disaggG, :, :, iCube)

            for pol = 1:2
                @views dest[:, pol] .+= factor .* src[:, pol]
            end
        end
    end

end

end
