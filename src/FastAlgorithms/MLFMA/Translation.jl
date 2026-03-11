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

Compute translation factors (alpha) for all levels.
"""
function compute_translation_factors!(
    nLevels::Int,
    levels::Dict{Int,LV},
    k::Real;
    near_range::Int = 4,
) where {LV<:AbstractLevel}

    # Precompute the far neighbor relative indices.
    # Theory: parent near_range = N → child far offset ≤ 2N+1 (proven by child-ID arithmetic).
    # Far condition: at least one dimension |Δ| > near_range (matches searchNearCubes).
    max_range = 2 * near_range + 1
    n_max = (2 * max_range + 1)^3

    all316FarNeighID = zeros(Int, 3, n_max)
    all343InFar316 = OffsetArray(
        zeros(Int, 2 * max_range + 1, 2 * max_range + 1, 2 * max_range + 1),
        -max_range:max_range,
        -max_range:max_range,
        -max_range:max_range,
    )

    indexfar316 = 0
    for kz = -max_range:max_range
        for ky = -max_range:max_range
            for kx = -max_range:max_range
                if (abs(kx) > near_range) || (abs(ky) > near_range) || (abs(kz) > near_range)
                    indexfar316 += 1
                    all343InFar316[kx, ky, kz] = indexfar316
                    all316FarNeighID[:, indexfar316] .= [kx, ky, kz]
                end
            end
        end
    end

    # Compute for each level (from 2 to nLevels)
    for iLevel = 2:nLevels
        level = levels[iLevel]
        cal_alpha_trans_on_level!(level, all316FarNeighID, all343InFar316, k, indexfar316)
    end
end

function cal_alpha_trans_on_level!(
    level::LevelInfo,
    all316FarNeighID::Matrix{Int},
    all343InFar316::OffsetArray,
    k::Real,
    nFar::Int,
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

    αTrans = zeros(CT, nPoles, nFar)

    # Loop over nFar directions
    Threads.@threads for iFarNei = 1:nFar
        # Relative vector * cube edge length
        RabVec = SVector{3,FT}(level.cubeEdgel .* all316FarNeighID[:, iFarNei])
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
    level.αTransIndex = all343InFar316
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
function translate!(level::LevelInfo)
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
        cube = cubes[iCube]
        farNeighborIDs = cube.farneighbors

        # Atomic add or just ignore race condition for debug
        # count_trans += length(farNeighborIDs)

        for iFarNei in farNeighborIDs
            farNeiCube = cubes[iFarNei]

            # Relative 3D ID: Target - Source
            relative3DID = cube.ID3D .- farNeiCube.ID3D

            # Get index in αTrans array
            idx = αTransIndex[relative3DID[1], relative3DID[2], relative3DID[3]]

            # Apply translation
            factor = view(αTrans, :, idx)

            # Debug magnitude
            # fmag = maximum(abs.(factor))
            # if fmag > max_factor; max_factor = fmag; end

            src = view(aggS, :, :, iFarNei)
            dest = view(disaggG, :, :, iCube)

            for pol = 1:2
                @views dest[:, pol] .+= factor .* src[:, pol]
            end
        end
    end

    # println("Debug: Level ", level.ID, " Translation Max Factor: ", maximum(abs.(αTrans)))
    # println("Debug: Level ", level.ID, " Total Interactions: ", count_trans)
end

end
