module Precomputations

using LinearAlgebra
using StaticArrays
using SparseArrays
using SpecialFunctions
using ....CoreModule
using ..Level
using ..Interpolation

export compute_shift_factors!, compute_interpolation_matrices!

"""
    compute_shift_factors!(nLevels::Int, levels::Dict{Int, LevelInfo}, k::Real)

Compute phase shift factors for all levels.
"""
function compute_shift_factors!(
    nLevels::Int,
    levels::Dict{Int,LV},
    k::Real,
) where {LV<:AbstractLevel}
    JK = im * k

    FT = eltype(levels[nLevels].cubeEdgel)

    # From level 2 to nLevels-1
    for iLevel = 2:(nLevels-1)
        level = levels[iLevel]

        # Offsets for 4 children (z < 0)
        # The other 4 are symmetric
        kidsOffsets = SMatrix{3,4,FT}(-1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1)

        # Vectors from parent center to child centers
        # cubeEdgel/4 * offset
        ΔCt2Cks = level.cubeEdgel / 4 .* kidsOffsets

        poles = level.poles
        nPoles = length(poles.r̂sθsϕs)

        # Allocate
        phaseShift2Kids = zeros(Complex{FT}, nPoles, 8)

        for iKid = 1:4
            ΔCt2Ck = ΔCt2Cks[:, iKid]
            for iPole = 1:nPoles
                r̂ = poles.r̂sθsϕs[iPole].r̂
                phaseShift2Kids[iPole, iKid] = exp(-JK * dot(r̂, ΔCt2Ck))
            end
        end

        # Symmetric children (conjugate)
        phaseShift2Kids[:, 8:-1:5] .= conj.(phaseShift2Kids[:, 1:4])

        # From kids to parent (conjugate)
        phaseShiftFromKids = conj.(phaseShift2Kids)

        level.phaseShift2Kids = phaseShift2Kids
        level.phaseShiftFromKids = phaseShiftFromKids
    end
end

"""
    compute_interpolation_matrices!(nLevels::Int, levels::Dict{Int, LevelInfo})

Compute interpolation matrices for all levels.
"""
function compute_interpolation_matrices!(
    nLevels::Int,
    levels::Dict{Int,LV},
) where {LV<:AbstractLevel}
    nInterp = 6 # Default interpolation points

    Threads.@threads for iLevel = nLevels:-1:2
        parentLevel = levels[iLevel-1]
        childLevel = levels[iLevel]

        # Compute interpolation matrix between Parent and Child.
        # FFTSpectral 路径：θ 方向 Lagrange 矩阵复用两步法构造，φ 方向改由 FFT 谱插值
        # （插值/反插值在 Aggregation/Disaggregation 调用 fft_interp_phi/fft_anterp_phi）。
        childLevel.interpWθϕ = if childLevel.poles isa FFTGLPolesInfo
            info = interpolationCSCMatCal(parentLevel.poles.inner, childLevel.poles.inner, nInterp)
            FFTInterpInfo(
                info.θCSC,
                info.θCSCT,
                length(childLevel.poles.Xθs),
                length(childLevel.poles.Xϕs),
                length(parentLevel.poles.Xϕs),
                childLevel.nCubes,
            )
        else
            interpolationCSCMatCal(parentLevel.poles, childLevel.poles, nInterp)
        end
    end
end

end
