module Level

using StaticArrays
using OffsetArrays
using ..Interpolation

export CubeInfo, AbstractLevel, LevelInfo, adaptive_near_range

mutable struct CubeInfo{IT<:Integer,FT<:Real}
    kidsInterval::UnitRange{IT}
    bfInterval::UnitRange{IT}
    kidsIn8::Vector{IT}
    geoIDs::Vector{IT}
    neighbors::Vector{IT}
    farneighbors::Vector{IT}
    ID3D::MVector{3,IT}
    center::MVector{3,FT}
end

abstract type AbstractLevel end

mutable struct LevelInfo{IT<:Integer,FT<:Real,IPT} <: AbstractLevel
    ID::IT
    isleaf::Bool
    L::IT
    nCubes::IT
    cubes::Vector{CubeInfo{IT,FT}}
    cubeEdgel::FT
    poles::AbstractPolesInfo{FT}
    interpWθϕ::IPT
    aggS::Array{Complex{FT},3}
    disaggG::Array{Complex{FT},3}
    phaseShift2Kids::Array{Complex{FT},2}
    phaseShiftFromKids::Array{Complex{FT},2}
    αTrans::Array{Complex{FT},2}
    αTransIndex::OffsetArray{IT,3,Array{IT,3}}

    function LevelInfo{IT,FT,IPT}() where {IT<:Integer,FT<:Real,IPT}
        new{IT,FT,IPT}()
    end
end

"""
    adaptive_near_range(λ, cubeEdgel, L; kr_factor = 0.55) -> Int

Per-level adaptive near-field radius (issue #22, problem 3).

The real-spectrum M2L series `Σ (2l+1)(-j)^l h_l⁽²⁾(kR) P_l(k̂·R̂)` contains
`h_l⁽²⁾(kR)` terms that grow with `l` until `l ≈ kR`; the series only converges
in floating point once the smallest far-pair distance satisfies
`kR_min ≳ kr_factor · L`. Empirically calibrated against gate GD2V
(λ = 1, w = 0.1 → L = 9): `kr_factor = 0.55` reproduces the passing
`near_range = 7`, while `near_range = 4` (kR_min/L ≈ 0.35) fails.

Since the smallest far-pair offset equals `near_range + 1`, the per-level
radius is

    near_range = max(1, ceil(kr_factor · L / (k·w)) − 1),   k = 2π/λ

The leaf level (smallest `w`, hence the strictest constraint) is binding.
"""
function adaptive_near_range(
    λ::Real,
    cubeEdgel::Real,
    L::Integer;
    kr_factor::Real = 0.55,
)
    k = 2π / λ
    kw = k * cubeEdgel
    kw > 0 || return 1
    return max(1, ceil(Int, kr_factor * L / kw) - 1)
end

end
