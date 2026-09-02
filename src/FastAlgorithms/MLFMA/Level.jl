module Level

using StaticArrays
using OffsetArrays
using ..Interpolation

export CubeInfo, AbstractLevel, LevelInfo, adaptive_near_range, scaled_translation_factor

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

    # ── 稳态 matvec 分配缓存（issue #22 第二轮：全部懒初始化，undef 合法）──
    # 叶层辐射/接收循环用的极点向量数组（构造一次）
    polevecs::Any
    # aggregate_upward! 的逐 chunk scratch / FFT 批量缓冲（含创建时的线程数，
    # 线程数变化时自动重建）
    uwScratch::Any
    # disaggregate_downward! 的逐 chunk scratch / FFT 批量缓冲
    dwScratch::Any
    # PMCHW 路径的 (tri_info, gq) 单元缓存（原每遍重建全网格几何）
    elemcache::Any

    function LevelInfo{IT,FT,IPT}() where {IT<:Integer,FT<:Real,IPT}
        new{IT,FT,IPT}()
    end
end

"""
    adaptive_near_range(λ, cubeEdgel, L; kr_factor = 0.55) -> Int

Tree-uniform adaptive near-field radius derived from the leaf level
(issue #22, problem 3).

The real-spectrum M2L series `Σ (2l+1)(-j)^l h_l⁽²⁾(kR) P_l(k̂·R̂)` contains
`h_l⁽²⁾(kR)` terms that grow with `l` until `l ≈ kR`; the series only converges
in floating point once the smallest far-pair distance satisfies
`kR_min ≳ kr_factor · L`. Empirically calibrated against gate GD2V
(λ = 1, w = 0.1 → L = 9): `kr_factor = 0.55` reproduces the passing
`near_range = 7`, while `near_range = 4` (kR_min/L ≈ 0.35) fails.

The near/far list tiling of the octree requires ONE radius for the whole tree
(the child-level exclusion window must fall inside the parent-level neighbor
window), so the radius is derived from the **leaf** level — the smallest `w`,
hence the binding constraint — and applied uniformly:

    near_range = max(1, ceil(kr_factor · L_leaf / (k·w_leaf)) − 1),   k = 2π/λ

Coarser levels then satisfy the criterion with a growing margin
(kR_min(ℓ) scales with `w_ℓ` while `L_ℓ` grows sublinearly).
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

"""
    scaled_translation_factor(k, cubeEdgel, near_range) -> Real

Per-level scaling factor `s` for the low-frequency-stable (scaled) M2L
diagonalization — Ergül & Karaosmanoğlu, URSI GA 2014 (issue #22, problem 3).

The scaled decomposition replaces the unstable Hankel series terms
`h_l⁽²⁾(kw)` by `s^l · h_l⁽²⁾(kw)` in the translation operator and rescales
the shift operators to `exp(ik k̂·v/s)`. Stability requires `s ≪ kw` (smallest
translation distance), while the exponential shift approximation requires
`s ≫ kv` (largest shift distance). The geometric mean of both bounds gives

    s = k·w·sqrt(2·near_range + 1)

with `s` proportional to the electrical box size, as recommended by the paper.
"""
function scaled_translation_factor(
    k::Real,
    cubeEdgel::Real,
    near_range::Integer,
)
    k > 0 || return 1.0
    return k * cubeEdgel * sqrt(2 * near_range + 1)
end

end
