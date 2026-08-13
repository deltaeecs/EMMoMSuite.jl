"""
PMCHWMLFMAOperator.jl — Phase 15 步骤 15.8–15.11

PMCHW 系统的 MLFMA 线性算子（自包含模块）。

设计依据：Gibson《MoM》Ch.11 Algorithm 14 (双八叉树 + 四遍远场)

矩阵结构（2N×2N）：
  Z = [Z^EJ   Z^EM]
      [Z^HJ   Z^HM]
其中：
  Z^EJ = L(k0,η0) + L(k1,η1)    因子: jk η / (16π)
  Z^EM = K(k0)    + K(k1)        K 算子（PMCHW 专用，无 n̂× 测试）
  Z^HJ = -Z^EM                   精确结构不变量
  Z^HM = Lₑ(k0,1/η0) + Lₑ(k1,1/η1)  因子: jk / (η·16π)

实现：
  近场 Z_near — assemble_near_field_pmchw → 2N×2N 稀疏矩阵
              （原生装配：仅计算 octree 近邻 cube 对的相互作用，不装配全稠密矩阵）
  远场 Z_far  — 4 遍 MLFMA（J×k0, J×k1, M×k0, M×k1）

四遍数据流：
  遍 1: J 系数 (x[1:N])    × k0 八叉树 → 解聚写 y[1:N](EJ) + y[N+1:2N](HJ)
  遍 2: J 系数 (x[1:N])    × k1 八叉树 → 解聚同上（累加）
  遍 3: M 系数 (x[N+1:2N]) × k0 八叉树 → 解聚写 y[1:N](EM) + y[N+1:2N](HM)
  遍 4: M 系数 (x[N+1:2N]) × k1 八叉树 → 解聚同上（累加）
"""
module PMCHWMLFMAOperatorModule

using LinearAlgebra
using StaticArrays
using SparseArrays
using Logging
using MPI
using ....CoreModule
using ....Geometry
using ....BasisFunctions
using ....IntegralEquations
import ....Solvers: BlockJacobiPreconditioner
using ..Octree
using ..OctreeBuilder
using ..Translation
using ..Disaggregation: disaggregate_downward!
using ..Aggregation: aggregate_upward!
using ..Level
using ..Interpolation: AbstractPolesInfo

using ....IntegralEquations.Impedance: get_triangle_info, get_triangles_info
using ....IntegralEquations.EFIEModule: efie_interaction!
using ....IntegralEquations.PMCHWModule: PMCHW
import ....IntegralEquations.PMCHWModule: _l_block_operator, calc_k_pmchw_term!

export PMCHWMLFMAErrorBudget, PMCHWMLFMAOperator, PMCHWMLFMAOperatorMPI, assemble_near_field_pmchw

"""
    PMCHWMLFMAErrorBudget{FT}

PMCHW-MLFMA 的近场/远场误差预算参数：
- `leaf_wavelength_divisor`：叶子边长占波长的比例；
- `near_range_scale` / `min_near_range` / `max_near_range`：近邻范围缩放与上下限；
- `L_min`：最小截断阶数；
- `fixed_near_range` / `fixed_leaf_size_eff`：固定近邻范围 / 固定有效叶子尺寸（`0` 表示自动）。

构造：`PMCHWMLFMAErrorBudget(::Type{FT}=Float64; kwargs...)`。
"""
struct PMCHWMLFMAErrorBudget{FT<:AbstractFloat}
    leaf_wavelength_divisor::FT
    near_range_scale::FT
    min_near_range::Int
    max_near_range::Int
    L_min::Int
    fixed_near_range::Int
    fixed_leaf_size_eff::FT
end

function PMCHWMLFMAErrorBudget(
    ::Type{FT} = Float64;
    leaf_wavelength_divisor::Real = 10.0,
    near_range_scale::Real = 8.0,
    min_near_range::Int = 8,
    max_near_range::Int = 32,
    L_min::Int = 0,
    fixed_near_range::Int = 0,
    fixed_leaf_size_eff::Real = 0.0,
) where {FT<:AbstractFloat}
    return PMCHWMLFMAErrorBudget{FT}(
        FT(leaf_wavelength_divisor),
        FT(near_range_scale),
        min_near_range,
        max_near_range,
        L_min,
        fixed_near_range,
        FT(fixed_leaf_size_eff),
    )
end

function _resolve_budget_parameters(
    budget::PMCHWMLFMAErrorBudget{FT},
    λ_min::FT,
    leaf_size::FT,
) where {FT<:AbstractFloat}
    leaf_size_eff =
        iszero(budget.fixed_leaf_size_eff) ? λ_min / budget.leaf_wavelength_divisor :
        budget.fixed_leaf_size_eff
    near_range =
        budget.fixed_near_range > 0 ? budget.fixed_near_range :
        clamp(
            round(Int, leaf_size / leaf_size_eff * budget.near_range_scale),
            budget.min_near_range,
            budget.max_near_range,
        )
    return leaf_size_eff, near_range
end

# ─────────────────────────────────────────────────────────────────────────────
# Struct
# ─────────────────────────────────────────────────────────────────────────────

"""
    PMCHWMLFMAOperator{FT,CT} <: AbstractIntegralOperator

PMCHW 系统的 MLFMA 线性算子，大小 2N×2N。
"""
struct PMCHWMLFMAOperator{FT<:AbstractFloat,CT<:Complex} <: AbstractIntegralOperator
    pmchw::PMCHW{FT,CT}
    basis::RWGBasis
    Z_near::SparseMatrixCSC{CT,Int}
    budget::PMCHWMLFMAErrorBudget{FT}
    leaf_size_eff::FT
    near_range::Int
    octree0::OctreeInfo
    octree1::OctreeInfo
    sorted_ids0::Vector{Int}
    inv_sorted_ids0::Vector{Int}
    sorted_ids1::Vector{Int}
    inv_sorted_ids1::Vector{Int}
    freq::FT
end

Base.size(op::PMCHWMLFMAOperator) = (2 * num_basis(op.basis), 2 * num_basis(op.basis))
Base.size(op::PMCHWMLFMAOperator, d::Int) = 2 * num_basis(op.basis)
Base.eltype(::PMCHWMLFMAOperator{FT,CT}) where {FT,CT} = CT

function _leaf_block_indices(op::PMCHWMLFMAOperator)
    leaf_level = op.octree0.levels[op.octree0.nLevels]
    N = num_basis(op.basis)
    blocks = Vector{Vector{Int}}()

    for cube in leaf_level.cubes
        isempty(cube.bfInterval) && continue
        leaf_ids = collect(op.sorted_ids0[cube.bfInterval])
        n_leaf = length(leaf_ids)
        block = Vector{Int}(undef, 2 * n_leaf)
        copyto!(block, 1, leaf_ids, 1, n_leaf)
        @inbounds for i = 1:n_leaf
            block[n_leaf+i] = leaf_ids[i] + N
        end
        push!(blocks, block)
    end

    return blocks
end

BlockJacobiPreconditioner(op::PMCHWMLFMAOperator) =
    BlockJacobiPreconditioner(op.Z_near, _leaf_block_indices(op))
BlockJacobiPreconditioner(op::PMCHWMLFMAOperator, ::Any) = BlockJacobiPreconditioner(op)

# ─────────────────────────────────────────────────────────────────────────────
# 构造函数
# ─────────────────────────────────────────────────────────────────────────────

function PMCHWMLFMAOperator(
    pmchw::PMCHW,
    basis::RWGBasis,
    leaf_size::Float64;
    budget = PMCHWMLFMAErrorBudget(Float64),
)
    N = num_basis(basis)
    FT = typeof(real(pmchw.freq))
    CT = Complex{FT}

    centers = reduce(hcat, [bf.center for bf in basis.functions])  # 3 × N

    λ0 = FT(2π) / real(pmchw.k0)
    λ1 = FT(2π) / real(pmchw.k1)

    λ0f = Float64(λ0)
    λ1f = Float64(λ1)
    λ_min = FT(min(λ0f, λ1f))
    budget_ft =
        budget isa PMCHWMLFMAErrorBudget{FT} ? budget :
        PMCHWMLFMAErrorBudget(
            FT;
            leaf_wavelength_divisor = budget.leaf_wavelength_divisor,
            near_range_scale = budget.near_range_scale,
            min_near_range = budget.min_near_range,
            max_near_range = budget.max_near_range,
            L_min = budget.L_min,
            fixed_near_range = budget.fixed_near_range,
            fixed_leaf_size_eff = budget.fixed_leaf_size_eff,
        )
    leaf_size_eff, near_range = _resolve_budget_parameters(budget_ft, λ_min, FT(leaf_size))
    octree0, sorted_ids0 = build_octree(
        centers,
        leaf_size_eff;
        λ = λ0f,
        near_range = near_range,
        L_min = budget_ft.L_min,
    )
    octree1, sorted_ids1 = build_octree(
        centers,
        leaf_size_eff;
        λ = λ1f,
        near_range = near_range,
        L_min = budget_ft.L_min,
    )

    inv_sorted_ids0 = Vector{Int}(undef, N)
    inv_sorted_ids1 = Vector{Int}(undef, N)
    for i = 1:N
        inv_sorted_ids0[sorted_ids0[i]] = i
        inv_sorted_ids1[sorted_ids1[i]] = i
    end

    Z_near = assemble_near_field_pmchw(pmchw, basis, octree0, sorted_ids0, inv_sorted_ids0)

    return PMCHWMLFMAOperator{FT,CT}(
        pmchw,
        basis,
        Z_near,
        budget_ft,
        leaf_size_eff,
        near_range,
        octree0,
        octree1,
        sorted_ids0,
        inv_sorted_ids0,
        sorted_ids1,
        inv_sorted_ids1,
        FT(pmchw.freq),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# mul!
# ─────────────────────────────────────────────────────────────────────────────

function LinearAlgebra.mul!(y::AbstractVector, A::PMCHWMLFMAOperator, x::AbstractVector)
    fill!(y, zero(eltype(y)))
    N = num_basis(A.basis)

    # 近场
    mul!(y, A.Z_near, x)

    # 远场
    y_far = zeros(eltype(y), 2N)

    function clear_agg!(oct)
        for (_, lv) in oct.levels
            isdefined(lv, :aggS) && fill!(lv.aggS, zero(eltype(lv.aggS)))
            isdefined(lv, :disaggG) && fill!(lv.disaggG, zero(eltype(lv.disaggG)))
        end
    end

    function up_translate_down!(oct)
        for levelID = (oct.nLevels-1):-1:2
            aggregate_upward!(oct.levels[levelID], oct.levels[levelID+1])
        end
        for levelID = 2:oct.nLevels
            translate!(oct.levels[levelID])
        end
        for levelID = 2:(oct.nLevels-1)
            disaggregate_downward!(oct.levels[levelID], oct.levels[levelID+1])
        end
    end

    k0 = A.pmchw.k0
    k1 = A.pmchw.k1

    # 遍 1: J×k0
    clear_agg!(A.octree0)
    aggregate_leaf_pmchw!(A.octree0, A.basis, x, A.sorted_ids0, 1:N, k0)
    up_translate_down!(A.octree0)
    disaggregate_leaf_pmchw_j!(A.octree0, A.basis, A.pmchw, y_far, A.sorted_ids0, :k0)

    # 遍 2: J×k1
    clear_agg!(A.octree1)
    aggregate_leaf_pmchw!(A.octree1, A.basis, x, A.sorted_ids1, 1:N, k1)
    up_translate_down!(A.octree1)
    disaggregate_leaf_pmchw_j!(A.octree1, A.basis, A.pmchw, y_far, A.sorted_ids1, :k1)

    # 遍 3: M×k0
    clear_agg!(A.octree0)
    aggregate_leaf_pmchw!(A.octree0, A.basis, x, A.sorted_ids0, (N+1):(2N), k0)
    up_translate_down!(A.octree0)
    disaggregate_leaf_pmchw_m!(A.octree0, A.basis, A.pmchw, y_far, A.sorted_ids0, :k0)

    # 遍 4: M×k1
    clear_agg!(A.octree1)
    aggregate_leaf_pmchw!(A.octree1, A.basis, x, A.sorted_ids1, (N+1):(2N), k1)
    up_translate_down!(A.octree1)
    disaggregate_leaf_pmchw_m!(A.octree1, A.basis, A.pmchw, y_far, A.sorted_ids1, :k1)

    y .+= y_far
    return y
end

Base.:*(A::PMCHWMLFMAOperator, x::AbstractVector) = (y = similar(x); mul!(y, A, x); y)

# ─────────────────────────────────────────────────────────────────────────────
# 15.9: aggregate_leaf_pmchw!
# ─────────────────────────────────────────────────────────────────────────────

function aggregate_leaf_pmchw!(
    octree::OctreeInfo,
    basis::RWGBasis,
    x::AbstractVector,
    sorted_ids::Vector{Int},
    x_range::AbstractUnitRange{Int},
    k::Number,
    ;
    cube_filter = nothing,
)
    leaf_level = octree.levels[octree.nLevels]
    offset = first(x_range) - 1

    FT = eltype(leaf_level.cubeEdgel)
    CT = Complex{FT}
    JK = CT(im * k)

    poles = leaf_level.poles
    nPoles = length(poles.r̂sθsϕs)
    nCubes = leaf_level.nCubes

    if !isdefined(leaf_level, :aggS)
        leaf_level.aggS = zeros(CT, nPoles, 2, nCubes)
    else
        fill!(leaf_level.aggS, zero(CT))
    end

    tri_info = get_triangles_info(basis.mesh, basis)
    gq = GaussQuadratureInfo(:Triangle, 4, FT)
    n_qp = length(gq.weight)

    poles_r̂ = [p.r̂ for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]

    Threads.@threads for iCube = 1:nCubes
        cube_filter !== nothing && !cube_filter(iCube) && continue
        cube = leaf_level.cubes[iCube]
        cubeCenter = cube.center

        for bfID_sorted in cube.bfInterval
            bfID_orig = sorted_ids[bfID_sorted]
            coeff = x[bfID_orig+offset]
            abs(coeff) < 1e-12 && continue

            bf = basis.functions[bfID_orig]

            for i_supp = 1:2
                tri_idx = bf.support[i_supp]
                tri_idx == 0 && continue

                tri = tri_info[tri_idx]
                v_all = tri.vertices

                local_edge = bf.local_edge_idx[i_supp]
                v_opp = v_all[:, local_edge]
                sign_supp = bf.signs[i_supp]

                for i_qp = 1:n_qp
                    L = gq.coordinate[:, i_qp]
                    r = v_all * L
                    rho = r - v_opp
                    r_local = r - cubeCenter
                    factor_vec = sign_supp * bf.edge_length / 2 * gq.weight[i_qp] * coeff

                    for iPole = 1:nPoles
                        r̂ = poles_r̂[iPole]
                        phase = exp(JK * dot(r̂, r_local))
                        vec = rho * factor_vec * phase

                        leaf_level.aggS[iPole, 1, iCube] += dot(poles_θhat[iPole], vec)
                        leaf_level.aggS[iPole, 2, iCube] += dot(poles_ϕhat[iPole], vec)
                    end
                end
            end
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.10: _receive_terms
# ─────────────────────────────────────────────────────────────────────────────

function _receive_terms(
    bf::RWG,
    basis::RWGBasis,
    field::AbstractArray,
    r0::AbstractVector,
    k::Number,
    poles::AbstractPolesInfo,
)
    CT = typeof(complex(one(real(k))))
    FT = real(CT)
    te = zero(CT)
    tm = zero(CT)
    JK = CT(im * k)

    poles_r̂ = [p.r̂ for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]
    nPoles = length(poles_r̂)

    tri_info = get_triangles_info(basis.mesh, basis)
    gq = GaussQuadratureInfo(:Triangle, 4, FT)
    n_qp = length(gq.weight)

    for i_supp = 1:2
        tri_idx = bf.support[i_supp]
        tri_idx == 0 && continue

        tri = tri_info[tri_idx]
        v_all = tri.vertices

        local_edge = bf.local_edge_idx[i_supp]
        v_opp = v_all[:, local_edge]
        sign_supp = bf.signs[i_supp]

        for i_qp = 1:n_qp
            L = gq.coordinate[:, i_qp]
            r = v_all * L
            rho = r - v_opp
            r_local = r - r0
            w_f = sign_supp * bf.edge_length / 2 * gq.weight[i_qp]

            for iPole = 1:nPoles
                r̂ = poles_r̂[iPole]
                θhat = poles_θhat[iPole]
                ϕhat = poles_ϕhat[iPole]
                phase = exp(-JK * dot(r̂, r_local))

                Eθ = CT(field[iPole, 1])
                Eϕ = CT(field[iPole, 2])

                E_inc = (Eθ * θhat + Eϕ * ϕhat) * phase
                te += dot(rho, E_inc) * w_f

                rhat_cross_E = (Eθ * ϕhat - Eϕ * θhat) * phase
                tm += dot(rho, rhat_cross_E) * w_f
            end
        end
    end
    return te, tm
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.10: disaggregate_leaf_pmchw_j!
# ─────────────────────────────────────────────────────────────────────────────

function disaggregate_leaf_pmchw_j!(
    octree::OctreeInfo,
    basis::RWGBasis,
    pmchw::PMCHW,
    y::AbstractVector,
    sorted_ids::Vector{Int},
    kmode::Symbol,
    ;
    cube_filter = nothing,
)
    N = num_basis(basis)
    k, η = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)

    factor_EJ = im * k * η / (4π)
    factor_HJ = im * k / (4π)

    leaf_level = octree.levels[octree.nLevels]
    isdefined(leaf_level, :disaggG) || return

    Threads.@threads for iCube = 1:leaf_level.nCubes
        cube_filter !== nothing && !cube_filter(iCube) && continue
        cube = leaf_level.cubes[iCube]
        field = view(leaf_level.disaggG, :, :, iCube)
        r0 = cube.center

        for bfID_sorted in cube.bfInterval
            bfID = sorted_ids[bfID_sorted]
            bf = basis.functions[bfID]

            te, tm = _receive_terms(bf, basis, field, r0, k, leaf_level.poles)
            y[bfID] += te * factor_EJ
            y[bfID+N] += tm * factor_HJ
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.10: disaggregate_leaf_pmchw_m!
# ─────────────────────────────────────────────────────────────────────────────

function disaggregate_leaf_pmchw_m!(
    octree::OctreeInfo,
    basis::RWGBasis,
    pmchw::PMCHW,
    y::AbstractVector,
    sorted_ids::Vector{Int},
    kmode::Symbol,
    ;
    cube_filter = nothing,
)
    N = num_basis(basis)
    k, η = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)

    factor_EM = -im * k / (4π)
    factor_HM = im * k / (η * 4π)

    leaf_level = octree.levels[octree.nLevels]
    isdefined(leaf_level, :disaggG) || return

    Threads.@threads for iCube = 1:leaf_level.nCubes
        cube_filter !== nothing && !cube_filter(iCube) && continue
        cube = leaf_level.cubes[iCube]
        field = view(leaf_level.disaggG, :, :, iCube)
        r0 = cube.center

        for bfID_sorted in cube.bfInterval
            bfID = sorted_ids[bfID_sorted]
            bf = basis.functions[bfID]

            te, tm = _receive_terms(bf, basis, field, r0, k, leaf_level.poles)
            y[bfID] += tm * factor_EM
            y[bfID+N] += te * factor_HM
        end
    end
    return nothing
end

"""
    assemble_near_field_pmchw(pmchw, basis, octree, sorted_ids, inv_sorted_ids)

装配 PMCHW-MLFMA 的完整 `2N×2N` 近场稀疏矩阵：
先显式装配完整 PMCHW 矩阵，再按八叉树叶子近邻关系提取近场分块
（`EJ` / `EM` / `HJ` / `HM` 四块同序）。
"""
function assemble_near_field_pmchw(
    pmchw::PMCHW,
    basis::RWGBasis,
    octree::OctreeInfo,
    sorted_ids::Vector{Int},
    inv_sorted_ids::Vector{Int},
    ;
    cube_filter::Function = i -> true,
)
    N = num_basis(basis)
    CT = Complex{typeof(pmchw.k0)}
    FT = eltype(basis.mesh.node)
    leaf_level = octree.levels[octree.nLevels]
    n_cubes = leaf_level.nCubes

    # ── 1. 三角形 → RWG 基函数映射（与 MLFMA 近场装配一致） ─────────────
    nt = num_elements(basis.mesh)
    tri_to_rwg = [Vector{Tuple{Int,Int,Float64}}() for _ = 1:nt]
    for (i, f) in enumerate(basis.functions)
        for k = 1:2
            tri_id = f.support[k]
            if tri_id > 0
                push!(tri_to_rwg[tri_id], (f.local_edge_idx[k], i, f.signs[k]))
            end
        end
    end
    all_tris = get_triangles_info(basis.mesh, basis)

    # ── 2. cube → 唯一三角形 ─────────────────────────────────────────────
    cube_tris_vec = [Int[] for _ = 1:n_cubes]
    for i_cube = 1:n_cubes
        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        tris_set = Set{Int}()
        for sorted_idx in cube.bfInterval
            global_id = sorted_ids[sorted_idx]
            f = basis.functions[global_id]
            for k = 1:2
                f.support[k] > 0 && push!(tris_set, f.support[k])
            end
        end
        cube_tris_vec[i_cube] = collect(tris_set)
    end

    # ── 3. 子块算子与核（与稠密 PMCHW 装配完全相同的因子/积分核） ────────
    k0 = pmchw.k0
    eta0 = pmchw.eta0
    k0_c = CT(k0)
    eta0_c = CT(eta0)
    k1_c = pmchw.k1
    eta1_c = pmchw.eta1

    efie_ej0 = _l_block_operator(k0, eta0, k0_c, eta0_c, :EJ)
    efie_ej1 = _l_block_operator(k1_c, eta1_c, k1_c, eta1_c, :EJ)
    efie_hm0 = _l_block_operator(k0, eta0, k0_c, eta0_c, :HM)
    efie_hm1 = _l_block_operator(k1_c, eta1_c, k1_c, eta1_c, :HM)

    gq_far = GaussQuadratureInfo(:Triangle, 4, FT)
    gq_near = GaussQuadratureInfo(:Triangle, 7, FT)
    mfie_k0 = (k = k0, eta = one(FT), gq_far = gq_far, gq_near = gq_near)
    mfie_k1 = (k = k1_c, eta = one(CT), gq_far = gq_far, gq_near = gq_near)

    # K 核的高斯点（与 _assemble_K_pmchw_offdiag! 相同）
    N_points = length(gq_far.weight)
    quad_points = Vector{SVector{N_points,SVector{3,FT}}}(undef, nt)
    Threads.@threads for t = 1:nt
        v_idx = basis.mesh.triangles[:, t]
        v1 = SVector{3,FT}(basis.mesh.node[:, v_idx[1]])
        v2 = SVector{3,FT}(basis.mesh.node[:, v_idx[2]])
        v3 = SVector{3,FT}(basis.mesh.node[:, v_idx[3]])
        quad_points[t] = SVector{N_points,SVector{3,FT}}(
            v1 * gq_far.coordinate[1, i] +
            v2 * gq_far.coordinate[2, i] +
            v3 * gq_far.coordinate[3, i] for i = 1:N_points
        )
    end

    function k_pmchw_interaction!(Z_local, op, t_test, t_src)
        t_test.triID == t_src.triID && return nothing
        r_test = quad_points[t_test.triID]
        r_src = quad_points[t_src.triID]
        calc_k_pmchw_term!(Z_local, op, t_test, t_src, r_test, r_src)
        return nothing
    end

    # ── 4. 原生装配循环（仅 octree 近邻 cube 对；线程并行，按线程 COO） ──
    n_threads = Threads.nthreads()
    max_tid = Threads.maxthreadid()

    estimated_nnz_per_thread = 0
    non_empty_cubes = 0
    for i_c = 1:n_cubes
        cube_c = leaf_level.cubes[i_c]
        isempty(cube_c.bfInterval) && continue
        non_empty_cubes += 1
        n_my = length(cube_tris_vec[i_c])
        for ni in cube_c.neighbors
            n_neigh = length(cube_tris_vec[ni])
            estimated_nnz_per_thread += n_my * n_neigh * 9 * 4
        end
    end
    estimated_nnz_per_thread = max(1024, div(estimated_nnz_per_thread, max(n_threads, 1)))
    estimated_nnz_per_thread = min(estimated_nnz_per_thread, 10_000_000)

    Is = [Vector{Int}(undef, estimated_nnz_per_thread) for _ = 1:max_tid]
    Js = [Vector{Int}(undef, estimated_nnz_per_thread) for _ = 1:max_tid]
    Vs = [Vector{CT}(undef, estimated_nnz_per_thread) for _ = 1:max_tid]
    counts = zeros(Int, max_tid)

    counter = Threads.Atomic{Int}(0)
    total_cubes = n_cubes

    println("Assembling PMCHW Near Field Matrix (native, octree neighbor pairs) with $n_threads threads...")
    println("  Non-empty Cubes: $non_empty_cubes, est nnz/thread: $estimated_nnz_per_thread")

    Threads.@threads :static for i_cube = 1:n_cubes
        tid = Threads.threadid()
        cube_filter(i_cube) || continue

        c = Threads.atomic_add!(counter, 1)
        if c % 200 == 0 || c == total_cubes
            print("\rProgress: $c / $total_cubes cubes")
        end

        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        my_tris = cube_tris_vec[i_cube]

        # efie_interaction! 末尾会乘各算子的 factor，因此 k0/k1 必须用独立缓冲
        Z_ej0 = zeros(CT, 3, 3)
        Z_ej1 = zeros(CT, 3, 3)
        Z_ej = zeros(CT, 3, 3)
        Z_hm0 = zeros(CT, 3, 3)
        Z_hm1 = zeros(CT, 3, 3)
        Z_hm = zeros(CT, 3, 3)
        # calc_k_pmchw_term! 末尾同样原地乘边长因子，k0/k1 必须独立缓冲
        Z_em0 = zeros(CT, 3, 3)
        Z_em1 = zeros(CT, 3, 3)
        Z_em = zeros(CT, 3, 3)

        for neighbor_idx in cube.neighbors
            neighbor_cube = leaf_level.cubes[neighbor_idx]
            isempty(neighbor_cube.bfInterval) && continue
            neigh_tris = cube_tris_vec[neighbor_idx]

            for t_test in my_tris
                tri_test = all_tris[t_test]
                for t_src in neigh_tris
                    tri_src = all_tris[t_src]

                    fill!(Z_ej0, zero(CT))
                    fill!(Z_ej1, zero(CT))
                    fill!(Z_hm0, zero(CT))
                    fill!(Z_hm1, zero(CT))
                    fill!(Z_em0, zero(CT))
                    fill!(Z_em1, zero(CT))
                    efie_interaction!(Z_ej0, efie_ej0, tri_test, tri_src)
                    efie_interaction!(Z_ej1, efie_ej1, tri_test, tri_src)
                    efie_interaction!(Z_hm0, efie_hm0, tri_test, tri_src)
                    efie_interaction!(Z_hm1, efie_hm1, tri_test, tri_src)
                    @. Z_ej = Z_ej0 + Z_ej1
                    @. Z_hm = Z_hm0 + Z_hm1
                    k_pmchw_interaction!(Z_em0, mfie_k0, tri_test, tri_src)
                    k_pmchw_interaction!(Z_em1, mfie_k1, tri_test, tri_src)
                    @. Z_em = Z_em0 + Z_em1

                    distribute_term_pmchw!(
                        Is, Js, Vs, counts, tid,
                        Z_ej, Z_em, Z_hm,
                        tri_to_rwg[t_test], tri_to_rwg[t_src],
                        cube.bfInterval, neighbor_cube.bfInterval,
                        inv_sorted_ids, N,
                    )
                end
            end
        end
    end

    # ── 5. 合并线程缓冲 → 稀疏矩阵 ───────────────────────────────────────
    for tid = 1:max_tid
        ct = counts[tid]
        resize!(Is[tid], ct)
        resize!(Js[tid], ct)
        resize!(Vs[tid], ct)
    end
    I_total = _concat_thread_buffers_pmchw(Is, counts)
    J_total = _concat_thread_buffers_pmchw(Js, counts)
    V_total = _concat_thread_buffers_pmchw(Vs, counts)

    Z_near = sparse(I_total, J_total, V_total, 2N, 2N)
    # 滤除机器精度噪声：K 块中多三角形组合可抵消到 ~1e-16×scale，这些条目
    # 相对误差度量失真且对 matvec 无意义。按稀疏合并后的总和值做绝对阈值。
    if nnz(Z_near) > 0
        scale = maximum(abs, nonzeros(Z_near))
        tol = 1e-12 * abs(scale)
        if tol > 0
            I_f, J_f, V_f = findnz(Z_near)
            keep = abs.(V_f) .> tol
            if !all(keep)
                Z_near = sparse(I_f[keep], J_f[keep], V_f[keep], 2N, 2N)
            end
        end
    end
    @info "  [PMCHWMLFMAOperator] 近场矩阵（原生装配）：nnz=$(nnz(Z_near)) / $(2N)×$(2N)"
    return Z_near
end

@inline function _ensure_capacity_pmchw!(
    Is::Vector{Vector{Int}},
    Js::Vector{Vector{Int}},
    Vs::Vector{Vector{CT}},
    counts::Vector{Int},
    tid::Int,
    needed::Int,
) where {CT}
    required = counts[tid] + needed
    cap = length(Is[tid])
    if required > cap
        new_cap = max(required, cap * 2)
        resize!(Is[tid], new_cap)
        resize!(Js[tid], new_cap)
        resize!(Vs[tid], new_cap)
    end
    return nothing
end

@inline function _push4!(
    Is, Js, Vs, counts, tid,
    i1, j1, v1, i2, j2, v2, i3, j3, v3, i4, j4, v4,
)
    _ensure_capacity_pmchw!(Is, Js, Vs, counts, tid, 4)
    @inbounds begin
        ct = counts[tid]
        Is[tid][ct + 1] = i1; Js[tid][ct + 1] = j1; Vs[tid][ct + 1] = v1
        Is[tid][ct + 2] = i2; Js[tid][ct + 2] = j2; Vs[tid][ct + 2] = v2
        Is[tid][ct + 3] = i3; Js[tid][ct + 3] = j3; Vs[tid][ct + 3] = v3
        Is[tid][ct + 4] = i4; Js[tid][ct + 4] = j4; Vs[tid][ct + 4] = v4
        counts[tid] = ct + 4
    end
    return nothing
end

"""
    distribute_term_pmchw!(Is, Js, Vs, counts, tid, Z_ej, Z_em, Z_hm,
                           test_bases, src_bases, test_interval, src_interval,
                           inv_sorted_ids, N)

把 3×3 元素相互作用写入 PMCHW 四个子块：
  EJ: (row=test, col=src)
  EM: (row=test, col=src+N)
  HJ: (row=test+N, col=src) = -EM（结构不变量）
  HM: (row=test+N, col=src+N)
行/列均为原始（非排序）基函数编号；`test_interval`/`src_interval` 为排序后的
cube 区间，用于保证每个 cube 只写出其 own 的行（与 MLFMA 近场装配一致）。
"""
@inline function distribute_term_pmchw!(
    Is, Js, Vs, counts, tid,
    Z_ej, Z_em, Z_hm,
    test_bases,
    src_bases,
    test_interval,
    src_interval,
    inv_sorted_ids,
    N,
)
    _ensure_capacity_pmchw!(Is, Js, Vs, counts, tid, 4 * length(test_bases) * length(src_bases))
    @inbounds for (loc_test, glob_test, sign_test) in test_bases
        sorted_idx_test = inv_sorted_ids[glob_test]
        sorted_idx_test in test_interval || continue

        for (loc_src, glob_src, sign_src) in src_bases
            sorted_idx_src = inv_sorted_ids[glob_src]
            sorted_idx_src in src_interval || continue

            st = sign_test * sign_src
            v_ej = Z_ej[loc_test, loc_src] * st
            v_em = Z_em[loc_test, loc_src] * st
            v_hm = Z_hm[loc_test, loc_src] * st
            _push4!(
                Is, Js, Vs, counts, tid,
                glob_test, glob_src, v_ej,
                glob_test, glob_src + N, v_em,
                glob_test + N, glob_src, -v_em,
                glob_test + N, glob_src + N, v_hm,
            )
        end
    end
    return nothing
end

function _concat_thread_buffers_pmchw(buffers::Vector{Vector{T}}, counts::Vector{Int}) where {T}
    total_count = sum(counts)
    merged = Vector{T}(undef, total_count)
    offset = 0
    @inbounds for tid in eachindex(buffers)
        ct = counts[tid]
        if ct == 0
            continue
        end
        copyto!(merged, offset + 1, buffers[tid], 1, ct)
        offset += ct
    end
    return merged
end

"""
    PMCHWMLFMAOperatorMPI(pmchw, basis, leaf_size; budget=..., comm=MPI.COMM_WORLD)

PMCHW 的 MPI 混合算子：双八叉树（各秩一致），近场按 cube 分区原生装配
（rank 拥有 `(i_cube-1)%P` 号叶 cube 的 J/M 行，不装配全稠密矩阵），
远场四遍按 cube 分区（聚合/转移/反聚合 + 每层 Allreduce），秩内 `@threads`。
每遍使用独立 `y_pass` 缓冲，Allreduce 后累加，避免跨遍重复求和。
"""
struct PMCHWMLFMAOperatorMPI{FT,CT} <: AbstractIntegralOperator
    pmchw::PMCHW
    basis::RWGBasis
    Z_near_local::SparseMatrixCSC{CT,Int}
    budget::PMCHWMLFMAErrorBudget{FT}
    leaf_size_eff::FT
    near_range::Int
    rows::Vector{Int}
    octree0::OctreeInfo
    octree1::OctreeInfo
    sorted_ids0::Vector{Int}
    inv_sorted_ids0::Vector{Int}
    sorted_ids1::Vector{Int}
    inv_sorted_ids1::Vector{Int}
    freq::FT
    comm
end

function PMCHWMLFMAOperatorMPI(
    pmchw::PMCHW,
    basis::RWGBasis,
    leaf_size::Float64;
    budget = PMCHWMLFMAErrorBudget(Float64),
    comm = MPI.COMM_WORLD,
)
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    N = num_basis(basis)
    FT = typeof(real(pmchw.freq))
    CT = Complex{FT}

    centers = reduce(hcat, [bf.center for bf in basis.functions])
    λ0 = FT(2π) / real(pmchw.k0)
    λ1 = FT(2π) / real(pmchw.k1)
    λ0f = Float64(λ0)
    λ1f = Float64(λ1)
    λ_min = FT(min(λ0f, λ1f))
    budget_ft =
        budget isa PMCHWMLFMAErrorBudget{FT} ? budget :
        PMCHWMLFMAErrorBudget(
            FT;
            leaf_wavelength_divisor = budget.leaf_wavelength_divisor,
            near_range_scale = budget.near_range_scale,
            min_near_range = budget.min_near_range,
            max_near_range = budget.max_near_range,
            L_min = budget.L_min,
            fixed_near_range = budget.fixed_near_range,
            fixed_leaf_size_eff = budget.fixed_leaf_size_eff,
        )
    leaf_size_eff, near_range = _resolve_budget_parameters(budget_ft, λ_min, FT(leaf_size))
    octree0, sorted_ids0 = build_octree(
        centers, leaf_size_eff;
        λ = λ0f, near_range = near_range, L_min = budget_ft.L_min,
    )
    octree1, sorted_ids1 = build_octree(
        centers, leaf_size_eff;
        λ = λ1f, near_range = near_range, L_min = budget_ft.L_min,
    )
    inv_sorted_ids0 = Vector{Int}(undef, N)
    inv_sorted_ids1 = Vector{Int}(undef, N)
    for i = 1:N
        inv_sorted_ids0[sorted_ids0[i]] = i
        inv_sorted_ids1[sorted_ids1[i]] = i
    end

    # 近场：按叶 cube 分区原生装配（仅计算近邻对，无全稠密矩阵）。
    rank_filter = i_cube -> (i_cube - 1) % P == rank
    Z_near_local = assemble_near_field_pmchw(
        pmchw, basis, octree0, sorted_ids0, inv_sorted_ids0;
        cube_filter = rank_filter,
    )
    rows = Int[]
    leaf_level = octree0.levels[octree0.nLevels]
    for (i_cube, cube) in enumerate(leaf_level.cubes)
        (i_cube - 1) % P == rank || continue
        for s in cube.bfInterval
            push!(rows, sorted_ids0[s])
            push!(rows, sorted_ids0[s] + N)
        end
    end
    sort!(rows)

    return PMCHWMLFMAOperatorMPI{FT,CT}(
        pmchw,
        basis,
        Z_near_local,
        budget_ft,
        leaf_size_eff,
        near_range,
        rows,
        octree0,
        octree1,
        sorted_ids0,
        inv_sorted_ids0,
        sorted_ids1,
        inv_sorted_ids1,
        FT(pmchw.freq),
        comm,
    )
end

Base.size(op::PMCHWMLFMAOperatorMPI) = (2 * num_basis(op.basis), 2 * num_basis(op.basis))
Base.size(op::PMCHWMLFMAOperatorMPI, d::Int) = 2 * num_basis(op.basis)
Base.eltype(::PMCHWMLFMAOperatorMPI{FT,CT}) where {FT,CT} = CT

function _pmchw_pass_mpi!(
    oct,
    A::PMCHWMLFMAOperatorMPI,
    sorted_ids,
    x,
    x_range,
    k,
    kmode,
    kind,
    cube_filter,
    y_pass,
    comm,
)
    for (_, lv) in oct.levels
        isdefined(lv, :aggS) && fill!(lv.aggS, zero(eltype(lv.aggS)))
        isdefined(lv, :disaggG) && fill!(lv.disaggG, zero(eltype(lv.disaggG)))
    end
    nL = oct.nLevels
    aggregate_leaf_pmchw!(oct, A.basis, x, sorted_ids, x_range, k; cube_filter = cube_filter)
    MPI.Allreduce!(oct.levels[nL].aggS, +, comm)
    for lv in (nL - 1):-1:2
        aggregate_upward!(oct.levels[lv], oct.levels[lv + 1]; cube_filter = cube_filter)
        MPI.Allreduce!(oct.levels[lv].aggS, +, comm)
    end
    for lv in 2:nL
        translate!(oct.levels[lv]; cube_filter = cube_filter)
    end
    if nL >= 2
        MPI.Allreduce!(oct.levels[2].disaggG, +, comm)
    end
    for lv in 2:(nL - 1)
        disaggregate_downward!(oct.levels[lv], oct.levels[lv + 1]; child_filter = cube_filter)
        if lv + 1 < nL
            MPI.Allreduce!(oct.levels[lv + 1].disaggG, +, comm)
        end
    end
    if kind === :j
        disaggregate_leaf_pmchw_j!(
            oct, A.basis, A.pmchw, y_pass, sorted_ids, kmode; cube_filter = cube_filter
        )
    else
        disaggregate_leaf_pmchw_m!(
            oct, A.basis, A.pmchw, y_pass, sorted_ids, kmode; cube_filter = cube_filter
        )
    end
    MPI.Allreduce!(y_pass, +, comm)
    return nothing
end

function LinearAlgebra.mul!(y::AbstractVector, A::PMCHWMLFMAOperatorMPI, x::AbstractVector)
    N = num_basis(A.basis)
    M = 2N
    rank = MPI.Comm_rank(A.comm)
    n_procs = MPI.Comm_size(A.comm)
    cube_filter = n_procs > 1 ? (i -> (i - 1) % n_procs == rank) : nothing

    # 近场（本地 cube 行）+ Allreduce：Z_near_local 只含本秩拥有的行，
    # 乘完后 y 仅这些行非零，Allreduce 得到完整近场乘积。
    fill!(y, 0)
    mul!(y, A.Z_near_local, x)
    MPI.Allreduce!(y, +, A.comm)

    # 远场：四遍，每遍独立 y_pass（Allreduce 后累加，避免跨遍重复求和）
    y_far = zeros(ComplexF64, M)
    y_pass = zeros(ComplexF64, M)
    k0 = A.pmchw.k0
    k1 = A.pmchw.k1
    _pmchw_pass_mpi!(A.octree0, A, A.sorted_ids0, x, 1:N, k0, :k0, :j, cube_filter, y_pass, A.comm)
    y_far .+= y_pass
    fill!(y_pass, 0)
    _pmchw_pass_mpi!(A.octree1, A, A.sorted_ids1, x, 1:N, k1, :k1, :j, cube_filter, y_pass, A.comm)
    y_far .+= y_pass
    fill!(y_pass, 0)
    _pmchw_pass_mpi!(A.octree0, A, A.sorted_ids0, x, (N+1):(2N), k0, :k0, :m, cube_filter, y_pass, A.comm)
    y_far .+= y_pass
    fill!(y_pass, 0)
    _pmchw_pass_mpi!(A.octree1, A, A.sorted_ids1, x, (N+1):(2N), k1, :k1, :m, cube_filter, y_pass, A.comm)
    y_far .+= y_pass

    y .+= y_far
    return y
end

Base.:*(A::PMCHWMLFMAOperatorMPI, x::AbstractVector) = (y = similar(x); mul!(y, A, x); y)

end # module PMCHWMLFMAOperatorModule
