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
              （策略：计算全矩阵后依据 octree 近邻关系稀疏化）
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
using ....CoreModule
using ....Geometry
using ....BasisFunctions
using ....IntegralEquations
using ..Octree
using ..OctreeBuilder
using ..Translation
using ..Disaggregation: disaggregate_downward!
using ..Aggregation: aggregate_upward!
using ..Level
using ..Interpolation: AbstractPolesInfo

using ....IntegralEquations.Impedance: get_triangle_info, get_triangles_info
using ....IntegralEquations.PMCHWModule: PMCHW

using ....IntegralEquations.PMCHWModule: assemble_impedance_matrix
const pmchw_assemble_full = assemble_impedance_matrix

export PMCHWMLFMAOperator, assemble_near_field_pmchw

# ─────────────────────────────────────────────────────────────────────────────
# Struct
# ─────────────────────────────────────────────────────────────────────────────

"""
    PMCHWMLFMAOperator{FT,CT} <: AbstractIntegralOperator

PMCHW 系统的 MLFMA 线性算子，大小 2N×2N。
"""
struct PMCHWMLFMAOperator{FT<:AbstractFloat,CT<:Complex} <: AbstractIntegralOperator
    pmchw           :: PMCHW{FT,CT}
    basis           :: RWGBasis
    Z_near          :: SparseMatrixCSC{CT,Int}
    octree0         :: OctreeInfo
    octree1         :: OctreeInfo
    sorted_ids0     :: Vector{Int}
    inv_sorted_ids0 :: Vector{Int}
    sorted_ids1     :: Vector{Int}
    inv_sorted_ids1 :: Vector{Int}
    freq            :: FT
end

Base.size(op::PMCHWMLFMAOperator)         = (2 * num_basis(op.basis), 2 * num_basis(op.basis))
Base.size(op::PMCHWMLFMAOperator, d::Int) = 2 * num_basis(op.basis)
Base.eltype(::PMCHWMLFMAOperator{FT,CT}) where {FT,CT} = CT

# ─────────────────────────────────────────────────────────────────────────────
# 构造函数
# ─────────────────────────────────────────────────────────────────────────────

function PMCHWMLFMAOperator(pmchw::PMCHW, basis::RWGBasis, leaf_size::Float64)
    N  = num_basis(basis)
    FT = typeof(real(pmchw.freq))
    CT = Complex{FT}

    centers = reduce(hcat, [bf.center for bf in basis.functions])  # 3 × N

    λ0 = FT(2π) / real(pmchw.k0)
    λ1 = FT(2π) / real(pmchw.k1)

    octree0, sorted_ids0 = build_octree(centers, leaf_size; λ = Float64(λ0))
    octree1, sorted_ids1 = build_octree(centers, leaf_size; λ = Float64(λ1))

    inv_sorted_ids0 = Vector{Int}(undef, N)
    inv_sorted_ids1 = Vector{Int}(undef, N)
    for i in 1:N
        inv_sorted_ids0[sorted_ids0[i]] = i
        inv_sorted_ids1[sorted_ids1[i]] = i
    end

    Z_near = assemble_near_field_pmchw(pmchw, basis, octree0, sorted_ids0, inv_sorted_ids0)

    return PMCHWMLFMAOperator{FT,CT}(
        pmchw, basis, Z_near,
        octree0, octree1,
        sorted_ids0, inv_sorted_ids0,
        sorted_ids1, inv_sorted_ids1,
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
            isdefined(lv, :aggS)    && fill!(lv.aggS,    zero(eltype(lv.aggS)))
            isdefined(lv, :disaggG) && fill!(lv.disaggG, zero(eltype(lv.disaggG)))
        end
    end

    function up_translate_down!(oct)
        for levelID in (oct.nLevels - 1):-1:2
            aggregate_upward!(oct.levels[levelID], oct.levels[levelID + 1])
        end
        for levelID in 2:oct.nLevels
            translate!(oct.levels[levelID])
        end
        for levelID in 2:(oct.nLevels - 1)
            disaggregate_downward!(oct.levels[levelID], oct.levels[levelID + 1])
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
    octree     :: OctreeInfo,
    basis      :: RWGBasis,
    x          :: AbstractVector,
    sorted_ids :: Vector{Int},
    x_range    :: AbstractUnitRange{Int},
    k          :: Number,
)
    leaf_level = octree.levels[octree.nLevels]
    offset     = first(x_range) - 1

    FT  = eltype(leaf_level.cubeEdgel)
    CT  = Complex{FT}
    JK  = CT(im * k)

    poles    = leaf_level.poles
    nPoles   = length(poles.r̂sθsϕs)
    nCubes   = leaf_level.nCubes

    if !isdefined(leaf_level, :aggS)
        leaf_level.aggS = zeros(CT, nPoles, 2, nCubes)
    else
        fill!(leaf_level.aggS, zero(CT))
    end

    tri_info   = get_triangles_info(basis.mesh, basis)
    gq         = GaussQuadratureInfo(:Triangle, 3, FT)
    n_qp       = length(gq.weight)

    poles_r̂    = [p.r̂    for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]

    Threads.@threads for iCube in 1:nCubes
        cube       = leaf_level.cubes[iCube]
        cubeCenter = cube.center

        for bfID_sorted in cube.bfInterval
            bfID_orig = sorted_ids[bfID_sorted]
            coeff     = x[bfID_orig + offset]
            abs(coeff) < 1e-12 && continue

            bf = basis.functions[bfID_orig]

            for i_supp in 1:2
                tri_idx = bf.support[i_supp]
                tri_idx == 0 && continue

                tri    = tri_info[tri_idx]
                v_all  = tri.vertices

                local_edge = bf.local_edge_idx[i_supp]
                v_opp      = v_all[:, local_edge]
                sign_supp  = bf.signs[i_supp]

                for i_qp in 1:n_qp
                    L          = gq.coordinate[:, i_qp]
                    r          = v_all * L
                    rho        = r - v_opp
                    r_local    = r - cubeCenter
                    factor_vec = sign_supp * bf.edge_length / 2 * gq.weight[i_qp] * coeff

                    for iPole in 1:nPoles
                        r̂    = poles_r̂[iPole]
                        phase = exp(JK * dot(r̂, r_local))
                        vec   = rho * factor_vec * phase

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
    bf    :: RWG,
    basis :: RWGBasis,
    field :: AbstractArray,
    r0    :: AbstractVector,
    k     :: Number,
    η     :: Number,
    poles :: AbstractPolesInfo,
)
    CT        = typeof(complex(one(real(k))))
    FT        = real(CT)
    te        = zero(CT)
    tm        = zero(CT)
    JK        = CT(im * k)
    inv_η     = CT(1 / η)

    poles_r̂    = [p.r̂    for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]
    nPoles     = length(poles_r̂)

    tri_info = get_triangles_info(basis.mesh, basis)
    gq       = GaussQuadratureInfo(:Triangle, 3, FT)
    n_qp     = length(gq.weight)

    for i_supp in 1:2
        tri_idx = bf.support[i_supp]
        tri_idx == 0 && continue

        tri   = tri_info[tri_idx]
        v_all = tri.vertices

        v1 = v_all[:, 1]; v2 = v_all[:, 2]; v3 = v_all[:, 3]
        normal = normalize(cross(v2 - v1, v3 - v1))

        local_edge = bf.local_edge_idx[i_supp]
        v_opp      = v_all[:, local_edge]
        sign_supp  = bf.signs[i_supp]

        for i_qp in 1:n_qp
            L        = gq.coordinate[:, i_qp]
            r        = v_all * L
            rho      = r - v_opp
            r_local  = r - r0
            w_f      = sign_supp * bf.edge_length / 2 * gq.weight[i_qp]

            for iPole in 1:nPoles
                r̂    = poles_r̂[iPole]
                θhat = poles_θhat[iPole]
                ϕhat = poles_ϕhat[iPole]
                phase = exp(-JK * dot(r̂, r_local))

                Eθ    = CT(field[iPole, 1])
                Eϕ    = CT(field[iPole, 2])

                E_inc = (Eθ * θhat + Eϕ * ϕhat) * phase
                te   += dot(rho, E_inc) * w_f

                H_inc  = (Eθ * ϕhat - Eϕ * θhat) * (phase * inv_η)
                tm    += dot(cross(rho, normal), H_inc) * w_f
            end
        end
    end
    return te, tm
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.10: disaggregate_leaf_pmchw_j!
# ─────────────────────────────────────────────────────────────────────────────

function disaggregate_leaf_pmchw_j!(
    octree     :: OctreeInfo,
    basis      :: RWGBasis,
    pmchw      :: PMCHW,
    y          :: AbstractVector,
    sorted_ids :: Vector{Int},
    kmode      :: Symbol,
)
    N    = num_basis(basis)
    k, η = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)

    factor_EJ = im * k * η / (4π)
    factor_HJ = -im * k      / (4π)

    leaf_level = octree.levels[octree.nLevels]
    isdefined(leaf_level, :disaggG) || return

    Threads.@threads for iCube in 1:leaf_level.nCubes
        cube  = leaf_level.cubes[iCube]
        field = view(leaf_level.disaggG, :, :, iCube)
        r0    = cube.center

        for bfID_sorted in cube.bfInterval
            bfID = sorted_ids[bfID_sorted]
            bf   = basis.functions[bfID]

            te, tm = _receive_terms(bf, basis, field, r0, k, η, leaf_level.poles)
            y[bfID]     += te * factor_EJ
            y[bfID + N] += tm * factor_HJ
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.10: disaggregate_leaf_pmchw_m!
# ─────────────────────────────────────────────────────────────────────────────

function disaggregate_leaf_pmchw_m!(
    octree     :: OctreeInfo,
    basis      :: RWGBasis,
    pmchw      :: PMCHW,
    y          :: AbstractVector,
    sorted_ids :: Vector{Int},
    kmode      :: Symbol,
)
    N    = num_basis(basis)
    k, η = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)

    factor_EM =  im * k / (4π)
    factor_HM =  im * k / (η * 4π)

    leaf_level = octree.levels[octree.nLevels]
    isdefined(leaf_level, :disaggG) || return

    Threads.@threads for iCube in 1:leaf_level.nCubes
        cube  = leaf_level.cubes[iCube]
        field = view(leaf_level.disaggG, :, :, iCube)
        r0    = cube.center

        for bfID_sorted in cube.bfInterval
            bfID = sorted_ids[bfID_sorted]
            bf   = basis.functions[bfID]

            te, tm = _receive_terms(bf, basis, field, r0, k, η, leaf_level.poles)
            y[bfID]     += tm * factor_EM
            y[bfID + N] += te * factor_HM
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.8: assemble_near_field_pmchw
# ─────────────────────────────────────────────────────────────────────────────

function assemble_near_field_pmchw(
    pmchw          :: PMCHW,
    basis          :: RWGBasis,
    octree         :: OctreeInfo,
    sorted_ids     :: Vector{Int},
    inv_sorted_ids :: Vector{Int},
)
    N  = num_basis(basis)
    CT = Complex{typeof(pmchw.k0)}

    println("  [PMCHWMLFMAOperator] 装配完整 2N×2N PMCHW 矩阵（N=$N）...")
    Z_full = pmchw_assemble_full(pmchw, basis)

    leaf_level = octree.levels[octree.nLevels]
    nCubes     = leaf_level.nCubes

    # 构建近邻对集合
    near_pairs = Set{Tuple{Int,Int}}()
    for i_cube in 1:nCubes
        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        my_bfs = [sorted_ids[s] for s in cube.bfInterval]

        for neigh_idx in cube.neighbors
            neigh_cube = leaf_level.cubes[neigh_idx]
            isempty(neigh_cube.bfInterval) && continue
            neigh_bfs = [sorted_ids[s] for s in neigh_cube.bfInterval]

            for i in my_bfs, j in neigh_bfs
                push!(near_pairs, (i,     j    ))   # EJ 块
                push!(near_pairs, (i,     j + N))   # EM 块
                push!(near_pairs, (i + N, j    ))   # HJ 块
                push!(near_pairs, (i + N, j + N))   # HM 块
            end
        end
    end

    nz  = length(near_pairs)
    IIs = Vector{Int}(undef, nz)
    JJs = Vector{Int}(undef, nz)
    VVs = Vector{CT}(undef, nz)

    idx = 1
    for (row, col) in near_pairs
        IIs[idx] = row
        JJs[idx] = col
        VVs[idx] = CT(Z_full[row, col])
        idx += 1
    end

    Z_near = sparse(IIs, JJs, VVs, 2N, 2N)
    println("  [PMCHWMLFMAOperator] 近场矩阵：nnz=$(nnz(Z_near)) / $(2N)×$(2N)")
    return Z_near
end

end # module PMCHWMLFMAOperatorModule
