module BlockEvaluatorModule

using StaticArrays
using ...CoreModule: AbstractIntegralOperator, num_elements, num_basis
using ...IntegralEquations
using ...IntegralEquations.Impedance: get_triangles_info
using ...IntegralEquations.EFIEModule: efie_interaction!, EFIE
using ...IntegralEquations.MFIEModule: mfie_interaction!, MFIE
using ...IntegralEquations.CFIEModule: CFIE
using ...IntegralEquations.PMCHWModule: PMCHW, calc_k_pmchw_term!
using ...IntegralEquations.EFIEModule: efie_from_keta
using ...BasisFunctions: RWGBasis
using ...Geometry: GaussQuadratureInfo, get_global_quad_points

export BlockEvaluator, PMCHWBlockEvaluator, eval_block

"""
    BlockEvaluator(op, basis)

预计算三角形信息与基函数→三角形映射，供 `eval_block` 按行/列集合求值阻抗块。
当前支持 RWG 单基函数表面算子（EFIE/MFIE/CFIE）。

# 字段
- `tri_info`：每个三角形的几何/求积信息（`TriangleInfo`）。
- `tri_to_rwg`：三角形 → 支撑其上的 RWG 基函数 `(local_edge, global_id, sign)`。
- `basis_tris`：基函数 → 支撑三角形 id 列表。
- `basis_local`：基函数 → `(tri_id, local_edge, sign)`（预留给多基函数扩展）。
"""
struct BlockEvaluator{OP,BF}
    op::OP
    basis::BF
    tri_info::Vector
    tri_to_rwg::Vector{Vector{Tuple{Int,Int,Float64}}}
    basis_tris::Vector{Vector{Int}}
    basis_local::Vector{Vector{Tuple{Int,Int,Float64}}}
end

function BlockEvaluator(op::AbstractIntegralOperator, basis::RWGBasis)
    nt = num_elements(basis.mesh)
    tri_to_rwg = [Vector{Tuple{Int,Int,Float64}}() for _ in 1:nt]
    basis_tris = [Int[] for _ in 1:num_basis(basis)]
    basis_local = [Tuple{Int,Int,Float64}[] for _ in 1:num_basis(basis)]
    for (i, f) in enumerate(basis.functions)
        for k in 1:2
            t = f.support[k]
            if t > 0
                push!(tri_to_rwg[t], (f.local_edge_idx[k], i, f.signs[k]))
                push!(basis_tris[i], t)
                push!(basis_local[i], (t, f.local_edge_idx[k], f.signs[k]))
            end
        end
    end
    tri_info = get_triangles_info(basis.mesh, basis)
    return BlockEvaluator(op, basis, tri_info, tri_to_rwg, basis_tris, basis_local)
end

"""
    eval_block(ev::BlockEvaluator, rows::Vector{Int}, cols::Vector{Int}) -> Matrix{ComplexF64}

求值 `Z[rows, cols]`（全局基函数索引）。按三角形对去重计算后分发，
支持 EFIE/MFIE/CFIE 组合。EFIE/MFIE 的近场规范序转置修复由
`efie_interaction!`/`mfie_interaction!` 内部处理，与稠密装配路径一致。
"""
function eval_block(ev::BlockEvaluator, rows::Vector{Int}, cols::Vector{Int})
    op = ev.op
    m, n = length(rows), length(cols)
    Z = zeros(ComplexF64, m, n)
    rowmap = Dict(r => i for (i, r) in enumerate(rows))
    colmap = Dict(c => j for (j, c) in enumerate(cols))

    row_tris = unique(vcat((ev.basis_tris[i] for i in rows)...))
    col_tris = unique(vcat((ev.basis_tris[j] for j in cols)...))

    Zloc = zeros(ComplexF64, 3, 3)
    needs_mfie = op isa CFIE
    Zmfie = needs_mfie ? zeros(ComplexF64, 3, 3) : Zloc

    for t1 in row_tris
        for t2 in col_tris
            fill!(Zloc, 0)
            if op isa CFIE
                efie_interaction!(Zloc, op.efie, ev.tri_info[t1], ev.tri_info[t2])
                fill!(Zmfie, 0)
                mfie_interaction!(Zmfie, op.mfie, ev.tri_info[t1], ev.tri_info[t2])
                @inbounds for q in eachindex(Zloc)
                    Zloc[q] = op.alpha * Zloc[q] + (1 - op.alpha) * Zmfie[q]
                end
            elseif op isa MFIE
                mfie_interaction!(Zloc, op, ev.tri_info[t1], ev.tri_info[t2])
            else
                efie_interaction!(Zloc, op, ev.tri_info[t1], ev.tri_info[t2])
            end

            # 分发到 (rows, cols) 块
            for (lt1, g1, s1) in ev.tri_to_rwg[t1]
                haskey(rowmap, g1) || continue
                for (lt2, g2, s2) in ev.tri_to_rwg[t2]
                    haskey(colmap, g2) || continue
                    Z[rowmap[g1], colmap[g2]] += Zloc[lt1, lt2] * s1 * s2
                end
            end
        end
    end
    return Z
end

"""
    PMCHWBlockEvaluator(op::PMCHW, basis::RWGBasis)

PMCHW 系统的 2N×2N 块求值器。全局索引 `1:N` 为 J 通道、`N+1:2N` 为 M 通道。
四个子块：
- Z^EJ（J 行 × J 列）：L(k₀)+L(k₁)，EFIE-like（factor=jkη/16π）；
- Z^HM（M 行 × M 列）：Lₑ(k₀)+Lₑ(k₁)（factor=jk/(η·16π)）；
- Z^EM（J 行 × M 列）：K^PMCHW(k₀)+K^PMCHW(k₁)；
- Z^HJ（M 行 × J 列）：-Z^EM。
"""
struct PMCHWBlockEvaluator{OP,BF}
    op::OP
    basis::BF
    N::Int
    tri_info::Vector
    tri_to_rwg::Vector{Vector{Tuple{Int,Int,Float64}}}
    basis_tris::Vector{Vector{Int}}
    ej0
    ej1
    hm0
    hm1
    k0::ComplexF64
    k1::ComplexF64
    gq_k
    quad_points::Vector
end

function PMCHWBlockEvaluator(op::PMCHW, basis::RWGBasis)
    nt = num_elements(basis.mesh)
    tri_to_rwg = [Vector{Tuple{Int,Int,Float64}}() for _ in 1:nt]
    basis_tris = [Int[] for _ in 1:num_basis(basis)]
    for (i, f) in enumerate(basis.functions)
        for k in 1:2
            t = f.support[k]
            if t > 0
                push!(tri_to_rwg[t], (f.local_edge_idx[k], i, f.signs[k]))
                push!(basis_tris[i], t)
            end
        end
    end
    tri_info = get_triangles_info(basis.mesh, basis)

    # L 块算子（复制 PMCHW.jl `_l_block_operator` 的 factor 约定）
    CT = ComplexF64
    k0 = CT(op.k0)
    eta0 = CT(op.eta0)
    k1 = CT(op.k1)
    eta1 = CT(op.eta1)
    ej0 = efie_from_keta(op.k0, op.eta0, im * k0 * eta0 / (16π))
    ej1 = efie_from_keta(op.k1, op.eta1, im * k1 * eta1 / (16π))
    hm0 = efie_from_keta(op.k0, op.eta0, im * k0 / (eta0 * 16π))
    hm1 = efie_from_keta(op.k1, op.eta1, im * k1 / (eta1 * 16π))

    # K 块求积点（与 _assemble_K_pmchw_offdiag! 一致：4 点）
    FT = eltype(basis.mesh.node)
    gq_k = GaussQuadratureInfo(:Triangle, 4, FT)
    N_points = length(gq_k.weight)
    quad_points = Vector{SVector{N_points,SVector{3,FT}}}(undef, nt)
    for t in 1:nt
        v_idx = basis.mesh.triangles[:, t]
        v1 = SVector{3,FT}(basis.mesh.node[:, v_idx[1]])
        v2 = SVector{3,FT}(basis.mesh.node[:, v_idx[2]])
        v3 = SVector{3,FT}(basis.mesh.node[:, v_idx[3]])
        quad_points[t] = SVector{N_points,SVector{3,FT}}(
            v1 * gq_k.coordinate[1, i] + v2 * gq_k.coordinate[2, i] + v3 * gq_k.coordinate[3, i]
            for i in 1:N_points
        )
    end

    return PMCHWBlockEvaluator(
        op,
        basis,
        num_basis(basis),
        tri_info,
        tri_to_rwg,
        basis_tris,
        ej0,
        ej1,
        hm0,
        hm1,
        k0,
        k1,
        gq_k,
        quad_points,
    )
end

"""
    eval_block(ev::PMCHWBlockEvaluator, rows::Vector{Int}, cols::Vector{Int})

求值 PMCHW 系统 `Z[rows, cols]`（全局 2N 索引，行/列可混合 J/M 通道）。
与稠密装配一致：EJ/HM 用 L 算子，EM 用 +K、HJ 用 -K；K 自对（同三角形）跳过。
"""
function eval_block(ev::PMCHWBlockEvaluator, rows::Vector{Int}, cols::Vector{Int})
    N = ev.N
    m, n = length(rows), length(cols)
    Z = zeros(ComplexF64, m, n)

    row_pass = [r <= N ? 1 : 2 for r in rows]
    col_pass = [c <= N ? 1 : 2 for c in cols]
    lrows = [mod1(r, N) for r in rows]
    lcols = [mod1(c, N) for c in cols]

    row_tris = unique(vcat((ev.basis_tris[li] for li in lrows)...))
    col_tris = unique(vcat((ev.basis_tris[lj] for lj in lcols)...))

    # 注意：efie_interaction! 末尾会整体乘以 efie.factor，同一缓冲连续调用
    # 不同算子会把前一次贡献再乘一次因子 → 每个 L 算子必须用独立缓冲再求和。
    Z_ej0 = zeros(ComplexF64, 3, 3)
    Z_ej1 = zeros(ComplexF64, 3, 3)
    Z_hm0 = zeros(ComplexF64, 3, 3)
    Z_hm1 = zeros(ComplexF64, 3, 3)
    Z_k0 = zeros(ComplexF64, 3, 3)
    Z_k1 = zeros(ComplexF64, 3, 3)
    mfie0 = (k = ev.k0, eta = one(ev.k0), gq_far = ev.gq_k, gq_near = ev.gq_k)
    mfie1 = (k = ev.k1, eta = one(ev.k1), gq_far = ev.gq_k, gq_near = ev.gq_k)

    for t1 in row_tris
        ti1 = ev.tri_info[t1]
        for t2 in col_tris
            ti2 = ev.tri_info[t2]

            # L 块：EJ/HM
            fill!(Z_ej0, 0)
            efie_interaction!(Z_ej0, ev.ej0, ti1, ti2)
            fill!(Z_ej1, 0)
            efie_interaction!(Z_ej1, ev.ej1, ti1, ti2)
            fill!(Z_hm0, 0)
            efie_interaction!(Z_hm0, ev.hm0, ti1, ti2)
            fill!(Z_hm1, 0)
            efie_interaction!(Z_hm1, ev.hm1, ti1, ti2)

            # K 块：EM/HJ（自对跳过，与 _assemble_K_pmchw_offdiag! 一致）
            fill!(Z_k0, 0)
            fill!(Z_k1, 0)
            if t1 != t2
                q1 = ev.quad_points[t1]
                q2 = ev.quad_points[t2]
                calc_k_pmchw_term!(Z_k0, mfie0, ti1, ti2, q1, q2)
                calc_k_pmchw_term!(Z_k1, mfie1, ti1, ti2, q1, q2)
            end

            for (lt1, g1, s1) in ev.tri_to_rwg[t1]
                for (lt2, g2, s2) in ev.tri_to_rwg[t2]
                    for (i, r) in enumerate(rows)
                        lrows[i] == g1 || continue
                        for (j, c) in enumerate(cols)
                            lcols[j] == g2 || continue
                            val = if row_pass[i] == 1 && col_pass[j] == 1
                                Z_ej0[lt1, lt2] + Z_ej1[lt1, lt2]
                            elseif row_pass[i] == 2 && col_pass[j] == 2
                                Z_hm0[lt1, lt2] + Z_hm1[lt1, lt2]
                            else
                                v = Z_k0[lt1, lt2] + Z_k1[lt1, lt2]
                                row_pass[i] == 1 ? v : -v
                            end
                            Z[i, j] += val * s1 * s2
                        end
                    end
                end
            end
        end
    end
    return Z
end

end # module BlockEvaluatorModule
