module BlockEvaluatorModule

using ...CoreModule: AbstractIntegralOperator, num_elements, num_basis
using ...IntegralEquations
using ...IntegralEquations.Impedance: get_triangles_info
using ...IntegralEquations.EFIEModule: efie_interaction!, EFIE
using ...IntegralEquations.MFIEModule: mfie_interaction!, MFIE
using ...IntegralEquations.CFIEModule: CFIE
using ...BasisFunctions: RWGBasis

export BlockEvaluator, eval_block

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

end # module BlockEvaluatorModule
