module Aggregation

using LinearAlgebra
using StaticArrays
using ....CoreModule
using ....CoreModule: Constants
using ....Geometry
using ....BasisFunctions
using ....IntegralEquations
using ..Level
using ..Octree
using ..Interpolation

using ....IntegralEquations.Impedance: get_triangle_info, get_triangles_info
using ....IntegralEquations.VEFIEModule: get_tetrahedra_info

export aggregate!, aggregate_leaf!, aggregate_upward!

function get_k(op::AbstractIntegralOperator)
    if hasfield(typeof(op), :k)
        return op.k
    elseif hasfield(typeof(op), :efie)
        return op.efie.k
    elseif hasfield(typeof(op), :freq)
        return 2π * op.freq / Constants.c0
    else
        error("Operator $(typeof(op)) does not have field k or efie or freq")
    end
end

"""
    aggregate!(octree::OctreeInfo, basis::AbstractBasisFunction, operator::AbstractIntegralOperator, x::AbstractVector, sorted_ids::Vector{Int})

Compute aggregation for all levels using coefficients x.
"""
function aggregate!(
    octree::OctreeInfo,
    basis::AbstractBasisFunction,
    operator::AbstractIntegralOperator,
    x::AbstractVector,
    sorted_ids::Vector{Int},
)
    aggregate!(octree, AbstractBasisFunction[basis], [num_basis(basis)], operator, x, sorted_ids)
end

function aggregate!(
    octree::OctreeInfo,
    bases::Vector{<:AbstractBasisFunction},
    offsets::Vector{Int},
    operator::AbstractIntegralOperator,
    x::AbstractVector,
    sorted_ids::Vector{Int};
    element_cache = nothing,
)
    # 1. Leaf level aggregation (Radiation pattern of basis functions)
    leafLevel = octree.levels[octree.nLevels]
    aggregate_leaf!(
        leafLevel,
        bases,
        offsets,
        operator,
        x,
        sorted_ids;
        element_cache = element_cache,
    )

    # 2. Upward pass (Child to Parent)
    for levelID = (octree.nLevels-1):-1:2
        parentLevel = octree.levels[levelID]
        childLevel = octree.levels[levelID+1]
        aggregate_upward!(parentLevel, childLevel)
    end
end

function get_basis_index(global_idx::Int, offsets::Vector{Int})
    idx = searchsortedfirst(offsets, global_idx)
    local_idx = idx == 1 ? global_idx : global_idx - offsets[idx-1]
    return idx, local_idx
end

"""
    aggregate_leaf!(level::LevelInfo, bases::Vector{<:AbstractBasisFunction}, offsets::Vector{Int}, operator::AbstractIntegralOperator, x::AbstractVector, sorted_ids::Vector{Int})

Compute radiation patterns of basis functions at the leaf level.
"""
function aggregate_leaf!(
    level::LevelInfo,
    bases::Vector{<:AbstractBasisFunction},
    offsets::Vector{Int},
    operator::AbstractIntegralOperator,
    x::AbstractVector,
    sorted_ids::Vector{Int};
    cube_filter = nothing,
    element_cache = nothing,
)
    FT = eltype(level.cubeEdgel)
    CT = Complex{FT}

    k = get_k(operator)
    JK = im * k

    poles = level.poles
    nPoles = length(poles.r̂sθsϕs)
    nCubes = level.nCubes

    if !isdefined(level, :aggS)
        level.aggS = zeros(CT, nPoles, 2, nCubes)
    else
        fill!(level.aggS, zero(CT))
    end

    # Precompute element info and quadrature (or reuse the operator-level cache
    # built by `_build_element_cache` — issue #22: rebuilding these on every matvec
    # dominated per-iteration heap allocations)
    local element_infos, gqs
    if element_cache === nothing
        element_infos = Any[]
        gqs = Any[]

        for b in bases
            if b isa RWGBasis
                push!(element_infos, get_triangles_info(b.mesh, b))
                push!(gqs, GaussQuadratureInfo(:Triangle, 4, FT))
            elseif b isa SWGBasis
                if hasfield(typeof(operator), :permittivities)
                    push!(element_infos, get_tetrahedra_info(b.mesh, b, operator.permittivities))
                else
                    # Fallback
                    push!(
                        element_infos,
                        get_tetrahedra_info(b.mesh, b, fill(1.0 + 0im, num_elements(b.mesh))),
                    )
                end
                push!(gqs, GaussQuadratureInfo(:Tetrahedron, 5, FT))
            else
                error("Unsupported basis type")
            end
        end
    else
        element_infos, gqs = element_cache
    end

    # 极点向量数组只构造一次并缓存在层上（issue #22：原实现每次调用重建 3 个数组）
    if !isdefined(level, :polevecs)
        level.polevecs = (
            [p.r̂ for p in poles.r̂sθsϕs],
            [p.θhat for p in poles.r̂sθsϕs],
            [p.ϕhat for p in poles.r̂sθsϕs],
        )
    end
    poles_r̂, poles_θhat, poles_ϕhat = level.polevecs

    # Loop over cubes
    Threads.@threads for iCube = 1:nCubes
        cube_filter !== nothing && !cube_filter(iCube) && continue
        cube = level.cubes[iCube]
        if isempty(cube.bfInterval)
            continue
        end

        cubeCenter = cube.center
        start_idx = first(cube.bfInterval)
        end_idx = last(cube.bfInterval)

        for i_sorted = start_idx:end_idx
            i_orig = sorted_ids[i_sorted]

            # Identify basis
            b_idx, i_local = get_basis_index(i_orig, offsets)
            basis = bases[b_idx]
            bf = basis.functions[i_local]

            coef = x[i_orig]
            if abs(coef) < 1e-12
                continue
            end

            elem_info = element_infos[b_idx]
            gq = gqs[b_idx]

            if basis isa RWGBasis
                add_radiation_pattern_rwg!(
                    level.aggS,
                    iCube,
                    basis,
                    bf,
                    elem_info,
                    gq,
                    k,
                    cubeCenter,
                    poles_r̂,
                    poles_θhat,
                    poles_ϕhat,
                    coef,
                )
            elseif basis isa SWGBasis
                add_radiation_pattern_swg!(
                    level.aggS,
                    iCube,
                    basis,
                    bf,
                    elem_info,
                    gq,
                    k,
                    cubeCenter,
                    poles_r̂,
                    poles_θhat,
                    poles_ϕhat,
                    coef,
                )
            end
        end
    end
end

"""
    add_radiation_pattern_rwg!(aggS, iCube, basis, bf, elem_info, gq, k, cubeCenter, poles_r̂, poles_θhat, poles_ϕhat, coef)

RWG 基函数的叶层辐射积分（论文式 (2-48)：`F_S(f_n^S, k̂)`）。

对第 `n` 个 RWG 基函数（支撑三角形 `T_{n,i}`，`i = 1, 2`），其辐射方向图按定义
为切向投影后的傅里叶积分：

```math
F_S(f_n^S, \\hat{k}) = \\int_{T_n} \\left(\\overline{I} - \\hat{k}\\hat{k}\\right) \\cdot
f_n(r')\\, e^{{\\rm j}k \\hat{k} \\cdot (r' - r_b)}\\, dS'
```

其中 `r_b` 为基函数所在叶层盒子中心。把统一基函数形式
`f_n(r) = a_n^{±}/(2A^{±}) ρ`（`a_n^{±} = ±l_n` 为带符号边长）代入并对支撑
三角形求和，只保留球坐标 `θ`/`φ` 分量（等价于 `(Ī − k̂k̂)` 投影）：

```math
\\mathcal{F}_p(\\hat{k}) = \\sum_{i=1}^{2} s_{n,i}\\, \\frac{l_n}{2}
\\int_{T_{n,i}} \\hat{e}_p(\\hat{k}) \\cdot \\rho_i(r')\\,
e^{{\\rm j}k \\hat{k} \\cdot (r' - r_b)}\\, dS', \\qquad p \\in \\{\\theta, \\phi\\}
```

代码实现：对每个支撑三角形的高斯点 `r_q`，累加
`aggS[iPole, p, iCube] += θ̂/φ̂ · (ρ * s*l/2 * w_q * x_n * e^{jk k̂·(r_q − r_b)})`，
其中 `coef = x_n` 为基函数系数，`gq.weight` 为参考三角形求积权重
（三角形雅可比已含于求积节点构造中）。

# Arguments
- `aggS`: 叶层聚合方向图，尺寸 `(nPoles, 2, nCubes)`（第 2 维为 θ/φ 分量）。
- `iCube`: 基函数所在盒子的全局编号。
- `basis`/`bf`: 基函数集合与当前基函数（含 `support`、`local_edge_idx`、`signs`、`edge_length`）。
- `elem_info`: 预计算的单元几何信息（`TriangleInfo`）。
- `gq`: 三角形高斯求积规则（默认 3 点）。
- `k`: 波数（`r_local` 以 λ 为单位时相位因子为 `e^{jk k̂·(r−r_b)}`）。
- `cubeCenter`: 盒子中心 `r_b`。
- `poles_r̂/θhat/ϕhat`: 球面采样点方向与对应的单位切向量。
- `coef`: 该基函数的系数 `x_n`。
"""
function add_radiation_pattern_rwg!(
    aggS,
    iCube,
    basis,
    bf,
    elem_info,
    gq,
    k,
    cubeCenter,
    poles_r̂,
    poles_θhat,
    poles_ϕhat,
    coef,
)
    JK = im * k
    FT = eltype(gq.coordinate)
    n_qp = length(gq.weight)
    nPoles = length(poles_r̂)

    for i_supp = 1:2
        tri_idx = bf.support[i_supp]
        if tri_idx == 0
            continue
        end

        tri = elem_info[tri_idx]
        tri_vertices = tri.vertices

        local_edge = bf.local_edge_idx[i_supp]
        opp_vertex_idx = local_edge
        # SVector 取列/运算全部栈分配（issue #22：原实现每 qp 产生 3-4 个
        # 堆 Vector，每 pole 再一个堆 Vector，百万级小分配累积成 GC 灾难）
        v_opp = SVector{3,FT}(
            tri_vertices[1, opp_vertex_idx],
            tri_vertices[2, opp_vertex_idx],
            tri_vertices[3, opp_vertex_idx],
        )

        sign = bf.signs[i_supp]

        for i_qp = 1:n_qp
            L = SVector{3,FT}(
                gq.coordinate[1, i_qp],
                gq.coordinate[2, i_qp],
                gq.coordinate[3, i_qp],
            )
            r = tri_vertices * L
            rho = r - v_opp
            r_local = r - SVector{3,FT}(cubeCenter)

            factor_vec = sign * bf.edge_length / 2 * gq.weight[i_qp] * coef

            for iPole = 1:nPoles
                r̂ = poles_r̂[iPole]
                phase = exp(JK * dot(r̂, r_local))
                vec = rho * (factor_vec * phase)

                aggS[iPole, 1, iCube] += dot(poles_θhat[iPole], vec)
                aggS[iPole, 2, iCube] += dot(poles_ϕhat[iPole], vec)
            end
        end
    end
end

"""
    add_radiation_pattern_swg!(aggS, iCube, basis, bf, elem_info, gq, k, cubeCenter, poles_r̂, poles_θhat, poles_ϕhat, coef)

SWG 基函数的叶层辐射积分（VEFIE/VSIE 路径，论文式 (2-48) 的体扩展）。

对 SWG 基函数，其统一形式为 `f_n(r) = a_n^{±}/(3V^{±}) ρ`（`a_n^{±} = ±A_n`
为带符号公共面面积），体等效电流还需乘以对比度 `κ = ε/ε_0`。辐射方向图：

```math
\\mathcal{F}_p(\\hat{k}) = \\kappa \\sum_{i=1}^{2} s_{n,i}\\, \\frac{A_n}{3}
\\int_{V_{n,i}} \\hat{e}_p(\\hat{k}) \\cdot \\rho_i(r')\\, e^{{\\rm j}k \\hat{k} \\cdot (r' - r_b)}\\, dV'
```

代码逐项累加 `aggS[iPole, p, iCube] += θ̂/φ̂ · (ρ * s*A/3 * w_q * x_n * κ * e^{jk k̂·(r_q−r_b)})`。

# Arguments
- `aggS`: 叶层聚合方向图，尺寸 `(nPoles, 2, nCubes)`。
- `basis`/`bf`: 基函数集合与当前 SWG 基函数（`support`、`local_face_idx`、`signs`、`area`）。
- `elem_info`: 预计算的四面体几何信息（`TetrahedraInfo`，含对比度 `κ`）。
- `gq`: 四面体高斯求积规则（默认 5 点）。
- 其余参数同 [`add_radiation_pattern_rwg!`](@ref)。
"""
function add_radiation_pattern_swg!(
    aggS,
    iCube,
    basis,
    bf,
    elem_info,
    gq,
    k,
    cubeCenter,
    poles_r̂,
    poles_θhat,
    poles_ϕhat,
    coef,
)
    JK = im * k
    FT = eltype(gq.coordinate)
    n_qp = length(gq.weight)
    nPoles = length(poles_r̂)

    for i_supp = 1:2
        tet_idx = bf.support[i_supp]
        if tet_idx == 0
            continue
        end

        tet = elem_info[tet_idx]
        tet_vertices = tet.vertices

        # Kappa
        kappa = tet.κ

        local_face = bf.local_face_idx[i_supp]
        v_free = SVector{3,FT}(
            tet_vertices[1, local_face],
            tet_vertices[2, local_face],
            tet_vertices[3, local_face],
        )

        sign = bf.signs[i_supp]

        for i_qp = 1:n_qp
            # 四面体重心坐标是 4 维（SMatrix{4,N}），L 需为 SVector{4}
            L = SVector{4,FT}(
                gq.coordinate[1, i_qp],
                gq.coordinate[2, i_qp],
                gq.coordinate[3, i_qp],
                gq.coordinate[4, i_qp],
            )
            r = tet_vertices * L
            rho = sign * (r - v_free)
            r_local = r - SVector{3,FT}(cubeCenter)

            factor_vec = bf.area / 3 * gq.weight[i_qp] * coef * kappa

            for iPole = 1:nPoles
                r̂ = poles_r̂[iPole]
                phase = exp(JK * dot(r̂, r_local))
                vec = rho * (factor_vec * phase)

                aggS[iPole, 1, iCube] += dot(poles_θhat[iPole], vec)
                aggS[iPole, 2, iCube] += dot(poles_ϕhat[iPole], vec)
            end
        end
    end
end

"""
    aggregate_upward!(parentLevel::LevelInfo, childLevel::LevelInfo)

把子层方向图聚合到父层：插值上采样 + 相移（论文式 (4-2)~(4-3)）。

```math
\\mathcal{F}(\\hat{k}^{l-1}) = \\sum_{c \\in \\mathrm{kids}(p)}
e^{{\\rm j}k \\hat{k} \\cdot (r_c - r_p)}\\, \\bm{\\Gamma}^{l-1,l}\\,
\\mathcal{F}_c(\\hat{k}^{l})
```

其中 `Γ^{l-1,l}` 为子层（第 `l` 层）到父层（第 `l-1` 层）的插值矩阵，
`e^{jk k̂·(r_c − r_p)}` 为子盒到父盒的相移（`phaseShiftFromKids`）。
对 Lebedev 一步插值，`Γ` 的 `θϕCSC` 直接作用在展平的 (θ, φ) 分量上；
对传统两段式，先 `ϕCSC` 再 `θCSC`。最终写入 `parentLevel.aggS`。
"""
function aggregate_upward!(
    parentLevel::LevelInfo,
    childLevel::LevelInfo;
    cube_filter = nothing,
)
    FT = eltype(parentLevel.cubeEdgel)
    CT = Complex{FT}

    nPolesParent = length(parentLevel.poles.r̂sθsϕs)
    nCubesParent = parentLevel.nCubes

    if !isdefined(parentLevel, :aggS)
        parentLevel.aggS = zeros(CT, nPolesParent, 2, nCubesParent)
    else
        fill!(parentLevel.aggS, zero(CT))
    end

    parentAggS = parentLevel.aggS
    childAggS = childLevel.aggS

    interp = childLevel.interpWθϕ

    phaseShift = parentLevel.phaseShiftFromKids

    # Per-chunk scratch buffers: one chunk of cubes per task, each chunk owns its
    # buffers (indexed by the *loop* variable, not Threads.threadid() — the
    # latter can exceed Threads.nthreads() on multi-pool runtimes). Buffers are
    # cached on the level and rebuilt only when the thread count changes —
    # issue #22 第二轮（原实现每次调用/每层重建）。
    nthreads = Threads.nthreads()
    sc = nothing
    if isdefined(parentLevel, :uwScratch) && parentLevel.uwScratch.nthreads == nthreads
        sc = parentLevel.uwScratch
    end
    nScratch = max(1, min(nthreads, nCubesParent))
    chunks = collect(Iterators.partition(1:nCubesParent, cld(nCubesParent, nScratch)))
    if sc === nothing
        if interp isa FFTInterpInfo
            fftbuf = Array{ComplexF64}(undef, interp.nθ * interp.M2, 2, childLevel.nCubes)
            sc = (nthreads = nthreads, fftbuf = fftbuf,
                scratch_agg = [Matrix{CT}(undef, nPolesParent, 2) for _ in 1:length(chunks)])
        elseif hasfield(typeof(interp), :θϕCSC)
            # Lebedev one-step: result length = size(θϕCSC, 1) = 2·nPolesParent
            sc = (nthreads = nthreads, fftbuf = nothing,
                scratch_flat = [Vector{CT}(undef, size(interp.θϕCSC, 1)) for _ in 1:length(chunks)])
        else
            # Traditional two-step: ϕ-interp result, then θ-interp result
            sc = (nthreads = nthreads, fftbuf = nothing,
                scratch_ϕ = [Matrix{CT}(undef, size(interp.ϕCSC, 1), 2) for _ in 1:length(chunks)],
                scratch_agg = [Matrix{CT}(undef, nPolesParent, 2) for _ in 1:length(chunks)])
        end
        parentLevel.uwScratch = sc
    end
    if interp isa FFTInterpInfo
        childPhi = fft_interp_phi_batch!(sc.fftbuf::Array{ComplexF64,3}, childAggS, interp)
        scratch_agg = sc.scratch_agg::Vector{Matrix{CT}}

        # 三种插值路径各自特化整个线程循环体：原先把三分支塞进 if-表达式
        # 使 aggInterp 为 Union 类型，逐 (kid, pol) 装箱广播临时
        # （issue #22 第二轮：PMCHW 实测每 matvec 数 MB）。
        Threads.@threads for ci in 1:length(chunks)
            for iCube in chunks[ci]
                cube_filter !== nothing && !cube_filter(iCube) && continue
                parentCube = parentLevel.cubes[iCube]

                for iKid = 1:length(parentCube.kidsInterval)
                    childID = parentCube.kidsInterval[iKid]
                    childIn8 = parentCube.kidsIn8[iKid]

                    aggInterp = scratch_agg[ci]
                    mul!(aggInterp, interp.θCSC, view(childPhi, :, :, childID))

                    shift = view(phaseShift, :, childIn8)
                    for pol = 1:2
                        @views parentAggS[:, pol, iCube] .+= shift .* aggInterp[:, pol]
                    end
                end
            end
        end
    elseif hasfield(typeof(interp), :θϕCSC)
        scratch_flat = sc.scratch_flat::Vector{Vector{CT}}

        Threads.@threads for ci in 1:length(chunks)
            for iCube in chunks[ci]
                cube_filter !== nothing && !cube_filter(iCube) && continue
                parentCube = parentLevel.cubes[iCube]

                for iKid = 1:length(parentCube.kidsInterval)
                    childID = parentCube.kidsInterval[iKid]
                    childIn8 = parentCube.kidsIn8[iKid]

                    aggChild = view(childAggS, :, :, childID)
                    # Lebedev 一步插值：θϕCSC 直接作用在展平的 (θ,ϕ) 分量上
                    mul!(scratch_flat[ci], interp.θϕCSC, vec(aggChild))
                    shift = view(phaseShift, :, childIn8)

                    fi = 0
                    for pol = 1:2, iP = 1:nPolesParent
                        fi += 1
                        parentAggS[iP, pol, iCube] += shift[iP] * scratch_flat[ci][fi]
                    end
                end
            end
        end
    else
        scratch_ϕ = sc.scratch_ϕ::Vector{Matrix{CT}}
        scratch_agg = sc.scratch_agg::Vector{Matrix{CT}}

        # 传统两段式：先 ϕ 后 θ
        Threads.@threads for ci in 1:length(chunks)
            for iCube in chunks[ci]
                cube_filter !== nothing && !cube_filter(iCube) && continue
                parentCube = parentLevel.cubes[iCube]

                for iKid = 1:length(parentCube.kidsInterval)
                    childID = parentCube.kidsInterval[iKid]
                    childIn8 = parentCube.kidsIn8[iKid]

                    aggChild = view(childAggS, :, :, childID)
                    aggInterp = scratch_agg[ci]
                    mul!(scratch_ϕ[ci], interp.ϕCSC, aggChild)
                    mul!(aggInterp, interp.θCSC, scratch_ϕ[ci])

                    shift = view(phaseShift, :, childIn8)
                    for pol = 1:2
                        @views parentAggS[:, pol, iCube] .+= shift .* aggInterp[:, pol]
                    end
                end
            end
        end
    end
end

end
