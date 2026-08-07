"""
test_pmchw_multilevel_quadrature_regression.jl — PMCHW 多层 source/receive quadrature 回归

目的：
  - 用一个明确 `nLevels > 3` 的 medium 球面夹具锁定 PMCHW `M-pass` 的多层回归；
  - 直接对照 far-only dense `EM/HM` 行块，避免只看整块 matvec 时掩盖 source/receive quadrature 失配；
  - 覆盖 default / loose 两档 near-range，防止 leaf source aggregation 再退回旧的 3 点规则。

说明：
  - 本测试聚焦 `M_only` probe 的 `M×k0` / `M×k1` 两条 pass；
  - medium 夹具运行成本高于常规单元测试，因此提供独立 batch 入口，不并入默认 `runtests.jl`。
"""

using Test
using EMMoMSuite
using LinearAlgebra
using Random
using SparseArrays

using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator, PMCHWMLFMAErrorBudget, aggregate_leaf_pmchw!, disaggregate_leaf_pmchw_m!
using EMMoMSuite.FastAlgorithms.MLFMA.Aggregation: aggregate_upward!
using EMMoMSuite.FastAlgorithms.MLFMA.Disaggregation: disaggregate_downward!
using EMMoMSuite.FastAlgorithms.MLFMA.Translation: translate!
using EMMoMSuite.IntegralEquations.PMCHWModule: assemble_K_pmchw_offdiag

function make_multilevel_pmchw_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    return pmchw, basis
end

function make_probe(N, which::Symbol; seed = 43)
    Random.seed!(seed)
    x = zeros(ComplexF64, 2N)
    if which === :M
        x[(N + 1):(2N)] .= randn(ComplexF64, N)
    else
        error("unsupported probe selector: $which")
    end
    x ./= norm(x)
    return x
end

function clear_agg!(oct)
    for (_, lv) in oct.levels
        isdefined(lv, :aggS) && fill!(lv.aggS, zero(eltype(lv.aggS)))
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

function near_pairs_basis(octree, sorted_ids)
    leaf_level = octree.levels[octree.nLevels]
    pairs = Set{Tuple{Int,Int}}()

    for i_cube in 1:leaf_level.nCubes
        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        test_ids = [sorted_ids[s] for s in cube.bfInterval]

        for neigh_idx in cube.neighbors
            neigh_cube = leaf_level.cubes[neigh_idx]
            isempty(neigh_cube.bfInterval) && continue
            src_ids = [sorted_ids[s] for s in neigh_cube.bfInterval]
            for i in test_ids, j in src_ids
                push!(pairs, (i, j))
            end
        end
    end

    return pairs
end

function direct_m_pass_refs(pmchw, basis, op)
    k1_real = abs(real(pmchw.k1))
    eta1_real = abs(real(pmchw.eta1))
    em_k0 = assemble_K_pmchw_offdiag(basis, pmchw.k0)
    em_k1 = assemble_K_pmchw_offdiag(basis, k1_real)
    hm_k0 = assemble_impedance_matrix(efie_from_keta(pmchw.k0, pmchw.eta0, im * pmchw.k0 / (pmchw.eta0 * 16pi)), basis)
    hm_k1 = assemble_impedance_matrix(efie_from_keta(k1_real, eta1_real, im * pmchw.k1 / (pmchw.eta1 * 16pi)), basis)

    for (i, j) in near_pairs_basis(op.octree0, op.sorted_ids0)
        em_k0[i, j] = 0
        hm_k0[i, j] = 0
    end
    for (i, j) in near_pairs_basis(op.octree1, op.sorted_ids1)
        em_k1[i, j] = 0
        hm_k1[i, j] = 0
    end

    return (
        k0 = (EM = em_k0, HM = hm_k0),
        k1 = (EM = em_k1, HM = hm_k1),
    )
end

function relcorr(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, corr
end

function split_metrics(y_fast, y_dense, N)
    rel_total, corr_total = relcorr(y_fast, y_dense)
    rel_e, corr_e = relcorr(view(y_fast, 1:N), view(y_dense, 1:N))
    rel_h, corr_h = relcorr(view(y_fast, (N + 1):(2N)), view(y_dense, (N + 1):(2N)))
    return (
        rel_total = rel_total,
        corr_total = corr_total,
        rel_e = rel_e,
        corr_e = corr_e,
        rel_h = rel_h,
        corr_h = corr_h,
    )
end

function evaluate_m_pass_case(pmchw, basis, budget)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    refs = direct_m_pass_refs(pmchw, basis, op)

    N = num_basis(basis)
    x_phys = make_probe(N, :M; seed = 43)
    x_m = view(x_phys, (N + 1):(2N))

    pass_rows = Dict{Symbol,NamedTuple}()
    for kmode in (:k0, :k1)
        octree = kmode === :k0 ? op.octree0 : op.octree1
        sorted_ids = kmode === :k0 ? op.sorted_ids0 : op.sorted_ids1
        k = kmode === :k0 ? pmchw.k0 : pmchw.k1
        ref = kmode === :k0 ? refs.k0 : refs.k1

        y_fast = zeros(ComplexF64, 2N)
        clear_agg!(octree)
        aggregate_leaf_pmchw!(octree, basis, x_phys, sorted_ids, (N + 1):(2N), k)
        up_translate_down!(octree)
        disaggregate_leaf_pmchw_m!(octree, basis, pmchw, y_fast, sorted_ids, kmode)

        y_dense = vcat(ref.EM * x_m, ref.HM * x_m)
        pass_rows[kmode] = split_metrics(y_fast, y_dense, N)
    end

    return (
        near_density = nnz(op.Z_near) / (2N)^2,
        nlevels0 = op.octree0.nLevels,
        nlevels1 = op.octree1.nLevels,
        m_k0 = pass_rows[:k0],
        m_k1 = pass_rows[:k1],
    )
end

@testset "PMCHW multilevel quadrature regression" begin
    pmchw, basis = make_multilevel_pmchw_fixture()
    N = num_basis(basis)

    default_case = evaluate_m_pass_case(pmchw, basis, PMCHWMLFMAErrorBudget(Float64))
    loose_case = evaluate_m_pass_case(
        pmchw,
        basis,
        PMCHWMLFMAErrorBudget(Float64; near_range_scale = 4.0, min_near_range = 4, max_near_range = 16),
    )

    @info "PMCHW multilevel quadrature regression" N default_near_density=default_case.near_density loose_near_density=loose_case.near_density default_nlevels=(default_case.nlevels0, default_case.nlevels1) loose_nlevels=(loose_case.nlevels0, loose_case.nlevels1) default_m_k0=default_case.m_k0 default_m_k1=default_case.m_k1 loose_m_k0=loose_case.m_k0 loose_m_k1=loose_case.m_k1

    @test N == 540
    @test default_case.nlevels0 > 3
    @test default_case.nlevels1 > 3
    @test loose_case.nlevels0 > 3
    @test loose_case.nlevels1 > 3

    @test default_case.m_k0.rel_total < 1e-3
    @test default_case.m_k1.rel_total < 5e-3
    @test default_case.m_k0.rel_e < 1e-3
    @test default_case.m_k1.rel_e < 5e-3
    @test default_case.m_k0.rel_h < 1e-2
    @test default_case.m_k1.rel_h < 1e-2

    @test loose_case.m_k0.rel_total < 1e-2
    @test loose_case.m_k1.rel_total < 5e-3
    @test loose_case.m_k0.rel_e < 1e-2
    @test loose_case.m_k1.rel_e < 5e-3
    @test loose_case.m_k0.rel_h < 1e-2
    @test loose_case.m_k1.rel_h < 1e-2
end