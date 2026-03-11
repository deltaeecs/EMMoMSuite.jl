"""
test_pmchw_block_fidelity_medium.jl — Phase 15 PMCHW medium block/operator fidelity gate

目的：
  - 固化 medium `N=540` 下 block/operator fidelity 的最新诊断结论；
  - 明确 strong-form 在映回物理空间后与 weak-form 具有相同的 fidelity 指标；
  - 锁定当前主导边界为 `M_only -> E-row`，并继续收敛到 `M×k0` 与 `M×k1` 两条 pass 各自都已失真。

说明：
  - 该回归运行成本较高，不并入默认 runtests；
  - 它是 Phase 15 medium long-Krylov backend fidelity 的专门结构化验收入口。
"""

using Test
using EMSuite
using LinearAlgebra
using SparseArrays
using Random
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator, aggregate_leaf_pmchw!, disaggregate_leaf_pmchw_m!
using EMSuite.FastAlgorithms.MLFMA.Aggregation: aggregate_upward!
using EMSuite.FastAlgorithms.MLFMA.Disaggregation: disaggregate_downward!
using EMSuite.FastAlgorithms.MLFMA.Translation: translate!
using EMSuite.IntegralEquations.PMCHWModule: assemble_K_pmchw_offdiag

function make_block_fidelity_medium_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    pairing = pmchw_block_pairing_matrix(basis)
    dense_shell = PMCHWBlockOperator(
        pmchw,
        basis,
        DensePMCHWBackend(Z_dense);
        test_pairing = pairing,
        trial_pairing = pairing,
    )
    return pmchw, basis, Z_dense, pairing, dense_shell
end

function make_probe(N, which::Symbol; seed = 42)
    Random.seed!(seed)
    x = zeros(ComplexF64, 2N)
    if which === :J
        x[1:N] .= randn(ComplexF64, N)
    elseif which === :M
        x[(N + 1):(2N)] .= randn(ComplexF64, N)
    else
        error("unsupported probe selector: $which")
    end
    x ./= norm(x)
    return x
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

function operator_for_mode(shell, mode::Symbol)
    mode === :weak && return weak_form(shell)
    mode === :strong && return strong_form(shell)
    error("unsupported operator mode: $mode")
end

function mode_probe(shell, x_phys, mode::Symbol)
    mode === :weak && return x_phys
    mode === :strong && return shell.trial_pairing \ x_phys
    error("unsupported operator mode: $mode")
end

function recover_mode_output(shell, y_mode, mode::Symbol)
    mode === :weak && return y_mode
    mode === :strong && return shell.test_pairing * y_mode
    error("unsupported operator mode: $mode")
end

function clear_agg!(oct)
    for (_, lv) in oct.levels
        if isdefined(lv, :aggS)
            fill!(lv.aggS, zero(eltype(lv.aggS)))
        end
        if isdefined(lv, :disaggG)
            fill!(lv.disaggG, zero(eltype(lv.disaggG)))
        end
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

function evaluate_budget_case(pmchw, basis, Z_dense, pairing, dense_shell, name, budget, probes)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    fast_shell = PMCHWBlockOperator(
        pmchw,
        basis,
        MatrixFreePMCHWBackend(op);
        test_pairing = pairing,
        trial_pairing = pairing,
        block_source = Z_dense,
    )

    N = num_basis(basis)
    mode_rows = Dict{Tuple{Symbol,String},NamedTuple}()
    for mode in (:weak, :strong)
        dense_op = operator_for_mode(dense_shell, mode)
        fast_op = operator_for_mode(fast_shell, mode)
        for (probe_name, x_phys) in probes
            x_mode_dense = mode_probe(dense_shell, x_phys, mode)
            x_mode_fast = mode_probe(fast_shell, x_phys, mode)
            y_dense = recover_mode_output(dense_shell, dense_op * x_mode_dense, mode)
            y_fast = recover_mode_output(fast_shell, fast_op * x_mode_fast, mode)
            mode_rows[(mode, probe_name)] = split_metrics(y_fast, y_dense, N)
        end
    end

    refs = direct_m_pass_refs(pmchw, basis, op)
    pass_rows = Dict{Symbol,NamedTuple}()
    x_m = probes["M_only"]
    x_m_coeff = view(x_m, (N + 1):(2N))
    for kmode in (:k0, :k1)
        octree = kmode === :k0 ? op.octree0 : op.octree1
        sorted_ids = kmode === :k0 ? op.sorted_ids0 : op.sorted_ids1
        k = kmode === :k0 ? pmchw.k0 : pmchw.k1
        y_fast = zeros(ComplexF64, 2N)
        clear_agg!(octree)
        aggregate_leaf_pmchw!(octree, basis, x_m, sorted_ids, (N + 1):(2N), k)
        up_translate_down!(octree)
        disaggregate_leaf_pmchw_m!(octree, basis, pmchw, y_fast, sorted_ids, kmode)
        ref = kmode === :k0 ? refs.k0 : refs.k1
        y_dense = vcat(ref.EM * x_m_coeff, ref.HM * x_m_coeff)
        pass_rows[kmode] = split_metrics(y_fast, y_dense, N)
    end

    return (
        near_density = nnz(op.Z_near) / (2N)^2,
        weak_j = mode_rows[(:weak, "J_only")],
        weak_m = mode_rows[(:weak, "M_only")],
        strong_j = mode_rows[(:strong, "J_only")],
        strong_m = mode_rows[(:strong, "M_only")],
        m_k0 = pass_rows[:k0],
        m_k1 = pass_rows[:k1],
    )
end

@testset "15.BF1 PMCHW medium block fidelity gate" begin
    pmchw, basis, Z_dense, pairing, dense_shell = make_block_fidelity_medium_fixture()
    probes = Dict(
        "J_only" => make_probe(num_basis(basis), :J; seed = 42),
        "M_only" => make_probe(num_basis(basis), :M; seed = 43),
    )

    default_case = evaluate_budget_case(
        pmchw,
        basis,
        Z_dense,
        pairing,
        dense_shell,
        "default",
        PMCHWMLFMAErrorBudget(Float64),
        probes,
    )
    loose_case = evaluate_budget_case(
        pmchw,
        basis,
        Z_dense,
        pairing,
        dense_shell,
        "loose_near",
        PMCHWMLFMAErrorBudget(Float64;
            near_range_scale = 4.0,
            min_near_range = 4,
            max_near_range = 16,
        ),
        probes,
    )

    @info "PMCHW medium block fidelity gate" N=num_basis(basis) default_near_density=default_case.near_density loose_near_density=loose_case.near_density default_weak_j=default_case.weak_j default_weak_m=default_case.weak_m default_m_k0=default_case.m_k0 default_m_k1=default_case.m_k1 loose_weak_m=loose_case.weak_m loose_m_k0=loose_case.m_k0 loose_m_k1=loose_case.m_k1

    @test num_basis(basis) == 540
    @test default_case.near_density > 0.85
    @test loose_case.near_density < 0.35

    @test abs(default_case.weak_j.rel_total - default_case.strong_j.rel_total) < 1e-12
    @test abs(default_case.weak_m.rel_total - default_case.strong_m.rel_total) < 1e-12
    @test abs(loose_case.weak_j.rel_total - loose_case.strong_j.rel_total) < 1e-12
    @test abs(loose_case.weak_m.rel_total - loose_case.strong_m.rel_total) < 1e-12

    @test default_case.weak_j.rel_total < 1e-3
    @test default_case.weak_j.rel_e < 1e-3

    @test default_case.weak_m.rel_total < 1e-3
    @test default_case.weak_m.rel_e < 1e-3
    @test default_case.weak_m.rel_h < 1e-3

    @test loose_case.weak_m.rel_total < 5e-3
    @test loose_case.weak_m.rel_e < 5e-3
    @test loose_case.weak_m.rel_h < 1e-2
    @test loose_case.weak_m.rel_e > 3 * default_case.weak_m.rel_e

    @test default_case.m_k0.rel_e < 1e-3
    @test default_case.m_k1.rel_e < 2e-3
    @test default_case.m_k0.rel_h < 1e-2
    @test default_case.m_k1.rel_h < 1e-2

    @test loose_case.m_k0.rel_e < 2e-3
    @test loose_case.m_k1.rel_e < 5e-3
    @test loose_case.m_k0.rel_h < 1e-2
    @test loose_case.m_k1.rel_h < 1e-2
    @test loose_case.m_k0.rel_e > default_case.m_k0.rel_e
    @test loose_case.m_k1.rel_e > default_case.m_k1.rel_e
end