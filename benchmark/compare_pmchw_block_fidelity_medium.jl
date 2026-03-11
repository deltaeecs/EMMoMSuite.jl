using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using SparseArrays
using Random
using Printf
using CSV
using DataFrames
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator, aggregate_leaf_pmchw!, disaggregate_leaf_pmchw_m!
using EMSuite.FastAlgorithms.MLFMA.Aggregation: aggregate_upward!
using EMSuite.FastAlgorithms.MLFMA.Disaggregation: disaggregate_downward!
using EMSuite.FastAlgorithms.MLFMA.Translation: translate!
using EMSuite.IntegralEquations.PMCHWModule: assemble_K_pmchw_offdiag

function make_fixture()
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
    return pmchw, basis, Z_dense, dense_shell
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
        norm_dense_e = norm(view(y_dense, 1:N)),
        norm_dense_h = norm(view(y_dense, (N + 1):(2N))),
        norm_fast_e = norm(view(y_fast, 1:N)),
        norm_fast_h = norm(view(y_fast, (N + 1):(2N))),
    )
end

function operator_for_mode(shell, mode::Symbol)
    if mode === :weak
        return weak_form(shell)
    elseif mode === :strong
        return strong_form(shell)
    end
    error("unsupported operator mode: $mode")
end

function mode_probe(shell, x_phys, mode::Symbol)
    if mode === :weak
        return x_phys
    elseif mode === :strong
        return shell.trial_pairing \ x_phys
    end
    error("unsupported operator mode: $mode")
end

function recover_mode_output(shell, y_mode, mode::Symbol)
    if mode === :weak
        return y_mode
    elseif mode === :strong
        return shell.test_pairing * y_mode
    end
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

function direct_m_pass_refs(pmchw, basis, op)
    k1_real = abs(real(pmchw.k1))
    eta1_real = abs(real(pmchw.eta1))
    em_k0 = assemble_K_pmchw_offdiag(basis, pmchw.k0)
    em_k1 = assemble_K_pmchw_offdiag(basis, k1_real)
    hm_k0 = assemble_impedance_matrix(efie_from_keta(pmchw.k0, pmchw.eta0, im * pmchw.k0 / (pmchw.eta0 * 16pi)), basis)
    hm_k1 = assemble_impedance_matrix(efie_from_keta(k1_real, eta1_real, im * pmchw.k1 / (pmchw.eta1 * 16pi)), basis)

    near_k0 = near_pairs_basis(op.octree0, op.sorted_ids0)
    near_k1 = near_pairs_basis(op.octree1, op.sorted_ids1)
    for (i, j) in near_k0
        em_k0[i, j] = 0
        hm_k0[i, j] = 0
    end
    for (i, j) in near_k1
        em_k1[i, j] = 0
        hm_k1[i, j] = 0
    end

    return (
        k0 = (EM = em_k0, HM = hm_k0),
        k1 = (EM = em_k1, HM = hm_k1),
    )
end

function evaluate_m_pass_rows(name, op, pmchw, basis, x_phys, refs)
    N = num_basis(basis)
    rows = NamedTuple[]

    for kmode in (:k0, :k1)
        octree = kmode === :k0 ? op.octree0 : op.octree1
        sorted_ids = kmode === :k0 ? op.sorted_ids0 : op.sorted_ids1
        k = kmode === :k0 ? pmchw.k0 : pmchw.k1
        y_fast = zeros(ComplexF64, 2N)
        clear_agg!(octree)
        aggregate_leaf_pmchw!(octree, basis, x_phys, sorted_ids, (N + 1):(2N), k)
        up_translate_down!(octree)
        disaggregate_leaf_pmchw_m!(octree, basis, pmchw, y_fast, sorted_ids, kmode)

        x_m = view(x_phys, (N + 1):(2N))
        ref = kmode === :k0 ? refs.k0 : refs.k1
        y_dense = vcat(ref.EM * x_m, ref.HM * x_m)
        metrics = split_metrics(y_fast, y_dense, N)

        push!(rows, (
            budget = name,
            mode = "weak",
            probe = "M_$(kmode)_pass",
            leaf_size_eff = op.leaf_size_eff,
            near_range = op.near_range,
            near_density = nnz(op.Z_near) / (length(x_phys)^2),
            rel_total = metrics.rel_total,
            corr_total = metrics.corr_total,
            rel_e = metrics.rel_e,
            corr_e = metrics.corr_e,
            rel_h = metrics.rel_h,
            corr_h = metrics.corr_h,
            norm_dense_e = metrics.norm_dense_e,
            norm_dense_h = metrics.norm_dense_h,
            norm_fast_e = metrics.norm_fast_e,
            norm_fast_h = metrics.norm_fast_h,
        ))
    end

    return rows
end

function evaluate_case(name, budget, pmchw, basis, Z_dense, dense_shell, probes)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    pairing = pmchw_block_pairing_matrix(basis)
    fast_shell = PMCHWBlockOperator(
        pmchw,
        basis,
        MatrixFreePMCHWBackend(op);
        test_pairing = pairing,
        trial_pairing = pairing,
        block_source = Z_dense,
    )
    N = num_basis(basis)
    rows = NamedTuple[]
    m_pass_refs = direct_m_pass_refs(pmchw, basis, op)

    for mode in (:weak, :strong)
        dense_op = operator_for_mode(dense_shell, mode)
        fast_op = operator_for_mode(fast_shell, mode)
        for (probe_name, x_phys) in probes
            x_mode_dense = mode_probe(dense_shell, x_phys, mode)
            x_mode_fast = mode_probe(fast_shell, x_phys, mode)
            y_dense = recover_mode_output(dense_shell, dense_op * x_mode_dense, mode)
            y_fast = recover_mode_output(fast_shell, fast_op * x_mode_fast, mode)
            metrics = split_metrics(y_fast, y_dense, N)

            push!(rows, (
                budget = name,
                mode = String(mode),
                probe = probe_name,
                leaf_size_eff = op.leaf_size_eff,
                near_range = op.near_range,
                near_density = nnz(op.Z_near) / (length(x_phys)^2),
                rel_total = metrics.rel_total,
                corr_total = metrics.corr_total,
                rel_e = metrics.rel_e,
                corr_e = metrics.corr_e,
                rel_h = metrics.rel_h,
                corr_h = metrics.corr_h,
                norm_dense_e = metrics.norm_dense_e,
                norm_dense_h = metrics.norm_dense_h,
                norm_fast_e = metrics.norm_fast_e,
                norm_fast_h = metrics.norm_fast_h,
            ))
        end
    end

    for (probe_name, x_phys) in probes
        probe_name == "M_only" || continue
        append!(rows, evaluate_m_pass_rows(name, op, pmchw, basis, x_phys, m_pass_refs))
    end

    return rows
end

function main()
    pmchw, basis, Z_dense, dense_shell = make_fixture()
    N = num_basis(basis)
    probes = [
        ("J_only", make_probe(N, :J; seed = 42)),
        ("M_only", make_probe(N, :M; seed = 43)),
    ]
    budgets = [
        ("default", PMCHWMLFMAErrorBudget(Float64)),
        (
            "loose_near",
            PMCHWMLFMAErrorBudget(Float64;
                near_range_scale = 4.0,
                min_near_range = 4,
                max_near_range = 16,
            ),
        ),
    ]

    println("PMCHW medium block/operator fidelity comparison")
    @printf("  N = %d (2N = %d)\n", N, 2N)

    rows = NamedTuple[]
    for (name, budget) in budgets
        case_rows = evaluate_case(name, budget, pmchw, basis, Z_dense, dense_shell, probes)
        append!(rows, case_rows)
        println()
        @printf("  budget = %s\n", name)
        for row in case_rows
            @printf(
                "    mode=%-6s probe=%-6s rel_total=%.6e rel_E=%.6e rel_H=%.6e corr_total=%.6f\n",
                row.mode,
                row.probe,
                row.rel_total,
                row.rel_e,
                row.rel_h,
                row.corr_total,
            )
        end
    end

    out_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(out_dir)
    out_path = joinpath(out_dir, "PMCHW_block_fidelity_medium.csv")
    CSV.write(out_path, DataFrame(rows))
    println()
    println("  -> saved $(out_path)")
end

main()