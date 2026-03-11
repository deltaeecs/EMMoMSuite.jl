using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using Printf
using CSV
using DataFrames
using SparseArrays
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    rhs = excitation_vector(pmchw, feed, basis)
    dense_shell = PMCHWBlockOperator(pmchw, basis)
    rhs_sf = strong_form_rhs(dense_shell, rhs)
    return pmchw, basis, Z_dense, rhs_sf, dense_shell
end

function relcorr(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, corr
end

function split_metrics(y_fast, y_dense)
    N2 = length(y_dense)
    N = N2 ÷ 2
    rel_total, corr_total = relcorr(y_fast, y_dense)
    rel_e, corr_e = relcorr(view(y_fast, 1:N), view(y_dense, 1:N))
    rel_h, corr_h = relcorr(view(y_fast, N + 1:N2), view(y_dense, N + 1:N2))
    return (
        rel_total = rel_total,
        corr_total = corr_total,
        rel_e = rel_e,
        corr_e = corr_e,
        rel_h = rel_h,
        corr_h = corr_h,
    )
end

function build_shell(pmchw, basis, Z_dense, budget)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    shell = PMCHWBlockOperator(pmchw, basis, MatrixFreePMCHWBackend(op); block_source = Z_dense)
    return op, shell
end

function arnoldi_compare(A_dense, A_fast, v0, steps)
    q = copy(v0)
    q ./= norm(q)
    basis = Vector{Vector{ComplexF64}}()
    rows = NamedTuple[]

    for step in 1:steps
        y_dense = A_dense * q
        y_fast = A_fast * q
        metrics = split_metrics(y_fast, y_dense)
        push!(rows, (
            step = step,
            rel_total = metrics.rel_total,
            corr_total = metrics.corr_total,
            rel_e = metrics.rel_e,
            corr_e = metrics.corr_e,
            rel_h = metrics.rel_h,
            corr_h = metrics.corr_h,
            norm_dense = norm(y_dense),
            norm_fast = norm(y_fast),
        ))

        push!(basis, q)
        w = copy(y_dense)
        for q_prev in basis
            w .-= dot(q_prev, w) * q_prev
        end
        beta = norm(w)
        beta <= 1e-12 && break
        q = w ./ beta
    end

    return rows
end

function main(; steps = 12)
    pmchw, basis, Z_dense, rhs_sf, dense_shell = make_fixture()
    A_dense = strong_form(dense_shell)

    budgets = [
        (
            name = "default",
            budget = PMCHWMLFMAErrorBudget(Float64),
        ),
        (
            name = "loose_near",
            budget = PMCHWMLFMAErrorBudget(Float64;
                near_range_scale = 4.0,
                min_near_range = 4,
                max_near_range = 16,
            ),
        ),
    ]

    println("PMCHW medium Krylov-subspace fidelity comparison")
    @printf("  N = %d (2N = %d)\n", num_basis(basis), 2 * num_basis(basis))
    @printf("  Arnoldi steps = %d\n", steps)

    results = NamedTuple[]
    for budget_case in budgets
        op, shell = build_shell(pmchw, basis, Z_dense, budget_case.budget)
        A_fast = strong_form(shell)
        rows = arnoldi_compare(A_dense, A_fast, rhs_sf, steps)

        rel_total_max = maximum(row.rel_total for row in rows)
        rel_e_max = maximum(row.rel_e for row in rows)
        rel_h_max = maximum(row.rel_h for row in rows)
        @printf("\n  budget = %s\n", budget_case.name)
        @printf("    near_range = %d, near_density = %.4f\n", op.near_range, nnz(op.Z_near) / (length(rhs_sf)^2))
        @printf("    max rel_total = %.6e, max rel_E = %.6e, max rel_H = %.6e\n", rel_total_max, rel_e_max, rel_h_max)

        for row in rows
            push!(results, merge((
                budget = budget_case.name,
                leaf_size_eff = op.leaf_size_eff,
                near_range = op.near_range,
                near_density = nnz(op.Z_near) / (length(rhs_sf)^2),
            ), row))
            @printf("    step=%2d rel_total=%.6e rel_E=%.6e rel_H=%.6e corr_total=%.6f\n", row.step, row.rel_total, row.rel_e, row.rel_h, row.corr_total)
        end
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_krylov_subspace_medium.csv")
    CSV.write(csv_path, DataFrame(results))
    println("\n  -> saved $(csv_path)")
end

main()