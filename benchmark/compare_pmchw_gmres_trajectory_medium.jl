using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMMoMSuite
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using Printf
using CSV
using DataFrames
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    rhs = excitation_vector(pmchw, feed, basis)
    Zin_lu = input_impedance(pmchw, feed, Z_dense \ rhs, basis)
    dense_shell = PMCHWBlockOperator(pmchw, basis)
    return pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell
end

function relcorr(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, corr
end

function solve_strong(A, rhs_sf, shell; restart, maxiter, reltol)
    coeffs_sf, hist = gmres(
        A,
        rhs_sf;
        restart = restart,
        maxiter = maxiter,
        reltol = reltol,
        log = true,
    )
    return recover_trial_coefficients(shell, coeffs_sf), hist.data[:resnorm][end], length(hist.data[:resnorm])
end

function build_shell(pmchw, basis, Z_dense, budget)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    shell = PMCHWBlockOperator(pmchw, basis, MatrixFreePMCHWBackend(op); block_source = Z_dense)
    return op, shell
end

function precompute_dense_trajectory(pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell, checkpoints; restart, reltol)
    dense_rows = Dict{Int, NamedTuple}()
    A_dense = strong_form(dense_shell)
    rhs_sf = strong_form_rhs(dense_shell, rhs)
    for checkpoint in checkpoints
        I_dense, res_dense, niters_dense = solve_strong(A_dense, rhs_sf, dense_shell; restart = restart, maxiter = checkpoint, reltol = reltol)
        Zin_dense = input_impedance(pmchw, feed, I_dense, basis)
        r_dense = rhs - Z_dense * I_dense
        dense_rows[checkpoint] = (
            I = I_dense,
            gmres_res = res_dense,
            niters = niters_dense,
            dense_res = norm(r_dense) / (norm(rhs) + 1e-30),
            Zin = Zin_dense,
            Zin_re_err_pct = abs(real(Zin_dense) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
            Zin_im_err_ohm = abs(imag(Zin_dense) - imag(Zin_lu)),
        )
    end
    return dense_rows
end

function evaluate_case(pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_rows, budget_name, budget, checkpoints; restart, reltol)
    op, shell = build_shell(pmchw, basis, Z_dense, budget)
    A_fast = strong_form(shell)
    rhs_sf = strong_form_rhs(shell, rhs)
    near_density = nnz(op.Z_near) / (length(rhs)^2)
    rows = NamedTuple[]

    for checkpoint in checkpoints
        I_fast, res_fast, niters_fast = solve_strong(A_fast, rhs_sf, shell; restart = restart, maxiter = checkpoint, reltol = reltol)
        dense_row = dense_rows[checkpoint]

        rel_I, corr_I = relcorr(I_fast, dense_row.I)
        r_dense_fast = rhs - Z_dense * I_fast
        dense_res_fast = norm(r_dense_fast) / (norm(rhs) + 1e-30)
        Zin_fast = input_impedance(pmchw, feed, I_fast, basis)

        push!(rows, (
            budget = budget_name,
            checkpoint = checkpoint,
            restart = restart,
            reltol = reltol,
            near_range = op.near_range,
            near_density = near_density,
            niters_dense = dense_row.niters,
            niters_fast = niters_fast,
            gmres_res_dense = dense_row.gmres_res,
            gmres_res_fast = res_fast,
            dense_res_dense = dense_row.dense_res,
            dense_res_fast = dense_res_fast,
            rel_I = rel_I,
            corr_I = corr_I,
            Zin_dense_re = real(dense_row.Zin),
            Zin_dense_im = imag(dense_row.Zin),
            Zin_fast_re = real(Zin_fast),
            Zin_fast_im = imag(Zin_fast),
            Zin_gap_re_ohm = abs(real(Zin_fast) - real(dense_row.Zin)),
            Zin_gap_im_ohm = abs(imag(Zin_fast) - imag(dense_row.Zin)),
            Zin_dense_re_err_pct = dense_row.Zin_re_err_pct,
            Zin_fast_re_err_pct = abs(real(Zin_fast) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
            Zin_dense_im_err_ohm = dense_row.Zin_im_err_ohm,
            Zin_fast_im_err_ohm = abs(imag(Zin_fast) - imag(Zin_lu)),
        ))
    end

    return rows
end

function main()
    pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell = make_fixture()
    checkpoints = [100, 300, 600]
    restart = 300
    reltol = 1e-6
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

    println("PMCHW medium GMRES trajectory comparison")
    @printf("  N = %d (2N = %d)\n", num_basis(basis), 2 * num_basis(basis))
    @printf("  restart = %d, reltol = %.1e\n", restart, reltol)
    @printf("  checkpoints = %s\n", join(string.(checkpoints), ", "))

    println("  precomputing dense trajectory...")
    dense_rows = precompute_dense_trajectory(pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell, checkpoints; restart = restart, reltol = reltol)

    results = NamedTuple[]
    for budget_case in budgets
        rows = evaluate_case(
            pmchw,
            basis,
            feed,
            Z_dense,
            rhs,
            Zin_lu,
            dense_rows,
            budget_case.name,
            budget_case.budget,
            checkpoints;
            restart = restart,
            reltol = reltol,
        )
        append!(results, rows)

        @printf("\n  budget = %s\n", budget_case.name)
        for row in rows
            @printf("    iter=%3d rel_I=%.6e dense_res_fast=%.6e gap=%.6fΩ/%.6fΩ gmres_res=%.6e\n",
                row.checkpoint,
                row.rel_I,
                row.dense_res_fast,
                row.Zin_gap_re_ohm,
                row.Zin_gap_im_ohm,
                row.gmres_res_fast,
            )
        end
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_gmres_trajectory_medium.csv")
    CSV.write(csv_path, DataFrame(results))
    println("\n  -> saved $(csv_path)")
end

main()