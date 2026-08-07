using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMMoMSuite
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using Printf
using CSV
using DataFrames
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function parse_cli(args)
    preset = isempty(args) ? "medium" : lowercase(args[1])
    if preset != "medium"
        error("unsupported preset: $preset (only medium is supported for the long-Krylov budget comparison)")
    end
    return (
        label = "medium",
        freq = 300e6,
        eps_r = 4.0,
        radius = 0.5,
        n_theta = 10,
        n_phi = 20,
        leaf_size = 0.10,
    )
end

function make_fixture(; freq, eps_r, radius, n_theta, n_phi)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r)
    feed = DeltaGapSource(freq, [1], 1.0 + 0im)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    rhs = excitation_vector(pmchw, feed, basis)
    Zin_lu = input_impedance(pmchw, feed, Z_dense \ rhs, basis)
    dense_shell = PMCHWBlockOperator(pmchw, basis)
    return pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell
end

function representative_budgets()
    return [
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
        (
            name = "fixed_leaf_0p04_nr9",
            budget = PMCHWMLFMAErrorBudget(Float64;
                fixed_leaf_size_eff = 0.04,
                fixed_near_range = 9,
            ),
        ),
    ]
end

function krylov_configs()
    return [
        (
            name = "short",
            restart = 100,
            maxiter = 100,
            reltol = 1e-4,
        ),
        (
            name = "long",
            restart = 300,
            maxiter = 600,
            reltol = 1e-6,
        ),
    ]
end

function solve_strong(shell, rhs; restart, maxiter, reltol)
    coeffs_sf, hist = gmres(
        strong_form(shell),
        strong_form_rhs(shell, rhs);
        restart = restart,
        maxiter = maxiter,
        reltol = reltol,
        log = true,
    )
    return recover_trial_coefficients(shell, coeffs_sf), hist.data[:resnorm][end], length(hist.data[:resnorm])
end

function run_comparison(; label, freq, eps_r, radius, n_theta, n_phi, leaf_size)
    pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell = make_fixture(; freq = freq, eps_r = eps_r, radius = radius, n_theta = n_theta, n_phi = n_phi)
    budgets = representative_budgets()
    configs = krylov_configs()
    results = NamedTuple[]

    println("PMCHW MLFMA budget long-Krylov comparison")
    @printf("  preset = %s\n", label)
    @printf("  N = %d (2N = %d)\n", num_basis(basis), 2 * num_basis(basis))
    @printf("  direct reference Zin = %+.6e + j(%+.6e)\n", real(Zin_lu), imag(Zin_lu))

    dense_by_config = Dict{String,NamedTuple}()
    for config in configs
        t_dense = @elapsed I_dense, res_dense, niters_dense = solve_strong(dense_shell, rhs; restart = config.restart, maxiter = config.maxiter, reltol = config.reltol)
        Zin_dense = input_impedance(pmchw, feed, I_dense, basis)
        dense_by_config[config.name] = (
            Zin = Zin_dense,
            resnorm = res_dense,
            niters = niters_dense,
            solve_s = t_dense,
        )
        @printf("\n  dense config = %s  restart=%d maxiter=%d reltol=%.1e\n", config.name, config.restart, config.maxiter, config.reltol)
        @printf("    Zin = %+.6e + j(%+.6e)\n", real(Zin_dense), imag(Zin_dense))
        @printf("    re err vs LU = %.3f %%, im err vs LU = %.6f Ω, resnorm = %.6e\n",
            abs(real(Zin_dense) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
            abs(imag(Zin_dense) - imag(Zin_lu)),
            res_dense,
        )
    end

    for config in configs
        dense_result = dense_by_config[config.name]
        for budget_case in budgets
            t_build = @elapsed op = PMCHWMLFMAOperator(pmchw, basis, leaf_size; budget = budget_case.budget)
            shell = PMCHWBlockOperator(pmchw, basis, MatrixFreePMCHWBackend(op); block_source = Z_dense)
            t_solve = @elapsed I_mlfma, res_mlfma, niters_mlfma = solve_strong(shell, rhs; restart = config.restart, maxiter = config.maxiter, reltol = config.reltol)
            Zin_mlfma = input_impedance(pmchw, feed, I_mlfma, basis)

            result = (
                config = config.name,
                budget = budget_case.name,
                restart = config.restart,
                maxiter = config.maxiter,
                reltol = config.reltol,
                leaf_size_eff = op.leaf_size_eff,
                near_range = op.near_range,
                near_density = nnz(op.Z_near) / (length(rhs)^2),
                Zin_dense_re = real(dense_result.Zin),
                Zin_dense_im = imag(dense_result.Zin),
                Zin_mlfma_re = real(Zin_mlfma),
                Zin_mlfma_im = imag(Zin_mlfma),
                Zin_dense_re_err_pct = abs(real(dense_result.Zin) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
                Zin_mlfma_re_err_pct = abs(real(Zin_mlfma) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
                Zin_dense_im_err_ohm = abs(imag(dense_result.Zin) - imag(Zin_lu)),
                Zin_mlfma_im_err_ohm = abs(imag(Zin_mlfma) - imag(Zin_lu)),
                Zin_gap_re_ohm = abs(real(Zin_mlfma) - real(dense_result.Zin)),
                Zin_gap_im_ohm = abs(imag(Zin_mlfma) - imag(dense_result.Zin)),
                res_dense = dense_result.resnorm,
                res_mlfma = res_mlfma,
                niters_dense = dense_result.niters,
                niters_mlfma = niters_mlfma,
                build_s = t_build,
                solve_s = t_solve,
            )
            push!(results, result)

            @printf("\n  config = %s, budget = %s\n", config.name, budget_case.name)
            @printf("    leaf_size_eff = %.6f, near_range = %d, near_density = %.4f\n", result.leaf_size_eff, result.near_range, result.near_density)
            @printf("    Zin_dense = %+.6e + j(%+.6e)\n", result.Zin_dense_re, result.Zin_dense_im)
            @printf("    Zin_mlfma = %+.6e + j(%+.6e)\n", result.Zin_mlfma_re, result.Zin_mlfma_im)
            @printf("    dense re/im err = %.3f %% / %.6f Ω\n", result.Zin_dense_re_err_pct, result.Zin_dense_im_err_ohm)
            @printf("    mlfma re/im err = %.3f %% / %.6f Ω\n", result.Zin_mlfma_re_err_pct, result.Zin_mlfma_im_err_ohm)
            @printf("    mlfma-dense gap = %.6f Ω / %.6f Ω, res = %.6e vs %.6e\n", result.Zin_gap_re_ohm, result.Zin_gap_im_ohm, result.res_dense, result.res_mlfma)
            @printf("    build = %.3fs, solve = %.3fs\n", result.build_s, result.solve_s)
        end
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_MLFMA_budget_krylov_medium.csv")
    CSV.write(csv_path, DataFrame(results))
    println("\n  -> saved $(csv_path)")
end

config = parse_cli(ARGS)
run_comparison(; config...)