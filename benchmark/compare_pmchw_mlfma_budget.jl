using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMMoMSuite
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using Random
using Printf
using Dates
using CSV
using DataFrames
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function parse_cli(args)
    preset = isempty(args) ? "small" : lowercase(args[1])
    case_mode = length(args) >= 2 ? lowercase(args[2]) : "representative"
    if preset == "small"
        return (
            label = "small",
            freq = 300e6,
            eps_r = 4.0,
            radius = 0.1,
            n_theta = 4,
            n_phi = 6,
            leaf_size = 0.10,
            gmres_reltol = 1e-4,
            gmres_maxiter = 100,
            case_mode = case_mode,
        )
    elseif preset == "mid"
        return (
            label = "mid",
            freq = 300e6,
            eps_r = 4.0,
            radius = 0.2,
            n_theta = 8,
            n_phi = 12,
            leaf_size = 0.10,
            gmres_reltol = 1e-4,
            gmres_maxiter = 100,
            case_mode = case_mode,
        )
    elseif preset == "medium"
        return (
            label = "medium",
            freq = 300e6,
            eps_r = 4.0,
            radius = 0.5,
            n_theta = 10,
            n_phi = 20,
            leaf_size = 0.10,
            gmres_reltol = 1e-4,
            gmres_maxiter = 100,
            case_mode = case_mode,
        )
    end

    error("unsupported preset: $preset (expected: small, mid, or medium)")
end

function make_fixture(; freq, eps_r, radius, n_theta, n_phi)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r)
    feed = DeltaGapSource(freq, [1], 1.0 + 0im)
    return pmchw, basis, feed
end

function build_budget_cases(case_mode = "representative")
    representative = [
        (
            name = "default",
            budget = PMCHWMLFMAErrorBudget(Float64),
            note = "heuristic baseline",
        ),
        (
            name = "loose_near",
            budget = PMCHWMLFMAErrorBudget(Float64;
                near_range_scale = 4.0,
                min_near_range = 4,
                max_near_range = 16,
            ),
            note = "smaller near span",
        ),
        (
            name = "fixed_leaf_0p04_nr9",
            budget = PMCHWMLFMAErrorBudget(Float64;
                fixed_leaf_size_eff = 0.04,
                fixed_near_range = 9,
            ),
            note = @sprintf("fixed leaf %.2f, fixed nr=9", 0.04),
        ),
    ]

    if case_mode == "representative"
        return representative
    elseif case_mode == "full"
        return vcat(
            representative[1:2],
            [
                (
                    name = "tight_near",
                    budget = PMCHWMLFMAErrorBudget(Float64;
                        near_range_scale = 12.0,
                        min_near_range = 12,
                        max_near_range = 32,
                    ),
                    note = "larger near span",
                ),
                (
                    name = "tight_near_L4",
                    budget = PMCHWMLFMAErrorBudget(Float64;
                        near_range_scale = 12.0,
                        min_near_range = 12,
                        max_near_range = 32,
                        L_min = 4,
                    ),
                    note = "larger near span + L_min=4",
                ),
            ],
            representative[3:3],
        )
    end

    error("unsupported case_mode: $case_mode (expected: representative or full)")
end

function relcorr(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, corr
end

function solve_budget_case(pmchw, basis, feed, Z_dense, V_rhs, budget_case; leaf_size, gmres_reltol, gmres_maxiter, y_ref, Zin_ref)
    t_build = @elapsed op = PMCHWMLFMAOperator(pmchw, basis, leaf_size; budget = budget_case.budget)

    y_mlfma = zeros(ComplexF64, length(V_rhs))
    t_matvec = @elapsed mul!(y_mlfma, op, y_ref)
    rel_matvec, corr_matvec = relcorr(y_mlfma, Z_dense * y_ref)

    shell = PMCHWBlockOperator(
        pmchw,
        basis,
        MatrixFreePMCHWBackend(op);
        block_source = Z_dense,
    )

    t_solve = @elapsed coeffs_sf, hist = gmres(
        strong_form(shell),
        strong_form_rhs(shell, V_rhs);
        reltol = gmres_reltol,
        maxiter = gmres_maxiter,
        log = true,
    )
    I_budget = recover_trial_coefficients(shell, coeffs_sf)
    Zin_budget = input_impedance(pmchw, feed, I_budget, basis)
    resnorm = hist.data[:resnorm][end]

    return (
        case = budget_case.name,
        note = budget_case.note,
        leaf_size_eff = op.leaf_size_eff,
        near_range = op.near_range,
        L_min = op.budget.L_min,
        nLevels0 = op.octree0.nLevels,
        nLevels1 = op.octree1.nLevels,
        nnz_near = nnz(op.Z_near),
        near_density = nnz(op.Z_near) / (length(V_rhs)^2),
        rel_matvec = rel_matvec,
        corr_matvec = corr_matvec,
        Zin_re = real(Zin_budget),
        Zin_im = imag(Zin_budget),
        Zin_re_err_pct = abs(real(Zin_budget) - real(Zin_ref)) / (max(abs(real(Zin_ref)), 1.0) + 1e-30) * 100,
        Zin_im_err_ohm = abs(imag(Zin_budget) - imag(Zin_ref)),
        resnorm = resnorm,
        build_s = t_build,
        matvec_s = t_matvec,
        solve_s = t_solve,
    )
end

function run_budget_sweep(; label, freq, eps_r, radius, n_theta, n_phi, leaf_size, gmres_reltol, gmres_maxiter, case_mode)
    pmchw, basis, feed = make_fixture(; freq = freq, eps_r = eps_r, radius = radius, n_theta = n_theta, n_phi = n_phi)
    N = num_basis(basis)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    V_rhs = excitation_vector(pmchw, feed, basis)
    I_ref = Z_dense \ V_rhs
    Zin_ref = input_impedance(pmchw, feed, I_ref, basis)

    Random.seed!(42)
    x_probe = randn(ComplexF64, 2N)
    x_probe ./= norm(x_probe)
    y_direct_probe = Z_dense * x_probe

    cases = build_budget_cases(case_mode)
    results = NamedTuple[]

    println("PMCHW MLFMA budget sweep")
    @printf("  preset = %s\n", label)
    @printf("  N = %d (2N = %d)\n", N, 2N)
    @printf("  freq = %.3f MHz, eps_r = %.3f, radius = %.3f m\n", freq / 1e6, eps_r, radius)
    @printf("  case mode = %s\n", case_mode)
    @printf("  direct reference Zin = %+.6e + j(%+.6e)\n", real(Zin_ref), imag(Zin_ref))
    @printf("  gmres config = reltol %.1e, maxiter %d\n", gmres_reltol, gmres_maxiter)

    for budget_case in cases
        result = solve_budget_case(
            pmchw,
            basis,
            feed,
            Z_dense,
            V_rhs,
            budget_case;
            leaf_size = leaf_size,
            gmres_reltol = gmres_reltol,
            gmres_maxiter = gmres_maxiter,
            y_ref = x_probe,
            Zin_ref = Zin_ref,
        )
        push!(results, result)

        @printf("\n  case = %s  (%s)\n", result.case, result.note)
        @printf("    leaf_size_eff = %.6f, near_range = %d, L_min = %d\n", result.leaf_size_eff, result.near_range, result.L_min)
        @printf("    nLevels = (%d, %d), nnz_near = %d, near_density = %.4f\n", result.nLevels0, result.nLevels1, result.nnz_near, result.near_density)
        @printf("    matvec rel = %.6e, corr = %.6f\n", result.rel_matvec, result.corr_matvec)
        @printf("    Zin = %+.6e + j(%+.6e)\n", result.Zin_re, result.Zin_im)
        @printf("    Zin re err = %.3f %%, im err = %.6f Ω\n", result.Zin_re_err_pct, result.Zin_im_err_ohm)
        @printf("    resnorm = %.6e, build = %.3fs, matvec = %.3fs, solve = %.3fs\n", result.resnorm, result.build_s, result.matvec_s, result.solve_s)
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_MLFMA_budget_sweep_$(label).csv")
    CSV.write(csv_path, DataFrame(results))
    println("\n  -> saved $(csv_path)")
end

config = parse_cli(ARGS)
run_budget_sweep(; config...)