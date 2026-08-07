using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMMoMSuite
using LinearAlgebra
using IterativeSolvers
using Printf
using CSV
using DataFrames

function make_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    shell = PMCHWBlockOperator(pmchw, basis)
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    rhs = excitation_vector(pmchw, source, basis)
    return pmchw, basis, shell, rhs
end

function solve_weak_against_lu(shell, rhs; reltol = 1e-5, maxiter = 250)
    A = weak_form(shell)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = norm(rhs - A * x_gmres) / (norm(rhs) + 1e-30),
        niters = length(hist.data[:resnorm]),
    )
end

function solve_strong_against_lu(shell, rhs; reltol = 1e-5, maxiter = 250)
    A = strong_form(shell)
    rhs_sf = strong_form_rhs(shell, rhs)
    coeffs_lu = A \ rhs_sf
    coeffs_gmres, hist = gmres(A, rhs_sf; reltol = reltol, maxiter = maxiter, log = true)
    x_lu = recover_trial_coefficients(shell, coeffs_lu)
    x_gmres = recover_trial_coefficients(shell, coeffs_gmres)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = norm(rhs - weak_form(shell) * x_gmres) / (norm(rhs) + 1e-30),
        niters = length(hist.data[:resnorm]),
    )
end

function main()
    pmchw, basis, shell, rhs = make_fixture()
    N = num_basis(basis)
    checkpoints = [50, 100, 150, 200, 250]

    println("PMCHW medium plane-wave dense weak/strong GMRES trajectory")
    @printf("  N = %d (2N = %d)\n", N, 2N)
    @printf("  checkpoints = %s\n", join(string.(checkpoints), ", "))
    @printf("  rhs norm     = %.6e\n", norm(rhs))

    rows = NamedTuple[]
    for checkpoint in checkpoints
        weak_row = solve_weak_against_lu(shell, rhs; reltol = 1e-5, maxiter = checkpoint)
        strong_row = solve_strong_against_lu(shell, rhs; reltol = 1e-5, maxiter = checkpoint)

        push!(rows, (
            checkpoint = checkpoint,
            weak_rel_err = weak_row.rel_err,
            weak_rel_res = weak_row.rel_res,
            weak_niters = weak_row.niters,
            strong_rel_err = strong_row.rel_err,
            strong_rel_res = strong_row.rel_res,
            strong_niters = strong_row.niters,
            rel_err_ratio = strong_row.rel_err / (weak_row.rel_err + 1e-30),
            rel_res_ratio = strong_row.rel_res / (weak_row.rel_res + 1e-30),
        ))

        @printf("  iter=%3d  weak rel=%.6e rres=%.6e iters=%3d | strong rel=%.6e rres=%.6e iters=%3d\n",
            checkpoint,
            weak_row.rel_err,
            weak_row.rel_res,
            weak_row.niters,
            strong_row.rel_err,
            strong_row.rel_res,
            strong_row.niters,
        )
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_gate_s_planewave_trajectory_medium.csv")
    CSV.write(csv_path, DataFrame(rows))
    println("\n  -> saved $(csv_path)")
end

main()