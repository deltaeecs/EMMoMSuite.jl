using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMMoMSuite
using LinearAlgebra
using IterativeSolvers
using Random
using Printf
using CSV
using DataFrames

function make_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    return pmchw, nmuller, basis
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; reltol = 1e-5, maxiter = 250)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    resnorms = hist.data[:resnorm]
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        resnorm = resnorms[end],
        niters = length(resnorms),
    )
end

function main()
    pmchw, nmuller, basis = make_fixture()
    N = num_basis(basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)
    checkpoints = [50, 100, 150, 200, 250]

    Random.seed!(42)
    rhs = ComplexF64.(randn(2N), randn(2N))
    rhs ./= norm(rhs)

    println("PMCHW vs N-Muller medium dense GMRES trajectory")
    @printf("  N = %d (2N = %d)\n", N, 2N)
    @printf("  checkpoints = %s\n", join(string.(checkpoints), ", "))

    rows = NamedTuple[]
    for checkpoint in checkpoints
        pmchw_row = solve_gmres_against_lu(Z_pmchw, rhs; reltol = 1e-5, maxiter = checkpoint)
        nmuller_row = solve_gmres_against_lu(Z_nmuller, rhs; reltol = 1e-5, maxiter = checkpoint)

        push!(rows, (
            checkpoint = checkpoint,
            pmchw_rel_err = pmchw_row.rel_err,
            pmchw_resnorm = pmchw_row.resnorm,
            pmchw_niters = pmchw_row.niters,
            nmuller_rel_err = nmuller_row.rel_err,
            nmuller_resnorm = nmuller_row.resnorm,
            nmuller_niters = nmuller_row.niters,
            rel_err_ratio = nmuller_row.rel_err / (pmchw_row.rel_err + 1e-30),
            resnorm_ratio = nmuller_row.resnorm / (pmchw_row.resnorm + 1e-30),
        ))

        @printf("  iter=%3d  PMCHW rel=%.6e res=%.6e | NMuller rel=%.6e res=%.6e\n",
            checkpoint,
            pmchw_row.rel_err,
            pmchw_row.resnorm,
            nmuller_row.rel_err,
            nmuller_row.resnorm,
        )
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_NMuller_gmres_trajectory_medium.csv")
    CSV.write(csv_path, DataFrame(rows))
    println("\n  -> saved $(csv_path)")
end

main()