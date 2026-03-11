using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using IterativeSolvers
using Printf
using CSV
using DataFrames

function make_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    return pmchw, nmuller, basis, source
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; reltol = 1e-5, maxiter = 250)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = norm(rhs - A * x_gmres) / (norm(rhs) + 1e-30),
        niters = length(hist.data[:resnorm]),
    )
end

function main()
    pmchw, nmuller, basis, source = make_fixture()
    N = num_basis(basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)
    rhs_pmchw = excitation_vector(pmchw, source, basis)
    rhs_nmuller = excitation_vector(nmuller, source, basis)
    checkpoints = [50, 100, 150, 200, 250]

    println("PMCHW vs N-Muller medium plane-wave dense GMRES trajectory")
    @printf("  N = %d (2N = %d)\n", N, 2N)
    @printf("  checkpoints = %s\n", join(string.(checkpoints), ", "))
    @printf("  rhs norms    = PMCHW %.6e | NMuller %.6e\n", norm(rhs_pmchw), norm(rhs_nmuller))

    rows = NamedTuple[]
    for checkpoint in checkpoints
        pmchw_row = solve_gmres_against_lu(Z_pmchw, rhs_pmchw; reltol = 1e-5, maxiter = checkpoint)
        nmuller_row = solve_gmres_against_lu(Z_nmuller, rhs_nmuller; reltol = 1e-5, maxiter = checkpoint)

        push!(rows, (
            checkpoint = checkpoint,
            pmchw_rel_err = pmchw_row.rel_err,
            pmchw_rel_res = pmchw_row.rel_res,
            pmchw_niters = pmchw_row.niters,
            nmuller_rel_err = nmuller_row.rel_err,
            nmuller_rel_res = nmuller_row.rel_res,
            nmuller_niters = nmuller_row.niters,
            rel_err_ratio = nmuller_row.rel_err / (pmchw_row.rel_err + 1e-30),
            rel_res_ratio = nmuller_row.rel_res / (pmchw_row.rel_res + 1e-30),
        ))

        @printf("  iter=%3d  PMCHW rel=%.6e rres=%.6e iters=%3d | NMuller rel=%.6e rres=%.6e iters=%3d\n",
            checkpoint,
            pmchw_row.rel_err,
            pmchw_row.rel_res,
            pmchw_row.niters,
            nmuller_row.rel_err,
            nmuller_row.rel_res,
            nmuller_row.niters,
        )
    end

    result_dir = joinpath(@__DIR__, "..", "test_results", "accuracy")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "PMCHW_NMuller_planewave_gmres_trajectory_medium.csv")
    CSV.write(csv_path, DataFrame(rows))
    println("\n  -> saved $(csv_path)")
end

main()