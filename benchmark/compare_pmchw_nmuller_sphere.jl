using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using IterativeSolvers
using Random
using Printf

function parse_cli_config(args)
    preset = isempty(args) ? "small" : lowercase(args[1])
    if preset == "small"
        return (
            label = "small",
            freq = 120e6,
            eps_r = 4.0,
            mu_r = 1.0,
            radius = 0.1,
            n_theta = 4,
            n_phi = 6,
            reltol = 1e-6,
            maxiter = 200,
        )
    elseif preset == "medium"
        return (
            label = "medium",
            freq = 120e6,
            eps_r = 4.0,
            mu_r = 1.0,
            radius = 0.1,
            n_theta = 6,
            n_phi = 10,
            reltol = 1e-5,
            maxiter = 250,
        )
    end

    error("unsupported preset: $preset (expected: small or medium)")
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; reltol = 1e-6, maxiter = 200)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    resnorms = hist.data[:resnorm]
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        resnorm = resnorms[end],
        niters = length(resnorms),
    )
end

function run_comparison(; label = "custom", freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 4, n_phi = 6, reltol = 1e-6, maxiter = 200)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    N = num_basis(basis)

    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)

    Random.seed!(42)
    rhs = ComplexF64.(randn(2N), randn(2N))
    rhs ./= norm(rhs)

    cond_pmchw = cond(Z_pmchw)
    cond_nmuller = cond(Z_nmuller)
    gmres_pmchw = solve_gmres_against_lu(Z_pmchw, rhs; reltol = reltol, maxiter = maxiter)
    gmres_nmuller = solve_gmres_against_lu(Z_nmuller, rhs; reltol = reltol, maxiter = maxiter)

    println("PMCHW vs N-Muller sphere comparison")
    @printf("  preset = %s\n", label)
    @printf("  N = %d (2N = %d)\n", N, 2N)
    @printf("  freq = %.3f MHz, eps_r = %.3f, radius = %.3f m\n", freq / 1e6, real(eps_r), radius)
    @printf("  gmres config = reltol %.1e, maxiter %d\n", reltol, maxiter)
    @printf("  cond(PMCHW)   = %.6e\n", cond_pmchw)
    @printf("  cond(NMuller) = %.6e\n", cond_nmuller)
    @printf("  cond ratio    = %.6e\n", cond_nmuller / cond_pmchw)
    @printf("  GMRES PMCHW   : iters=%d, res=%.6e, rel_vs_LU=%.6e\n", gmres_pmchw.niters, gmres_pmchw.resnorm, gmres_pmchw.rel_err)
    @printf("  GMRES NMuller : iters=%d, res=%.6e, rel_vs_LU=%.6e\n", gmres_nmuller.niters, gmres_nmuller.resnorm, gmres_nmuller.rel_err)
end

config = parse_cli_config(ARGS)
run_comparison(; config...)