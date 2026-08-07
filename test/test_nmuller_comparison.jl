using Test
using EMMoMSuite
using LinearAlgebra
using IterativeSolvers
using Random

function make_nmuller_comparison_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 4, n_phi = 6)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    return pmchw, nmuller, basis
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; reltol = 1e-6, maxiter = 200)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    resnorms = hist.data[:resnorm]
    return (
        x_lu = x_lu,
        x_gmres = x_gmres,
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        resnorm = resnorms[end],
        niters = length(resnorms),
    )
end

@testset "PMCHW vs N-Muller dense comparison" begin
    pmchw, nmuller, basis = make_nmuller_comparison_fixture()
    N = num_basis(basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)

    Random.seed!(42)
    rhs = ComplexF64.(randn(2N), randn(2N))
    rhs ./= norm(rhs)

    cond_pmchw = cond(Z_pmchw)
    cond_nmuller = cond(Z_nmuller)
    gmres_pmchw = solve_gmres_against_lu(Z_pmchw, rhs)
    gmres_nmuller = solve_gmres_against_lu(Z_nmuller, rhs)

    @info "PMCHW vs N-Muller dense comparison" N cond_pmchw cond_nmuller pmchw_iters=gmres_pmchw.niters nmuller_iters=gmres_nmuller.niters pmchw_rel_err=gmres_pmchw.rel_err nmuller_rel_err=gmres_nmuller.rel_err pmchw_resnorm=gmres_pmchw.resnorm nmuller_resnorm=gmres_nmuller.resnorm

    @test isfinite(cond_pmchw)
    @test isfinite(cond_nmuller)
    @test cond_pmchw > 0
    @test cond_nmuller > 0

    @test gmres_pmchw.niters > 0
    @test gmres_nmuller.niters > 0
    @test cond_nmuller < cond_pmchw / 10
    @test gmres_nmuller.rel_err < gmres_pmchw.rel_err
    @test gmres_nmuller.resnorm < gmres_pmchw.resnorm
end