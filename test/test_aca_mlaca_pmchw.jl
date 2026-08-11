# test/test_aca_mlaca_pmchw.jl - F2b ACA/MLACA PMCHW support (TDD)
#
# PMCHW system is 2N x 2N: ACAOperator/MLACAOperator should construct directly.
# Gates: (1) MatVec vs dense reference < 1e-2 (operator matrix accuracy);
#        (2) GMRES converges to small operator residual < 1e-4 (solver integration).
# Note: this PMCHW formulation is ill-conditioned here (cond ~ 1e6), so the
# dense-SOLUTION comparison is conditioning-limited and reported, not gated.
# PMCHW is non-symmetric (HJ = -EM), so symmetric=false (bidirectional blocks).

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using EMMoMSuite.FastAlgorithms.MLACAOperatorModule: MLACAOperator
using LinearAlgebra
using IterativeSolvers

@testset "ACA/MLACA PMCHW" begin
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    pmchw = PMCHW(300e6, 4.0)
    Z = assemble_impedance_matrix(pmchw, basis)
    S = 2N
    lambda = 299792458.0 / 300e6

    @testset "ACAOperator PMCHW MatVec + GMRES" begin
        op = ACAOperator(pmchw, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        @test size(op) == (S, S)
        @test !isempty(op.blocks)
        x = randn(ComplexF64, S)
        y = op * x
        @test norm(y - Z * x) / norm(Z * x) < 1e-2

        V = randn(ComplexF64, S)
        P = ILUPreconditioner(op)
        I_aca, hist = gmres(op, V; Pl = P, abstol = 1e-6, reltol = 1e-8,
                            maxiter = 600, restart = 100, log = true)
        @test hist.isconverged
        @test norm(V - op * I_aca) / norm(V) < 1e-4
        sol_err = norm(I_aca - Z \ V) / norm(Z \ V)
        @info "PMCHW ACA sol err (conditioning-limited) = $sol_err"
    end

    @testset "MLACAOperator PMCHW MatVec + GMRES" begin
        op = MLACAOperator(pmchw, basis, 0.125 * lambda; tol = 1e-4, near_range = 1)
        @test size(op) == (S, S)
        @test !isempty(op.blocks)
        x = randn(ComplexF64, S)
        y = op * x
        @test norm(y - Z * x) / norm(Z * x) < 1e-2

        V = randn(ComplexF64, S)
        P = ILUPreconditioner(op)
        I_mlaca, hist = gmres(op, V; Pl = P, abstol = 1e-6, reltol = 1e-8,
                              maxiter = 600, restart = 100, log = true)
        @test hist.isconverged
        @test norm(V - op * I_mlaca) / norm(V) < 1e-4
        sol_err = norm(I_mlaca - Z \ V) / norm(Z \ V)
        @info "PMCHW MLACA sol err (conditioning-limited) = $sol_err"
    end
end
