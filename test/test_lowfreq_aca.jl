# test/test_lowfreq_aca.jl - low-frequency ACA complement path (TDD)
#
# At 30 MHz (lambda=10 m, sphere ~0.1 lambda), standard plane-wave MLFMA
# suffers low-frequency breakdown (docs/theory/fast_algorithms.md Sec 3/6).
# ACA is kernel-agnostic and should stay accurate: MatVec and GMRES solution
# must match the dense reference within 1e-2.

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using LinearAlgebra
using IterativeSolvers

@testset "low-frequency ACA complement" begin
    freq = 30e6
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    lambda = 299792458.0 / freq

    Z = assemble_impedance_matrix(efie, basis)
    # 0.03λ 叶层：球体包围盒内出现多层盒子与非邻对，真正触发低频远块压缩
    op = ACAOperator(efie, basis, 0.03 * lambda; tol = 1e-4, near_range = 1)

    @test size(op) == (N, N)
    @test !isempty(op.blocks)
    x = randn(ComplexF64, N)
    y = op * x
    mv_err = norm(y - Z * x) / norm(Z * x)
    @test mv_err < 1e-2

    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    I_aca = gmres(op, V; abstol = 1e-6, reltol = 1e-8, maxiter = 300, restart = 50)
    I_dir = Z \ V
    @test norm(I_aca - I_dir) / norm(I_dir) < 1e-2
end
