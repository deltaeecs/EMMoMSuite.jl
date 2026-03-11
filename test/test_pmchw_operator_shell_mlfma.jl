using Test
using EMSuite
using IterativeSolvers
using LinearAlgebra
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_pmchw_shell_mlfma_fixture(; freq = 300e6, eps_r = 4.0, radius = 0.1, n_theta = 4, n_phi = 6)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r)
    Z = assemble_impedance_matrix(pmchw, basis)
    dense_shell = PMCHWBlockOperator(pmchw, basis, DensePMCHWBackend(Z))
    mlfma_op = PMCHWMLFMAOperator(pmchw, basis, 0.10)
    mlfma_shell = PMCHWBlockOperator(pmchw, basis, MatrixFreePMCHWBackend(mlfma_op); block_source = Z)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    V = excitation_vector(pmchw, feed, basis)
    return pmchw, basis, dense_shell, mlfma_shell, feed, V
end

@testset "PMCHW MLFMA backend shell" begin
    pmchw, basis, dense_shell, mlfma_shell, feed, V = make_pmchw_shell_mlfma_fixture()
    N = num_basis(basis)

    @testset "matrix-free weak form matches wrapped operator" begin
        x = ComplexF64.(randn(2N), randn(2N))
        y_shell = mlfma_shell * x
        y_backend = weak_form(mlfma_shell) * x
        @test y_shell ≈ y_backend rtol=1e-12 atol=1e-12
    end

    @testset "strong form wrapper is executable for MLFMA backend" begin
        A_strong = strong_form(mlfma_shell)
        rhs_strong = strong_form_rhs(mlfma_shell, V)
        coeffs_strong, hist_strong = gmres(A_strong, rhs_strong; reltol = 1e-4, maxiter = 100, log = true)
        I_strong = recover_trial_coefficients(mlfma_shell, coeffs_strong)
        @test length(I_strong) == 2N
        @test isfinite(norm(I_strong))
        @test hist_strong.data[:resnorm][end] < 1e-2
    end

    @testset "MLFMA shell weak/strong stay close to dense shell on small fixture" begin
        I_dense_weak, _ = gmres(weak_form(dense_shell), V; reltol = 1e-4, maxiter = 100, log = true)

        coeffs_dense_strong, _ = gmres(
            strong_form(dense_shell),
            strong_form_rhs(dense_shell, V);
            reltol = 1e-4,
            maxiter = 100,
            log = true,
        )
        I_dense_strong = recover_trial_coefficients(dense_shell, coeffs_dense_strong)

        I_mlfma_weak, _ = gmres(weak_form(mlfma_shell), V; reltol = 1e-4, maxiter = 100, log = true)
        coeffs_mlfma_strong, _ = gmres(
            strong_form(mlfma_shell),
            strong_form_rhs(mlfma_shell, V);
            reltol = 1e-4,
            maxiter = 100,
            log = true,
        )
        I_mlfma_strong = recover_trial_coefficients(mlfma_shell, coeffs_mlfma_strong)

        Zin_dense_weak = input_impedance(pmchw, feed, I_dense_weak, basis)
        Zin_dense_strong = input_impedance(pmchw, feed, I_dense_strong, basis)
        Zin_mlfma_weak = input_impedance(pmchw, feed, I_mlfma_weak, basis)
        Zin_mlfma_strong = input_impedance(pmchw, feed, I_mlfma_strong, basis)

        weak_gap = abs(Zin_mlfma_weak - Zin_dense_weak) / (abs(Zin_dense_weak) + 1e-30)
        strong_gap = abs(Zin_mlfma_strong - Zin_dense_strong) / (abs(Zin_dense_strong) + 1e-30)

        @test weak_gap < 5e-2
        @test strong_gap < 5e-2
    end
end