"""
test_pmchw_gate_s_mlfma_medium.jl — Phase 15 Gate S 中尺度四路对照回归

目的：在 shell 语义下，把下列四路放到同一中尺度夹具上做正式对照：
  - dense weak
  - dense strong
  - MLFMA weak
  - MLFMA strong

说明：
  - 本回归固定使用 N=540 球面 PMCHW 夹具。
  - 运行成本较高（分钟级），因此暂不并入默认 runtests；作为 Gate S 的专门回归入口保留。
"""

using Test
using EMMoMSuite
using IterativeSolvers
using LinearAlgebra
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_gate_s_mlfma_medium_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    dense_shell = PMCHWBlockOperator(pmchw, basis)
    mlfma_shell = PMCHWBlockOperator(
        pmchw,
        basis,
        MatrixFreePMCHWBackend(PMCHWMLFMAOperator(pmchw, basis, 0.10)),
    )
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    rhs = excitation_vector(pmchw, feed, basis)
    return pmchw, basis, dense_shell, mlfma_shell, feed, rhs
end

function _solve_weak(shell, rhs)
    coeffs, hist = gmres(weak_form(shell), rhs; reltol = 1e-4, maxiter = 100, log = true)
    return coeffs, hist.data[:resnorm][end]
end

function _solve_strong(shell, rhs)
    coeffs_sf, hist = gmres(
        strong_form(shell),
        strong_form_rhs(shell, rhs);
        reltol = 1e-4,
        maxiter = 100,
        log = true,
    )
    return recover_trial_coefficients(shell, coeffs_sf), hist.data[:resnorm][end]
end

@testset "15.S2 Medium full Gate S split" begin
    pmchw, basis, dense_shell, mlfma_shell, feed, rhs = make_gate_s_mlfma_medium_fixture()

    I_dense_weak, res_dense_weak = _solve_weak(dense_shell, rhs)
    I_dense_strong, res_dense_strong = _solve_strong(dense_shell, rhs)
    I_mlfma_weak, res_mlfma_weak = _solve_weak(mlfma_shell, rhs)
    I_mlfma_strong, res_mlfma_strong = _solve_strong(mlfma_shell, rhs)

    Zin_dense_weak = input_impedance(pmchw, feed, I_dense_weak, basis)
    Zin_dense_strong = input_impedance(pmchw, feed, I_dense_strong, basis)
    Zin_mlfma_weak = input_impedance(pmchw, feed, I_mlfma_weak, basis)
    Zin_mlfma_strong = input_impedance(pmchw, feed, I_mlfma_strong, basis)

    relw = norm(I_mlfma_weak - I_dense_weak) / (norm(I_dense_weak) + 1e-30)
    rels = norm(I_mlfma_strong - I_dense_strong) / (norm(I_dense_strong) + 1e-30)
    zgw = abs(Zin_mlfma_weak - Zin_dense_weak) / (abs(Zin_dense_weak) + 1e-30)
    zgs = abs(Zin_mlfma_strong - Zin_dense_strong) / (abs(Zin_dense_strong) + 1e-30)
    resgap_w = abs(res_mlfma_weak - res_dense_weak) / (abs(res_dense_weak) + 1e-30)
    resgap_s = abs(res_mlfma_strong - res_dense_strong) / (abs(res_dense_strong) + 1e-30)

    @info "Medium full Gate S split" relw rels zgw zgs res_dense_weak res_dense_strong res_mlfma_weak res_mlfma_strong resgap_w resgap_s

    @test num_basis(basis) == 540
    @test relw < 1e-2
    @test rels < 1e-2
    @test zgw < 1e-4
    @test zgs < 1e-4
    @test res_dense_strong < res_dense_weak
    @test res_mlfma_strong < res_mlfma_weak
    @test resgap_w < 1e-2
    @test resgap_s < 1e-2
end