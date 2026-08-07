using Test
using EMMoMSuite
using IterativeSolvers
using LinearAlgebra

function make_gate_s_medium_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    shell = PMCHWBlockOperator(pmchw, basis)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    V = excitation_vector(pmchw, feed, basis)
    return pmchw, basis, shell, feed, V
end

@testset "15.S1 Dense weak/strong Gate S split" begin
    pmchw, basis, shell, feed, V = make_gate_s_medium_fixture()

    Z = weak_form(shell)
    I_lu = Z \ V
    Zin_lu = input_impedance(pmchw, feed, I_lu, basis)

    I_weak, hist_weak = gmres(weak_form(shell), V; reltol = 1e-4, maxiter = 200, log = true)
    I_strong_coeff, hist_strong = gmres(
        strong_form(shell),
        strong_form_rhs(shell, V);
        reltol = 1e-4,
        maxiter = 200,
        log = true,
    )
    I_strong = recover_trial_coefficients(shell, I_strong_coeff)

    Zin_weak = input_impedance(pmchw, feed, I_weak, basis)
    Zin_strong = input_impedance(pmchw, feed, I_strong, basis)

    weak_err = abs(Zin_weak - Zin_lu) / (abs(Zin_lu) + 1e-30)
    strong_err = abs(Zin_strong - Zin_lu) / (abs(Zin_lu) + 1e-30)
    weak_res = hist_weak.data[:resnorm][end]
    strong_res = hist_strong.data[:resnorm][end]

    @info "Dense Gate S split" Zin_lu Zin_weak Zin_strong weak_err strong_err weak_res strong_res

    @test num_basis(basis) == 540
    @test strong_err < weak_err
    @test strong_res <= weak_res
end