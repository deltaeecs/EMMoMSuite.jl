using Test
using EMSuite
using EMSuite.IntegralEquations
using EMSuite.PostProcessing
using EMSuite.FastAlgorithms

@testset "Legacy Parity Constants" begin
    # Constants from MoM_Basics/src/ParametersSet.jl
    # JKeta_div_16pi = 1im * k0 * eta0 / (16pi)
    # div4pi = 1 / (4pi)
    # minus_jk_div_16pi2 = -1im * k0 / (16pi^2)

    freq = 300e6
    c0 = 299792458.0
    mu0 = 4 * pi * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2 * pi * freq / c0
    eta = sqrt(mu0 / eps0)

    # 1) EFIE factor parity
    expected_efie_factor = im * k * eta / (16 * pi)
    efie = EFIE(freq)
    @test efie.factor ≈ expected_efie_factor atol = 1e-10

    # 2) Far-field factor parity (documented expectation)
    expected_ff_factor = -im * k * eta / (4 * pi)
    println("Expected FarField Factor: ", expected_ff_factor)

    # 3) Translation factor parity (documented expectation)
    expected_trans_factor = -im * k / (16 * pi^2)
    println("Expected Translation Factor: ", expected_trans_factor)
end
