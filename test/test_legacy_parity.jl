using Testusing Test





































end    println("Expected Translation Factor: ", expected_trans_factor)    expected_trans_factor = -im * k / (16 * π^2)    # Legacy: -jk / (16 * pi^2)    # 3. Translation Factor    println("Expected FarField Factor: ", expected_ff_factor)    expected_ff_factor = -im * k * eta / (4 * π)    # Legacy: -jk * eta / (4 * pi)    # 2. FarField Factor        @test efie.factor ≈ expected_efie_factor atol=1e-10    efie = EFIE(freq)        expected_efie_factor = im * k * eta / (16 * π)    # Legacy: jk * eta / (16 * pi)    # 1. EFIE Factor    eta = sqrt(mu0 / eps0)    k = 2π * freq / c0    eps0 = 1.0 / (c0^2 * mu0)    mu0 = 4π * 1e-7    c0 = 299792458.0    freq = 300e6    # mjKdiv16π² = -Params.JK_0/(4π)^2    # div4π = 1/4π    # JKηdiv16π = 1im*K_0*η_0/16π    # Constants from MoM_Basics/src/ParametersSet.jl@testset "Legacy Parity Constants" beginusing EMSuite.FastAlgorithmsusing EMSuite.PostProcessingusing EMSuite.IntegralEquationsusing EMSuiteusing EMSuite
using EMSuite.IntegralEquations
using EMSuite.PostProcessing
using EMSuite.FastAlgorithms

@testset "Legacy Parity Constants" begin
    # Constants from MoM_Basics/src/ParametersSet.jl
    # JKηdiv16π = 1im*K_0*η_0/16π
    # div4π = 1/4π
    # mjKdiv16π² = -Params.JK_0/(4π)^2

    freq = 300e6
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)

    # 1. EFIE Factor
    # Legacy: jk * eta / (16 * pi)
    expected_efie_factor = im * k * eta / (16 * π)
    
    efie = EFIE(freq)
    @test efie.factor ≈ expected_efie_factor atol=1e-10
    
    # 2. FarField Factor
    # Legacy: -jk * eta / (4 * pi)
    # Note: FarField.jl uses local variable, so we can't check it directly from a struct.
    # But we can check the code or output.
    # Here we just document the expectation.
    expected_ff_factor = -im * k * eta / (4 * π)
    println("Expected FarField Factor: ", expected_ff_factor)

    # 3. Translation Factor
    # Legacy: -jk / (16 * pi^2)
    expected_trans_factor = -im * k / (16 * π^2)
    
    # We need to access the translation factor from MLFMA.
    # It's usually in the Translation operator or function.
    # Let's check Translation.jl to see if it's exposed.
    # It's a local constant in `translation_matrix!`.
    # We can't test it directly without exposing it.
    println("Expected Translation Factor: ", expected_trans_factor)

end
