using Test

using EMMoMSuite

@testset "EMMoMSuite package import" begin
    @test isdefined(EMMoMSuite, :set_frequency!)
    @test isdefined(EMMoMSuite, :PlaneWave)
    @test isdefined(EMMoMSuite, :EMMoMSuiteConfig)
end