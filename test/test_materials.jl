using Test
using EMSuite
using EMSuite.CoreModule.Materials
using EMSuite.CoreModule.Constants

@testset "Materials" begin
    @testset "PEC" begin
        m = PEC()
        @test permittivity(m) == Inf
        @test permeability(m) == Constants.mu0
    end

    @testset "Dielectric" begin
        # Vacuum
        d0 = Dielectric(1.0, 1.0)
        @test permittivity(d0) ≈ Constants.eps0
        @test permeability(d0) ≈ Constants.mu0
        @test impedance(d0) ≈ Constants.eta0

        # Custom material
        er = 4.0
        ur = 2.0
        d = Dielectric(er, ur)
        @test permittivity(d) ≈ er * Constants.eps0
        @test permeability(d) ≈ ur * Constants.mu0
        @test impedance(d) ≈ sqrt(ur/er) * Constants.eta0
    end
end
