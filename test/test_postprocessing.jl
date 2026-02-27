using Test
using EMSuite
using EMSuite.PostProcessing
using EMSuite.Geometry
using EMSuite.Utilities.Parameters
using EMSuite.BasisFunctions
using StaticArrays
using LinearAlgebra

@testset "PostProcessing" begin
    # Setup parameters
    freq = 300e6 # 300 MHz
    set_frequency!(freq)
    k0 = get_k0()
    
    # Create a proper mesh with 2 triangles sharing an edge
    nodes = [
        0.0 1.0 0.0 1.0;
        0.0 0.0 1.0 1.0;
        0.0 0.0 0.0 0.0
    ]
    elements = [1 2; 2 4; 3 3]
    tags = [1, 1]
    mesh = TriangleMesh(2, nodes, elements, tags)
    basis = RWGBasis(mesh)
    
    # Mock coefficients (one basis function)
    ICoeff = [1.0 + 0.0im]
    
    # Test geoElectricJCal
    Jtris = geoElectricJCal(ICoeff, basis)
    @test size(Jtris, 1) == 3  # 3D vector
    @test size(Jtris, 2) == num_elements(mesh)
    # Check if J is non-zero
    @test norm(Jtris) > 0
    
    # Test RCS
    θs = [0.0, pi/2]
    ϕs = [0.0, pi/2]
    
    RCS_data, RCS_total, RCS_dB = radarCrossSection(θs, ϕs, ICoeff, basis)
    
    @test size(RCS_data) == (2, 2, 2)
    @test size(RCS_total) == (2, 2)
    
    # Test FarField
    source = nothing
    
    farE_data = farField(θs, ϕs, ICoeff, basis, source)
    @test size(farE_data) == (2, 2, 2)

end
