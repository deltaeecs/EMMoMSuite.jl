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
# ─────────────────────────────────────────────────────────────────────────────
# NearField.jl 覆盖率
# ─────────────────────────────────────────────────────────────────────────────
@testset "NearField" begin
    freq = 300e6
    set_frequency!(freq)

    # 最小网格：2 个三角形共享一条边
    nodes    = [0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
    elements = [1 2; 2 4; 3 3]
    tags     = [1, 1]
    mesh     = TriangleMesh(2, nodes, elements, tags)
    basis    = RWGBasis(mesh)
    N        = num_basis(basis)

    I_coeffs = ones(ComplexF64, N)

    # 观测点：距离天线 1m 处
    FT = Float64
    points = SVector{3,FT}[SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.5)]

    E_field = calculate_near_field(points, basis, I_coeffs)

    @test length(E_field) == 2
    @test eltype(E_field) <: SVector{3}
    # 场应为有限值（非 NaN/Inf）
    @test all(all(isfinite, e) for e in E_field)
end

# ─────────────────────────────────────────────────────────────────────────────
# MieSeries.jl 覆盖率
# ─────────────────────────────────────────────────────────────────────────────
@testset "MieSeries" begin
    radius = 0.3      # 30 cm 球
    freq   = 300e6    # 300 MHz

    # 背散射 (theta=pi)
    theta_range = [0.0, pi/2, pi]
    rcs = calculate_mie_rcs_pec_sphere(radius, freq, theta_range)

    @test length(rcs) == 3
    @test all(isfinite, rcs)
    @test all(rcs .>= 0)

    # 后向散射近似: σ_back ≈ π r² for ka ≫ 1 (光学区)
    # 这里 ka ≈ 1.88，处于谐振区，但 RCS 应在 πr² 量级
    c0   = 299792458.0
    k    = 2π * freq / c0
    ka   = k * radius  # ≈ 1.88
    @test rcs[end] > 0  # 后向 RCS > 0
end