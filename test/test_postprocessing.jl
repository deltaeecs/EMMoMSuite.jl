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

# ─────────────────────────────────────────────────────────────────────────────
# FieldCut.jl
# ─────────────────────────────────────────────────────────────────────────────
@testset "FieldCut" begin
    using EMSuite.PostProcessing: field_cut_line, field_cut_plane

    freq = 300e6
    set_frequency!(freq)

    # 最小 RWG 网格
    nodes    = [0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
    elements = [1 2; 2 4; 3 3]
    mesh     = TriangleMesh(2, nodes, elements, [1,1])
    basis    = RWGBasis(mesh)
    I_coeffs = ones(ComplexF64, num_basis(basis))

    FT = Float64

    # --- field_cut_line ---
    p1 = SVector{3,FT}(2.0, 0.0, 0.0)
    p2 = SVector{3,FT}(2.0, 0.0, 2.0)
    N  = 5

    pts_line, E_line = field_cut_line(p1, p2, N, basis, I_coeffs)

    @test length(pts_line) == N
    @test length(E_line)   == N
    @test eltype(pts_line) <: SVector{3}
    @test eltype(E_line)   <: SVector{3}
    # 第一个点等于 p1，最后一个等于 p2
    @test pts_line[1]   ≈ p1
    @test pts_line[end] ≈ p2
    # 所有场值有限
    @test all(all(isfinite, e) for e in E_line)

    # --- field_cut_plane ---
    origin = SVector{3,FT}(2.0, -1.0, -1.0)
    u_vec  = SVector{3,FT}(0.0, 1.0, 0.0)   # y 方向
    v_vec  = SVector{3,FT}(0.0, 0.0, 1.0)   # z 方向
    Nu, Nv = 3, 4

    pts_plane, E_plane = field_cut_plane(origin, u_vec, v_vec, Nu, Nv, basis, I_coeffs)

    @test size(pts_plane) == (Nu, Nv)
    @test size(E_plane)   == (Nu, Nv)
    @test all(all(isfinite, e) for e in E_plane)
    # 验证角点坐标
    @test pts_plane[1, 1] ≈ origin
    @test pts_plane[end, 1] ≈ origin + u_vec
    @test pts_plane[1, end] ≈ origin + v_vec
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 17.3 — geoVolumeCurrentCal
# ─────────────────────────────────────────────────────────────────────────────
@testset "geoVolumeCurrentCal" begin
    using EMSuite.PostProcessing: geoVolumeCurrentCal

    freq = 300e6
    set_frequency!(freq)

    # 最小 Tet 网格：1×1×1 box，4×4×4 分割 → 6*4³=384 tets
    mesh = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
    FT   = Float64
    CT   = ComplexF64

    # ─── SWG 基 ──────────────────────────────────────────────────────────────
    basis_swg = SWGBasis(mesh)
    N_swg     = num_basis(basis_swg)
    I_swg     = ones(CT, N_swg)
    ε_r       = 2.0 .+ 0.0im
    perms_swg = fill(CT(ε_r), num_elements(mesh))

    J_swg = geoVolumeCurrentCal(I_swg, basis_swg, perms_swg)

    # 返回每个 tet 一个 SVector{3,CT}
    @test length(J_swg) == num_elements(mesh)
    @test eltype(J_swg) <: SVector{3}
    @test all(all(isfinite, j) for j in J_swg)
    # 对称网格上合力矩应接近 0（各向同性均匀 mesh 期望净 J≈0）
    J_mag = norm.(J_swg)
    @test maximum(J_mag) > 0.0   # 场不全为零

    # ─── PWC 基 ──────────────────────────────────────────────────────────────
    basis_pwc = PWCBasis(mesh)
    N_pwc     = num_basis(basis_pwc)  # 3 * num_tets
    I_pwc     = ones(CT, N_pwc)
    perms_pwc = fill(CT(ε_r), num_elements(mesh))

    J_pwc = geoVolumeCurrentCal(I_pwc, basis_pwc, perms_pwc)

    @test length(J_pwc) == num_elements(mesh)
    @test eltype(J_pwc) <: SVector{3}
    @test all(all(isfinite, j) for j in J_pwc)
    # 所有 I=1 → J = κ₀*(ε_r-1)*[1,1,1] 标准化后的量级
    # 实际上 PWC 中 I 直接对应 J*V，所以 J_tet = [I_x/V, I_y/V, I_z/V]
    # 这里每个分量 I=1，所以 J_x = J_y = J_z = 1/V
    # 验证各分量相等
    Jx = [real(j[1]) for j in J_pwc]
    Jy = [real(j[2]) for j in J_pwc]
    Jz = [real(j[3]) for j in J_pwc]
    @test all(Jx .≈ Jy)
    @test all(Jy .≈ Jz)
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 17.5 — AntennaMetrics
# ─────────────────────────────────────────────────────────────────────────────
@testset "AntennaMetrics" begin
    using EMSuite.PostProcessing: antenna_directivity, input_impedance, beam_metrics

    freq = 300e6
    set_frequency!(freq)

    # ─── beam_metrics ────────────────────────────────────────────────────────
    @testset "beam_metrics" begin
        # Synthetic 1-D pattern: main lobe at index 5, two side lobes
        θs  = collect(range(0.0, π, length=21))
        # Peak (0 dB) at center; parabolic main lobe, side lobes ≈ -15 dB
        pat = -30.0 .* ((θs .- π/2) / (π/4)).^2
        pat = clamp.(pat, -40.0, 0.0)
        # Add artificial side lobe
        pat[2] = -15.0

        bm = beam_metrics(θs, pat)

        @test bm isa NamedTuple
        @test haskey(bm, :peak_angle)
        @test haskey(bm, :HPBW)
        @test haskey(bm, :SLL_dB)

        # Peak should be at π/2
        @test bm.peak_angle ≈ π/2 atol=0.01

        # HPBW > 0
        @test bm.HPBW > 0

        # SLL_dB should be negative (side lobe below main lobe)
        @test bm.SLL_dB < 0

        # With only one lobe (monotone decay), SLL_dB = 0
        mono_pat = -collect(range(0.0, 30.0, length=10))
        bm2 = beam_metrics(collect(range(0.0, π, length=10)), mono_pat)
        @test bm2.SLL_dB == 0.0
    end

    # ─── input_impedance ─────────────────────────────────────────────────────
    @testset "input_impedance" begin
        # Minimal RWG mesh (2 triangles → 1 basis function)
        nodes    = [0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
        elements = [1 2; 2 4; 3 3]
        mesh     = TriangleMesh(2, nodes, elements, [1,1])
        basis    = RWGBasis(mesh)

        # DeltaGapSource on edge 1, voltage = 1.0 V
        source  = DeltaGapSource(freq, [1], 1.0 + 0.0im)
        ICoeff  = ComplexF64[0.1 - 0.01im]   # mock MoM solution

        Z_in = input_impedance(source, ICoeff, basis)

        # Z_in = 1.0 / (0.1 - 0.01im)
        expected = 1.0 / (0.1 - 0.01im)
        @test Z_in ≈ expected

        # Z_in is Complex
        @test Z_in isa Complex
    end

    # ─── antenna_directivity ─────────────────────────────────────────────────
    @testset "antenna_directivity" begin
        # Minimal RWG mesh
        nodes    = [0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
        elements = [1 2; 2 4; 3 3]
        mesh     = TriangleMesh(2, nodes, elements, [1,1])
        basis    = RWGBasis(mesh)
        ICoeff   = ComplexF64[1.0 + 1.0im]

        # Coarse θ/ϕ grid (needs ≥ 2 per axis)
        θs = collect(range(0.0, π, length=5))
        ϕs = collect(range(0.0, 2π, length=5))

        result = antenna_directivity(θs, ϕs, ICoeff, basis)

        @test result isa NamedTuple
        @test haskey(result, :D)
        @test haskey(result, :P_rad)
        @test !haskey(result, :G)        # P_input not supplied

        # Shape
        @test size(result.D) == (length(θs), length(ϕs))

        # D ≥ 0
        @test all(result.D .>= 0)
        @test all(isfinite, result.D)
        @test result.P_rad >= 0

        # ── With P_input → also get G and η_eff ──────────────────────────────
        result2 = antenna_directivity(θs, ϕs, ICoeff, basis; P_input=result.P_rad)

        @test haskey(result2, :G)
        @test haskey(result2, :η_eff)
        # η_eff should be 1.0 when P_input = P_rad
        @test result2.η_eff ≈ 1.0 atol=1e-10
        # G ≈ D when η_eff = 1
        @test result2.G ≈ result.D atol=1e-10
    end
end