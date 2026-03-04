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

        # Z_in = V / (I_coeff * edge_length)
        # Shared edge: node2=(1,0,0) to node3=(0,1,0) → length = √2
        edge_len = sqrt(2.0)
        expected = 1.0 / ((0.1 - 0.01im) * edge_len)
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

# ─────────────────────────────────────────────────────────────────────────────
# Phase 17.8 — absorbed_power / SAR
# ─────────────────────────────────────────────────────────────────────────────
@testset "Absorption" begin
    using EMSuite.PostProcessing: absorbed_power, sar

    freq = 300e6
    set_frequency!(freq)
    ω   = 2π * freq
    c0  = 299792458.0
    μ₀  = 4π * 1e-7
    ε₀  = 1.0 / (c0^2 * μ₀)
    CT  = ComplexF64

    # 最小 Tet 网格：2×2×2 节点 → 6 tets，每个体积 = 1/6
    mesh  = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
    N_tet = num_elements(mesh)

    # 辅助：提取每个 tet 的实际体积
    function tet_vols(m)
        vs = [abs(dot(m.node[:, m.tetras[2,t]] - m.node[:, m.tetras[1,t]],
                      cross(m.node[:, m.tetras[3,t]] - m.node[:, m.tetras[1,t]],
                            m.node[:, m.tetras[4,t]] - m.node[:, m.tetras[1,t]]))) / 6
              for t in 1:num_elements(m)]
        return vs
    end
    V = tet_vols(mesh)

    ε_r   = 2.0 - 1.0im   # lossy: ε_r'' = 1.0, |ε_r-1|² = |(1-j)|² = 2
    perms = fill(CT(ε_r), N_tet)

    # ─── PWC: 精确解析验证 ─────────────────────────────────────────────────
    basis_pwc = PWCBasis(mesh)
    N_pwc     = num_basis(basis_pwc)   # = 3 * N_tet
    I_coeffs  = ones(CT, N_pwc)

    # J_pol[t] = κ * [1,1,1] where κ = (ε_r-1)/ε_r
    # E_t = J_pol / (jωε₀(ε_r-1)) → |E_t|² = |J_pol|²/(ω²ε₀²|ε_r-1|²) = 3|κ|²/(ω²ε₀²|ε_r-1|²)
    # P_abs,t = (ω/2) ε₀ ε_r'' |E_t|² V_t
    #         = ε_r'' |J_pol|² V_t / (2ωε₀|ε_r-1|²)
    # = (-Im(ε_r)) * 3 * |κ|² V_t / (2ωε₀|ε_r-1|²)
    # = ε_r'' * 3 / |ε_r|² V_t / (2ωε₀)
    # with ε_r = 2-j: -Im=1, |ε_r|²=5 → P_abs,t = 3V/(10ωε₀)
    eps_pp      = -imag(ε_r)             # 1.0
    expected_Pd = eps_pp * 3.0 / (abs2(ε_r) * 2 * ω * ε₀)   # W/m³

    result = absorbed_power(basis_pwc, I_coeffs, perms)

    @test result isa NamedTuple
    @test haskey(result, :P_total)
    @test haskey(result, :P_density)
    @test length(result.P_density) == N_tet
    @test all(result.P_density .>= 0)

    # 各单元 P_density 应等于解析值
    @test result.P_density ≈ fill(expected_Pd, N_tet) rtol=1e-10

    # 能量守恒：P_total = sum(P_density * V)
    @test result.P_total ≈ sum(result.P_density .* V) rtol=1e-10

    # ─── 无损材料：P_total = 0 ───────────────────────────────────────────────
    perms_ll  = fill(CT(2.0 + 0.0im), N_tet)
    result_ll = absorbed_power(basis_pwc, I_coeffs, perms_ll)
    @test result_ll.P_total < 1e-30
    @test all(result_ll.P_density .< 1e-30)

    # ─── SAR ─────────────────────────────────────────────────────────────────
    ρ_body    = 1000.0           # kg/m³
    rhos      = fill(ρ_body, N_tet)
    sar_result = sar(basis_pwc, I_coeffs, perms, rhos)

    @test sar_result isa NamedTuple
    @test haskey(sar_result, :SAR_total)
    @test haskey(sar_result, :SAR_per_element)
    @test length(sar_result.SAR_per_element) == N_tet

    # SAR_per_element[t] = P_density[t] / ρ
    @test sar_result.SAR_per_element ≈ result.P_density ./ ρ_body rtol=1e-10

    # SAR_total = P_total / total_mass
    total_mass = sum(rhos .* V)
    @test sar_result.SAR_total ≈ result.P_total / total_mass rtol=1e-10

    # ─── SWG: 结构性测试 ─────────────────────────────────────────────────────
    basis_swg = SWGBasis(mesh)
    N_swg     = num_basis(basis_swg)
    I_swg     = ones(CT, N_swg)
    result_swg = absorbed_power(basis_swg, I_swg, perms)

    @test result_swg.P_total > 0
    @test length(result_swg.P_density) == N_tet
    @test all(isfinite, result_swg.P_density)
    @test all(result_swg.P_density .>= 0)
    @test result_swg.P_total ≈ sum(result_swg.P_density .* V) rtol=1e-10
end