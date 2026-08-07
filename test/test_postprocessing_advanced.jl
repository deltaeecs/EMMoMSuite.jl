"""
    test_postprocessing_advanced.jl

Phase 21 测试：高级后处理与快速算法
- SolverResult            (21.0)
- NearFieldGrid/Line      (21.1)
- FarFieldPattern + 极化  (21.2)
- RCSBatch               (21.3)
- MLFMACache             (21.4)
"""

using Test
using EMMoMSuite
using StaticArrays
using LinearAlgebra

@testset "Phase 21 — 高级后处理与快速算法" begin

    # ─────────────────────────────────────────────────────────────────────────
    # 21.0  SolverResult
    # ─────────────────────────────────────────────────────────────────────────
    @testset "SolverResult" begin
        # 用假的 basis-like 对象和系数构造 SolverResult
        dummy_coeffs  = ComplexF64[1.0+0im, 0.5-0.3im, -0.2+0.1im]
        dummy_basis   = nothing   # phase 21 allows nothing for unit tests
        sr = SolverResult(dummy_basis, dummy_coeffs, 1.0e9)

        @test sr.freq ≈ 1.0e9
        @test sr.I_coeffs === dummy_coeffs
        @test length(sr) == 3

        # 带 metadata 关键字参数构造
        sr2 = SolverResult(dummy_basis, dummy_coeffs, 2.4e9;
                           metadata = Dict("equation" => "EFIE"))
        @test sr2.metadata["equation"] == "EFIE"
        @test sr2.freq ≈ 2.4e9

        # 默认 metadata 为空 Dict
        @test isempty(sr.metadata)

        # 负频率应抛出 ArgumentError
        @test_throws ArgumentError SolverResult(dummy_basis, dummy_coeffs, -1.0)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 21.1  NearFieldGrid / NearFieldLine
    # ─────────────────────────────────────────────────────────────────────────
    @testset "NearFieldGrid" begin
        Nu, Nv = 5, 7
        coords_u = range(0.0, 1.0; length=Nu) |> collect
        coords_v = range(0.0, 2.0; length=Nv) |> collect
        E = [SVector{3,ComplexF64}(ComplexF64(i+j*im), 0, 0) for i in 1:Nu, j in 1:Nv]
        H = zeros(SVector{3, ComplexF64}, Nu, Nv)
        grid = NearFieldGrid(coords_u, coords_v, E, H, 1e9, Dict{String,Any}())

        @test size(grid.E_field) == (Nu, Nv)
        @test size(grid.H_field) == (Nu, Nv)
        @test grid.freq ≈ 1e9
        @test length(grid.coords_u) == Nu
        @test length(grid.coords_v) == Nv

        # 维度不匹配应抛出 ArgumentError
        bad_coords = [0.0, 1.0]   # length 2 ≠ Nu=5
        @test_throws ArgumentError NearFieldGrid(bad_coords, coords_v, E, H, 1e9, Dict{String,Any}())
    end

    @testset "NearFieldLine" begin
        N = 10
        arc = collect(range(0.0, 1.0; length=N))
        pts = [SA[Float64(t), 0.0, 0.0] for t in arc]
        E   = [SVector{3,ComplexF64}(ComplexF64(t), 0, 0) for t in arc]
        H   = zeros(SVector{3, ComplexF64}, N)
        line = NearFieldLine(arc, pts, E, H, 500e6)

        @test length(line.arc_length) == N
        @test length(line.E_field) == N
        @test line.freq ≈ 500e6

        # 长度不一致应抛出 ArgumentError
        @test_throws ArgumentError NearFieldLine(arc, pts[1:3], E, H, 1e9)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 21.2  FarFieldPattern + 极化
    # ─────────────────────────────────────────────────────────────────────────
    @testset "FarFieldPattern 构造" begin
        freqs   = [300e6, 600e6]
        theta   = collect(range(0, π; length=37))
        phi     = collect(range(0, 2π; length=73))
        Nf, Nθ, Nφ = 2, 37, 73
        E_theta = randn(ComplexF64, Nf, Nθ, Nφ)
        E_phi   = randn(ComplexF64, Nf, Nθ, Nφ)

        ff = FarFieldPattern(freqs, theta, phi, E_theta, E_phi)

        @test size(ff.E_theta) == (Nf, Nθ, Nφ)
        @test length(ff.freqs) == Nf
        @test length(ff.theta) == Nθ
        @test length(ff.phi)   == Nφ
        # 初始缓存为 nothing
        @test isnothing(ff._gain)
        @test isnothing(ff._axial_ratio)

        # 尺寸不匹配应抛出
        @test_throws ArgumentError FarFieldPattern(freqs, theta, phi,
            randn(ComplexF64, 3, Nθ, Nφ), E_phi)  # Nf mismatch
    end

    @testset "FarFieldPattern gain" begin
        # 均匀各向同性辐射体：所有方向 |E_θ|² + |E_φ|² = 1 → D ≈ 1 (0 dBi)
        freqs = [300e6]
        theta = collect(range(1e-3, π-1e-3; length=73))  # 避开 θ=0/π 奇点
        phi   = collect(range(0.0, 2π*(1-1/73); length=73))
        Nθ, Nφ = 73, 73
        # 均匀场：E_θ = 1/√2, E_φ = 0 → U∝ 1/(2η₀) 均匀
        E_theta = fill(ComplexF64(1/sqrt(2)), 1, Nθ, Nφ)
        E_phi   = zeros(ComplexF64, 1, Nθ, Nφ)

        ff = FarFieldPattern(freqs, theta, phi, E_theta, E_phi)
        g_matrix = gain(ff; freq_idx=1)   # Nθ × Nφ

        @test size(g_matrix) == (Nθ, Nφ)
        # 均匀辐射体增益 ≈ 0 dBi（允许数值积分误差 0.3 dB）
        g_max = maximum(g_matrix)
        g_min = minimum(g_matrix)
        @test abs(g_max - 0.0) < 0.3    # dBi
        @test abs(g_min - 0.0) < 0.3

        # 便捷函数
        g_pt = gain_db(ff, theta[37], phi[1]; freq_idx=1)
        @test abs(g_pt - 0.0) < 0.3
    end

    @testset "FarFieldPattern hpbw + sll" begin
        # 合成一个简单方向图：在 θ=π/2 附近有主瓣，理论 HPBW 可用解析估计
        freqs = [1e9]
        Nθ = 181
        Nφ = 73  # ≥ 2，避免单点 φ 积分退化
        theta = collect(range(0.0, π; length=Nθ))
        phi   = collect(range(0.0, 2π * (1 - 1/Nφ); length=Nφ))
        # 余弦包络（各向同性对 φ）：在 θ=π/2 时最大
        # 功率 ∝ |E_θ|² = cos²(θ-π/2)，半功率点在 ±45° → HPBW ≈ 90°
        E_theta = [cos(θ - π/2) + 0im for θ in theta, φ in phi]  # Nθ×Nφ
        E_theta = reshape(E_theta, 1, Nθ, Nφ)
        E_phi   = zeros(ComplexF64, 1, Nθ, Nφ)

        ff = FarFieldPattern(freqs, theta, phi, E_theta, E_phi)
        hp = hpbw(ff; plane=:E, freq_idx=1)   # E-plane cut（沿 θ 方向，phi=0）

        # cos²(θ-π/2)=0.5 → θ-π/2 = ±45° → HPBW ≈ 90°
        @test abs(hp - 90.0) < 5.0  # 允许 5° 误差（数值积分离散化）

        sll_val = side_lobe_level(ff; plane=:E, freq_idx=1)
        # cos² 无旁瓣，返回 -Inf 或很低
        @test sll_val < -20.0  # dB
    end

    @testset "axial_ratio 圆极化" begin
        # 纯左旋圆极化：E_θ = 1，E_φ = jm → AR = 1 (0 dB)
        freqs = [1e9]; Nθ=5; Nφ=5
        theta = collect(range(0.1, π-0.1; length=Nθ))
        phi   = collect(range(0.0, 2π*(1-1/Nφ); length=Nφ))
        E_theta = ones(ComplexF64, 1, Nθ, Nφ)
        E_phi   = fill(ComplexF64(0+1im), 1, Nθ, Nφ)   # 与 E_θ 幅度相等、相位差 90°

        ff = FarFieldPattern(freqs, theta, phi, E_theta, E_phi)
        ar = axial_ratio(ff; freq_idx=1)   # Nθ × Nφ  dB

        @test size(ar) == (Nθ, Nφ)
        # LHCP: AR = 0 dB（允许 0.05 dB 误差）
        @test all(abs.(ar) .< 0.05)
    end

    @testset "co_cross_decompose (Ludwig III)" begin
        # 纯 θ-极化：Ludwig III 垂直 = E_θ，水平 = E_φ（phi=0 切面）
        freqs = [1e9]; Nθ=5; Nφ=1
        theta = collect(range(0.1, π-0.1; length=Nθ))
        phi   = [0.0]
        E_theta = ones(ComplexF64, 1, Nθ, 1)
        E_phi   = zeros(ComplexF64, 1, Nθ, 1)

        ff = FarFieldPattern(freqs, theta, phi, E_theta, E_phi)
        decomp = co_cross_decompose(ff; freq_idx=1)

        @test haskey(decomp, :co)
        @test haskey(decomp, :cross)
        # 纯 co 极化：cross 功率应为零
        cross_power = sum(abs2, decomp.cross)
        co_power    = sum(abs2, decomp.co)
        @test cross_power < 1e-20
        @test co_power > 0

        # XPD（dB）趋向 Inf 或一个很大的值
        xpd_val = xpd(ff; freq_idx=1)
        @test xpd_val > 40.0  # dB
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 21.3  RCSBatch
    # ─────────────────────────────────────────────────────────────────────────
    @testset "RCSResult 构造与查询" begin
        freqs     = [300e6, 600e6, 900e6]
        theta_inc = [0.0, π/4, π/2]
        phi_inc   = [0.0]
        Nf, Nθ, Nφ = 3, 3, 1

        rcs_vv = ones(Float64, Nf, Nθ, Nφ)          # 1 m² = 0 dBsm
        rcs_hh = 2ones(Float64, Nf, Nθ, Nφ)
        rcs_vh = zeros(Float64, Nf, Nθ, Nφ)
        rcs_hv = zeros(Float64, Nf, Nθ, Nφ)
        result = RCSResult(freqs, theta_inc, phi_inc, rcs_vv, rcs_hh, rcs_vh, rcs_hv)

        @test length(result.freqs) == Nf
        @test size(result.rcs_vv)  == (Nf, Nθ, Nφ)
        @test size(result.rcs_hh)  == (Nf, Nθ, Nφ)

        # frequency response（theta_idx=1, phi_idx=1, VV）
        f_ax, rcs_ax = rcs_frequency_response(result; theta_idx=1, phi_idx=1, polarization=:VV)
        @test length(f_ax)   == Nf
        @test length(rcs_ax) == Nf
        @test all(rcs_ax .≈ 0.0)  # 10*log10(1) = 0 dBsm

        f_ax_hh, rcs_ax_hh = rcs_frequency_response(result; theta_idx=2, phi_idx=1, polarization=:HH)
        @test all(rcs_ax_hh .≈ 10*log10(2))   # ≈ 3 dBsm

        # angular pattern（freq_idx=1, VV, E-plane）
        ang_ax, rcs_ang = rcs_angular_pattern(result; freq_idx=1, polarization=:VV, cut_plane=:E)
        @test length(ang_ax)   == Nθ
        @test length(rcs_ang)  == Nθ

        # 维度不匹配应抛出
        bad_rcs = ones(Float64, 2, Nθ, Nφ)  # Nf=2 ≠ 3
        @test_throws ArgumentError RCSResult(freqs, theta_inc, phi_inc,
                                             bad_rcs, rcs_hh, rcs_vh, rcs_hv)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 21.4  MLFMACache
    # ─────────────────────────────────────────────────────────────────────────
    @testset "MLFMACache 结构体与有效性" begin
        # 构造一个假的 cache（不需要真正 MLFMA 树）
        cache = MLFMACache(;
            freq      = 300e6,
            is_valid  = true,
            metadata  = Dict("n_basis" => 100),
        )

        @test cache.freq ≈ 300e6
        @test cache.is_valid == true
        @test cache.metadata["n_basis"] == 100

        # invalidate_cache! 应将 is_valid 设为 false
        invalidate_cache!(cache)
        @test cache.is_valid == false

        # 频率变化时，validate_cache 应返回 false
        @test !validate_cache(cache, 300e6)   # 已被 invalidate

        cache2 = MLFMACache(; freq=300e6, is_valid=true)
        @test validate_cache(cache2, 300e6)   # same freq → valid
        @test !validate_cache(cache2, 600e6)  # 频率不同 → invalid
    end

    @testset "solve_multi_rhs! 等价性" begin
        # 构造一个小 2×2 系统，多 RHS 求解结果与逐次求解应一致
        A  = ComplexF64[2+im  0.5; 0.5  3-0.5im]
        B  = ComplexF64[1 0; 0 1; 1 -1]'  # 2×3 rhs 矩阵

        # 逐列参考解
        X_ref = hcat([A \ B[:, k] for k in 1:3]...)

        # via solve_multi_rhs! (stub: for non-MLFMA case, falls back to factorization)
        X_out = solve_multi_rhs(A, B)

        @test X_out ≈ X_ref  atol=1e-12
    end

end  # @testset Phase 21
