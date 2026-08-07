"""
test_pmchw_excitation.jl — Phase 15 步骤 15.1 / 15.3

测试覆盖：
  15.1/15.2  excitation_vector(op::PMCHW, source::DeltaGapSource, basis::RWGBasis)
             → 长度 2N 向量；前 N 行（E 方程）有 delta-gap；后 N 行（H 方程）全零
  15.3       input_impedance(op::PMCHW, source, I_2N, basis)
             → 仅使用 J 部分（前 N 个系数），M 部分被忽略
"""

using Test
using EMMoMSuite
using LinearAlgebra

@testset "15.1/15.2/15.3 PMCHW DeltaGap 激励 + input_impedance" begin

    # ── 测试用小网格（4×6 球面三角剖分，快速） ─────────────────────────
    mesh  = generate_sphere_mesh(0.1, 4, 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)
    @test N > 0

    freq  = 300e6
    pmchw = PMCHW(freq, 4.0)          # εᵣ = 4（无损介质）

    V_feed  = 1.0 + 0im
    feed    = DeltaGapSource(freq, [1], V_feed)   # 第 1 条 RWG 边作馈电边
    ℓ₁      = basis.functions[1].edge_length       # 第 1 条边的长度

    # ─────────────────────────────────────────────────────────────────────
    # 15.2: excitation_vector(PMCHW, DeltaGapSource, RWGBasis)
    # ─────────────────────────────────────────────────────────────────────
    @testset "vector length is 2N" begin
        V = excitation_vector(pmchw, feed, basis)
        @test length(V) == 2N
    end

    @testset "V_E feed edge has correct value" begin
        V = excitation_vector(pmchw, feed, basis)
        @test V[1] ≈ V_feed * ℓ₁
    end

    @testset "V_E non-feed entries are zero" begin
        V = excitation_vector(pmchw, feed, basis)
        if N > 1
            @test all(iszero, V[2:N])
        end
    end

    @testset "V_H part (back N entries) is all zeros" begin
        V = excitation_vector(pmchw, feed, basis)
        @test all(iszero, V[N+1:2N])
    end

    @testset "multiple feed edges all land in V_E part" begin
        feed2 = DeltaGapSource(freq, [1, 2], V_feed)
        V     = excitation_vector(pmchw, feed2, basis)
        ℓ₂    = basis.functions[2].edge_length
        @test V[1] ≈ V_feed * ℓ₁
        @test V[2] ≈ V_feed * ℓ₂
        @test all(iszero, V[N+1:2N])
    end

    @testset "out-of-range edge index is silently skipped" begin
        bad = DeltaGapSource(freq, [N + 1], V_feed)   # index N+1 超出 RWGBasis 范围
        V   = excitation_vector(pmchw, bad, basis)
        @test all(iszero, V)
    end

    # ─────────────────────────────────────────────────────────────────────
    # 15.3: input_impedance(PMCHW, DeltaGapSource, I_2N, RWGBasis)
    # ─────────────────────────────────────────────────────────────────────
    @testset "input_impedance uses J part only" begin
        # 构造合成解向量：
        #   I_J[1] = 2.0 / ℓ₁  → 加权电流 I_J[1] * ℓ₁ = 2.0
        #   I_M[1] = 999 + 999j （M 部分应被忽略）
        I_2N = zeros(ComplexF64, 2N)
        I_2N[1]   = 2.0 / ℓ₁
        I_2N[N+1] = 999.0 + 999.0im   # M 部分（应被忽略）
        Z_in = input_impedance(pmchw, feed, I_2N, basis)
        @test Z_in ≈ ComplexF64(V_feed) / 2.0   # V / (I_J[1] * ℓ₁) = 1.0 / 2.0
    end

    @testset "input_impedance correct with complex voltage" begin
        V_c    = 2.0 + 1.0im
        feed_c = DeltaGapSource(freq, [1], V_c)
        I_2N   = zeros(ComplexF64, 2N)
        I_2N[1] = 4.0 / ℓ₁   # 加权电流 = 4
        Z_in = input_impedance(pmchw, feed_c, I_2N, basis)
        @test Z_in ≈ V_c / 4.0
    end

    @testset "input_impedance throws when J feed current is zero" begin
        I_zero = zeros(ComplexF64, 2N)
        I_zero[N+1] = 1.0  # 只有 M 部分，J 部分为零
        @test_throws ErrorException input_impedance(pmchw, feed, I_zero, basis)
    end

    @testset "input_impedance Re(Z_in) positive for resistive load" begin
        # 简单正阻: I_J[1] = 1/Z_true, Re(Z_true)>0
        Z_true = 50.0 + 25.0im
        I_2N   = zeros(ComplexF64, 2N)
        I_2N[1] = V_feed / (Z_true * ℓ₁)  # I_J * ℓ₁ = V / Z_true
        Z_in = input_impedance(pmchw, feed, I_2N, basis)
        @test real(Z_in) > 0
        @test Z_in ≈ Z_true  rtol=1e-10
    end
end
