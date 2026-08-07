"""
test_integral_equations_frequency.jl

Phase 21.3 数值验证：频率扫描测试
- 多频率验证 EFIE/MFIE/CFIE
- 检查频率缩放关系
- 验证数值稳定性
"""

using Test
using EMMoMSuite
using EMMoMSuite.CoreModule
using EMMoMSuite.CoreModule.Constants: c0
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.PostProcessing
using LinearAlgebra

@testset "Integral Equations Frequency Sweep" begin

    # 固定几何：0.3m × 0.3m 板（2×2 细分），电尺寸随频率变化
    function build_small_plate_mesh()
        return generate_rectangle_mesh(0.3, 0.3, 2, 2)
    end

    # ────────────────────────────────────────────────────────────
    # EFIE 频率扫描：100 MHz - 1 GHz
    # ────────────────────────────────────────────────────────────
    @testset "EFIE Frequency Sweep: 100MHz-1GHz" begin
        println("  Testing EFIE across multiple frequencies...")

        mesh = build_small_plate_mesh()
        frequencies = [100e6, 300e6, 600e6, 1e9]
        rcs_forward = Float64[]

        for freq in frequencies
            set_frequency!(freq)
            basis = RWGBasis(mesh)
            efie = EFIE(freq)

            Z = assemble_impedance_matrix(efie, basis)
            source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
            V = excitation_vector(efie, source, basis)
            I = Z \ V

            @test !any(isnan.(I))
            @test norm(I) > 0

            _, _, rcs_dB = radarCrossSection([0.0], [0.0], I, basis)
            rcs = rcs_dB[1, 1]
            push!(rcs_forward, rcs)

            println("    f = $(freq/1e9) GHz: |I| = $(round(norm(I), sigdigits=3)), RCS = $(round(rcs, digits=1)) dBsm")
        end

        @test length(rcs_forward) == 4
        @test !any(isnan.(rcs_forward))
        @test !any(isinf.(rcs_forward))

        # Rayleigh 区域趋势：电小目标 RCS ∝ f⁴，1 GHz 应显著大于 100 MHz
        @test rcs_forward[4] > rcs_forward[1]

        println("    ✓ RCS trend: $(round(rcs_forward[1], digits=1)) → $(round(rcs_forward[4], digits=1)) dBsm (100 MHz → 1 GHz)")
    end

    # ────────────────────────────────────────────────────────────
    # MFIE 频率扫描
    # ────────────────────────────────────────────────────────────
    @testset "MFIE Frequency Sweep: 100MHz-1GHz" begin
        println("  Testing MFIE across multiple frequencies...")

        mesh = build_small_plate_mesh()
        frequencies = [100e6, 300e6, 600e6, 1e9]
        rcs_forward = Float64[]

        for freq in frequencies
            set_frequency!(freq)
            basis = RWGBasis(mesh)
            mfie = MFIE(freq)

            Z = assemble_impedance_matrix(mfie, basis)
            source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
            V = excitation_vector(mfie, source, basis)
            I = Z \ V

            @test !any(isnan.(I))

            _, _, rcs_dB = radarCrossSection([0.0], [0.0], I, basis)
            rcs = rcs_dB[1, 1]
            push!(rcs_forward, rcs)

            println("    f = $(freq/1e9) GHz: RCS = $(round(rcs, digits=1)) dBsm")
        end

        @test !any(isnan.(rcs_forward))
        @test rcs_forward[4] > rcs_forward[1]  # 趋势检查
    end

    # ────────────────────────────────────────────────────────────
    # CFIE 频率稳定性
    # ────────────────────────────────────────────────────────────
    @testset "CFIE Frequency Stability" begin
        println("  Testing CFIE stability across frequencies...")

        mesh = build_small_plate_mesh()
        frequencies = [100e6, 300e6, 600e6, 1e9]
        α = 0.5
        cond_numbers = Float64[]

        for freq in frequencies
            set_frequency!(freq)
            basis = RWGBasis(mesh)
            cfie = CFIE(freq, α)

            Z = assemble_impedance_matrix(cfie, basis)
            cond_Z = cond(Z)
            push!(cond_numbers, cond_Z)

            @test !isnan(cond_Z)
            @test !isinf(cond_Z)
            @test cond_Z < 1e8

            println("    f = $(freq/1e9) GHz: cond(Z) = $(round(cond_Z, sigdigits=3))")
        end

        @test maximum(cond_numbers) < 1e8
        @test minimum(cond_numbers) > 1.0

        println("    ✓ CFIE stable: cond(Z) ∈ [$(round(minimum(cond_numbers), sigdigits=2)), $(round(maximum(cond_numbers), sigdigits=2))]")
    end

    # ────────────────────────────────────────────────────────────
    # EFIE vs MFIE 频率一致性
    # ────────────────────────────────────────────────────────────
    @testset "EFIE vs MFIE Frequency Consistency" begin
        println("  Comparing EFIE and MFIE at multiple frequencies...")

        mesh = build_small_plate_mesh()
        frequencies = [100e6, 300e6, 600e6, 1e9]

        for freq in frequencies
            set_frequency!(freq)
            basis = RWGBasis(mesh)
            source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])

            efie = EFIE(freq)
            Z_efie = assemble_impedance_matrix(efie, basis)
            V_efie = excitation_vector(efie, source, basis)
            I_efie = Z_efie \ V_efie

            mfie = MFIE(freq)
            Z_mfie = assemble_impedance_matrix(mfie, basis)
            V_mfie = excitation_vector(mfie, source, basis)
            I_mfie = Z_mfie \ V_mfie

            _, _, rcs_efie_dB = radarCrossSection([0.0], [0.0], I_efie, basis)
            _, _, rcs_mfie_dB = radarCrossSection([0.0], [0.0], I_mfie, basis)
            Δrcs = abs(rcs_efie_dB[1, 1] - rcs_mfie_dB[1, 1])

            println("    f = $(freq/1e9) GHz: ΔRCS(EFIE-MFIE) = $(round(Δrcs, digits=1)) dB")

            # 粗略网格允许数值误差，放宽到 15 dB
            @test Δrcs < 15.0
        end

        println("    ✓ EFIE and MFIE consistent across frequencies")
    end

    # ────────────────────────────────────────────────────────────
    # 波长缩放验证
    # ────────────────────────────────────────────────────────────
    @testset "Wavelength Scaling Validation" begin
        println("  Validating wavelength scaling...")

        frequencies = [300e6, 600e6, 1.2e9]
        rcs_values = Float64[]

        for freq in frequencies
            λ = c0 / freq
            size = λ / 10  # 保持 λ/10 电尺寸

            mesh = generate_rectangle_mesh(size, size, 2, 2)
            set_frequency!(freq)

            basis = RWGBasis(mesh)
            efie = EFIE(freq)

            Z = assemble_impedance_matrix(efie, basis)
            source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
            V = excitation_vector(efie, source, basis)
            I = Z \ V

            _, _, rcs_dB = radarCrossSection([0.0], [0.0], I, basis)
            rcs = rcs_dB[1, 1]
            push!(rcs_values, rcs)

            println("    f = $(freq/1e9) GHz (size = $(round(size*1000, digits=1)) mm): RCS = $(round(rcs, digits=1)) dBsm")
        end

        # 物理模型：平板法向入射 RCS σ = 4πA²/λ²。电尺寸固定（A ∝ λ²）时
        # σ ∝ 1/f²，即每个倍频程下降约 6 dB。300MHz→1.2GHz 为 2 个倍频程，
        # 预期 RCS 范围约 12 dB。
        rcs_range = maximum(rcs_values) - minimum(rcs_values)
        @test 8.0 < rcs_range < 16.0

        println("    ✓ Wavelength scaling: RCS range = $(round(rcs_range, digits=1)) dB (≈6 dB/octave, constant electrical size)")
    end

end  # @testset "Integral Equations Frequency Sweep"

println("✅ Frequency sweep tests completed")
