"""
test_integral_equations_endtoend.jl

Phase 21.2 集成测试：端到端求解流程
- EFIE/MFIE/CFIE 完整求解链
- 激励 → 组装 → 求解 → RCS
- 使用真实网格（generate_rectangle_mesh）
"""

using Test
using EMMoMSuite
using EMMoMSuite.CoreModule
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using LinearAlgebra

@testset "Integral Equations End-to-End" begin

    # 1m × 1m 板（2×2 细分，8 个三角形）
    function build_plate_mesh()
        return generate_rectangle_mesh(1.0, 1.0, 2, 2)
    end

    # ────────────────────────────────────────────────────────────
    # EFIE 端到端：简单板 @ 300 MHz
    # ────────────────────────────────────────────────────────────
    @testset "EFIE End-to-End: Simple Plate 300MHz" begin
        println("  Testing EFIE end-to-end solve...")

        mesh = build_plate_mesh()
        freq = 300e6
        set_frequency!(freq)

        basis = RWGBasis(mesh)
        n_basis = num_basis(basis)
        @test n_basis > 0

        # 组装 EFIE 矩阵
        efie = EFIE(freq)
        Z = assemble_impedance_matrix(efie, basis)
        @test size(Z) == (n_basis, n_basis)
        @test !any(isnan.(Z))
        @test !any(isinf.(Z))

        # 平面波激励（z 方向入射，x 极化）
        source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
        V = excitation_vector(efie, source, basis)
        @test length(V) == n_basis
        @test !any(isnan.(V))
        @test norm(V) > 0

        # 求解电流（直接 LU，因为矩阵很小）
        I = Z \ V
        @test length(I) == n_basis
        @test !any(isnan.(I))
        @test norm(I) > 0

        # 计算 RCS（双站，θ ∈ [0, π]，每 5° 一个点）
        θ_obs = collect(range(0.0, π, length = 37))
        ϕ_obs = [0.0]
        _, rcs_total, rcs_dB = radarCrossSection(θ_obs, ϕ_obs, I, basis)

        @test size(rcs_dB) == (37, 1)
        @test !any(isnan.(rcs_dB))
        @test !any(isinf.(rcs_dB))

        # 物理验证：前向散射（θ=0）应有合理量级（dBsm）
        rcs_forward = rcs_dB[1, 1]
        rcs_backward = rcs_dB[end, 1]
        @test rcs_forward > -50.0
        @test rcs_backward > -50.0

        println("    ✓ EFIE solve: |I| = $(round(norm(I), sigdigits=3))")
        println("    ✓ RCS: forward = $(round(rcs_forward, digits=1)) dBsm, backward = $(round(rcs_backward, digits=1)) dBsm")
    end

    # ────────────────────────────────────────────────────────────
    # MFIE 端到端：简单板 @ 300 MHz
    # ────────────────────────────────────────────────────────────
    @testset "MFIE End-to-End: Simple Plate 300MHz" begin
        println("  Testing MFIE end-to-end solve...")

        mesh = build_plate_mesh()
        freq = 300e6
        set_frequency!(freq)

        basis = RWGBasis(mesh)
        n_basis = num_basis(basis)

        # 组装 MFIE 矩阵
        mfie = MFIE(freq)
        Z_mfie = assemble_impedance_matrix(mfie, basis)
        @test size(Z_mfie) == (n_basis, n_basis)
        @test !any(isnan.(Z_mfie))

        # 注：对共面平板网格，K-operator 组装结果恰好为普通转置对称（实测 0.0）；
        # 非对称性在 3D 网格上验证（见 test_mfie_decomposition.jl）。
        @test norm(Z_mfie - transpose(Z_mfie)) / norm(Z_mfie) < 1e-6

        # 平面波激励
        source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
        V_mfie = excitation_vector(mfie, source, basis)
        @test norm(V_mfie) > 0

        I_mfie = Z_mfie \ V_mfie
        @test !any(isnan.(I_mfie))
        @test norm(I_mfie) > 0

        θ_obs = collect(range(0.0, π, length = 37))
        ϕ_obs = [0.0]
        _, _, rcs_mfie = radarCrossSection(θ_obs, ϕ_obs, I_mfie, basis)
        @test !any(isnan.(rcs_mfie))

        println("    ✓ MFIE solve: |I| = $(round(norm(I_mfie), sigdigits=3))")
        println("    ✓ MFIE RCS: forward = $(round(rcs_mfie[1, 1], digits=1)) dBsm")
    end

    # ────────────────────────────────────────────────────────────
    # CFIE 端到端：α 参数扫描
    # ────────────────────────────────────────────────────────────
    @testset "CFIE End-to-End: α-Sweep" begin
        println("  Testing CFIE α-parameter sweep...")

        # CFIE 在开平板上于 α=0.5 时激励向量数学上精确为零（开面退化），
        # 因此 α 扫描使用闭合球体网格，验证解对 α 的不变性。
        mesh = generate_sphere_mesh(0.5, 6, 12)
        freq = 300e6
        set_frequency!(freq)

        basis = RWGBasis(mesh)
        source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])

        # α 参数扫描
        α_values = [0.3, 0.5, 0.7]
        rcs_results = Float64[]

        for α in α_values
            cfie = CFIE(freq, α)
            Z_cfie = assemble_impedance_matrix(cfie, basis)
            V_cfie = excitation_vector(cfie, source, basis)
            I_cfie = Z_cfie \ V_cfie

            @test !any(isnan.(I_cfie))
            @test norm(I_cfie) > 0

            _, _, rcs_dB = radarCrossSection([0.0], [0.0], I_cfie, basis)
            rcs = rcs_dB[1, 1]
            push!(rcs_results, rcs)

            println("    α = $α: |I| = $(round(norm(I_cfie), sigdigits=3)), RCS = $(round(rcs, digits=1)) dBsm")
        end

        @test length(rcs_results) == 3
        @test !any(isnan.(rcs_results))

        # CFIE 的解应对 α 近似不变（同一物理问题的不同组合）
        rcs_range = maximum(rcs_results) - minimum(rcs_results)
        @test rcs_range < 5.0

        println("    ✓ α-sweep: RCS range = $(round(rcs_range, digits=2)) dB across α (α-invariance)")
    end

    # ────────────────────────────────────────────────────────────
    # EFIE vs MFIE：电流对比
    # ────────────────────────────────────────────────────────────
    @testset "EFIE vs MFIE: Current Comparison" begin
        println("  Comparing EFIE and MFIE currents...")

        mesh = build_plate_mesh()
        freq = 300e6
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

        # EFIE 和 MFIE 电流应在同一量级（PEC 边界条件）
        ratio = norm(I_efie) / norm(I_mfie)
        @test 0.1 < ratio < 10.0
        println("    |I_EFIE| = $(round(norm(I_efie), sigdigits=3))")
        println("    |I_MFIE| = $(round(norm(I_mfie), sigdigits=3))")
        println("    Ratio = $(round(ratio, digits=2))")

        _, _, rcs_efie_dB = radarCrossSection([0.0], [0.0], I_efie, basis)
        _, _, rcs_mfie_dB = radarCrossSection([0.0], [0.0], I_mfie, basis)
        rcs_efie = rcs_efie_dB[1, 1]
        rcs_mfie = rcs_mfie_dB[1, 1]
        Δrcs = rcs_efie - rcs_mfie

        println("    RCS_EFIE = $(round(rcs_efie, digits=1)) dBsm")
        println("    RCS_MFIE = $(round(rcs_mfie, digits=1)) dBsm")
        println("    Δ = $(round(Δrcs, digits=1)) dB")

        # RCS 差异应在合理范围内（< 10 dB for simple geometry）
        @test abs(Δrcs) < 10.0
    end

end  # @testset "Integral Equations End-to-End"

println("✅ Integral equations end-to-end tests completed")
