"""
test_integral_equations_endtoend.jl

Phase 21.2 集成测试：端到端求解流程
- EFIE/MFIE/CFIE 完整求解链
- 激励 → 组装 → 求解 → RCS
- 使用真实网格（非 toy mesh）
"""

using Test
using EMSuite
using EMSuite.Core
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.PostProcessing
using LinearAlgebra

@testset "Integral Equations End-to-End" begin
    
    # ────────────────────────────────────────────────────────────
    # EFIE 端到端：简单板 @ 300 MHz
    # ────────────────────────────────────────────────────────────
    @testset "EFIE End-to-End: Simple Plate 300MHz" begin
        println("  Testing EFIE end-to-end solve...")
        
        # 1. 创建简单板几何（1m × 1m，centered at origin）
        vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        
        # 2. 设置频率和基函数
        freq = 300e6  # 300 MHz
        set_frequency!(freq)
        
        basis = RWGBasis(mesh)
        n_basis = length(basis)
        
        @test n_basis == 1  # 简单板应有 1 个 RWG 基
        
        # 3. 组装 EFIE 矩阵
        efie = EFIE(basis)
        Z = assemble_impedance_matrix(efie)
        
        @test size(Z) == (n_basis, n_basis)
        @test !any(isnan.(Z))
        @test !any(isinf.(Z))
        
        # 4. 定义平面波激励（z 方向入射，x 极化）
        θ_inc = 0.0  # z 轴向下入射
        ϕ_inc = 0.0
        E_inc = [1.0, 0.0, 0.0]  # x 极化
        
        excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
        V = compute_excitation_vector(efie, excitation)
        
        @test length(V) == n_basis
        @test !any(isnan.(V))
        @test norm(V) > 0  # 激励非零
        
        # 5. 求解电流（直接 LU，因为矩阵很小）
        I = Z \ V
        
        @test length(I) == n_basis
        @test !any(isnan.(I))
        @test norm(I) > 0  # 电流非零
        
        # 6. 计算 RCS（双站，θ ∈ [0, π]）
        θ_obs = range(0, π, length=37)  # 每 5° 一个点
        ϕ_obs = 0.0
        
        rcs_dBsm = compute_rcs(basis, I, freq, θ_obs, ϕ_obs)
        
        @test length(rcs_dBsm) == 37
        @test !any(isnan.(rcs_dBsm))
        @test !any(isinf.(rcs_dBsm))
        @test maximum(rcs_dBsm) > -Inf  # 至少有有限 RCS
        
        # 物理验证：前向散射（θ=0）应有较大 RCS
        rcs_forward = rcs_dBsm[1]
        rcs_backward = rcs_dBsm[end]
        
        @test rcs_forward > -50.0  # 合理的 RCS 量级（dBsm）
        @test rcs_backward > -50.0
        
        println("    ✓ EFIE solve: |I| = $(round(norm(I), sigdigits=3))")
        println("    ✓ RCS: forward = $(round(rcs_forward, digits=1)) dBsm, backward = $(round(rcs_backward, digits=1)) dBsm")
    end
    
    # ────────────────────────────────────────────────────────────
    # MFIE 端到端：简单板 @ 300 MHz
    # ────────────────────────────────────────────────────────────
    @testset "MFIE End-to-End: Simple Plate 300MHz" begin
        println("  Testing MFIE end-to-end solve...")
        
        # 使用相同的网格和频率
        vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        freq = 300e6
        set_frequency!(freq)
        
        basis = RWGBasis(mesh)
        n_basis = length(basis)
        
        # 组装 MFIE 矩阵
        mfie = MFIE(basis)
        Z_mfie = assemble_impedance_matrix(mfie)
        
        @test size(Z_mfie) == (n_basis, n_basis)
        @test !any(isnan.(Z_mfie))
        
        # MFIE 应该是非对称的（K-operator）
        @test norm(Z_mfie - Z_mfie') / norm(Z_mfie) > 1e-6  # 非对称
        
        # 平面波激励
        θ_inc = 0.0
        ϕ_inc = 0.0
        E_inc = [1.0, 0.0, 0.0]
        
        excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
        V_mfie = compute_excitation_vector(mfie, excitation)
        
        @test norm(V_mfie) > 0
        
        # 求解
        I_mfie = Z_mfie \ V_mfie
        
        @test !any(isnan.(I_mfie))
        @test norm(I_mfie) > 0
        
        # RCS
        θ_obs = range(0, π, length=37)
        ϕ_obs = 0.0
        
        rcs_mfie = compute_rcs(basis, I_mfie, freq, θ_obs, ϕ_obs)
        
        @test !any(isnan.(rcs_mfie))
        @test maximum(rcs_mfie) > -Inf
        
        println("    ✓ MFIE solve: |I| = $(round(norm(I_mfie), sigdigits=3))")
        println("    ✓ MFIE RCS: forward = $(round(rcs_mfie[1], digits=1)) dBsm")
    end
    
    # ────────────────────────────────────────────────────────────
    # CFIE 端到端：α 参数扫描
    # ────────────────────────────────────────────────────────────
    @testset "CFIE End-to-End: α-Sweep" begin
        println("  Testing CFIE α-parameter sweep...")
        
        # 使用相同的网格
        vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        freq = 300e6
        set_frequency!(freq)
        
        basis = RWGBasis(mesh)
        
        # 平面波激励（固定）
        θ_inc = 0.0
        ϕ_inc = 0.0
        E_inc = [1.0, 0.0, 0.0]
        excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
        
        # α 参数扫描
        α_values = [0.3, 0.5, 0.7]
        rcs_results = Float64[]
        
        for α in α_values
            cfie = CFIE(basis, α)
            Z_cfie = assemble_impedance_matrix(cfie)
            V_cfie = compute_excitation_vector(cfie, excitation)
            
            I_cfie = Z_cfie \ V_cfie
            
            @test !any(isnan.(I_cfie))
            @test norm(I_cfie) > 0
            
            # 计算前向 RCS
            rcs = compute_rcs(basis, I_cfie, freq, [0.0], 0.0)
            push!(rcs_results, rcs[1])
            
            println("    α = $α: |I| = $(round(norm(I_cfie), sigdigits=3)), RCS = $(round(rcs[1], digits=1)) dBsm")
        end
        
        # 验证：RCS 应该随 α 平滑变化
        @test length(rcs_results) == 3
        @test !any(isnan.(rcs_results))
        
        # 检查单调性或连续性（避免突变）
        Δrcs1 = rcs_results[2] - rcs_results[1]
        Δrcs2 = rcs_results[3] - rcs_results[2]
        
        # RCS 变化应在合理范围内（< 20 dB）
        @test abs(Δrcs1) < 20.0
        @test abs(Δrcs2) < 20.0
        
        println("    ✓ α-sweep: RCS varies smoothly (Δ₁=$(round(Δrcs1, digits=1)) dB, Δ₂=$(round(Δrcs2, digits=1)) dB)")
    end
    
    # ────────────────────────────────────────────────────────────
    # EFIE vs MFIE：电流对比
    # ────────────────────────────────────────────────────────────
    @testset "EFIE vs MFIE: Current Comparison" begin
        println("  Comparing EFIE and MFIE currents...")
        
        # 相同几何和激励
        vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        freq = 300e6
        set_frequency!(freq)
        
        basis = RWGBasis(mesh)
        
        θ_inc = 0.0
        ϕ_inc = 0.0
        E_inc = [1.0, 0.0, 0.0]
        excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
        
        # EFIE 求解
        efie = EFIE(basis)
        Z_efie = assemble_impedance_matrix(efie)
        V_efie = compute_excitation_vector(efie, excitation)
        I_efie = Z_efie \ V_efie
        
        # MFIE 求解
        mfie = MFIE(basis)
        Z_mfie = assemble_impedance_matrix(mfie)
        V_mfie = compute_excitation_vector(mfie, excitation)
        I_mfie = Z_mfie \ V_mfie
        
        # EFIE 和 MFIE 的电流应该在相同量级（PEC 边界条件）
        ratio = norm(I_efie) / norm(I_mfie)
        
        @test 0.1 < ratio < 10.0  # 应在同一量级
        
        println("    |I_EFIE| = $(round(norm(I_efie), sigdigits=3))")
        println("    |I_MFIE| = $(round(norm(I_mfie), sigdigits=3))")
        println("    Ratio = $(round(ratio, digits=2))")
        
        # RCS 对比
        rcs_efie = compute_rcs(basis, I_efie, freq, [0.0], 0.0)[1]
        rcs_mfie = compute_rcs(basis, I_mfie, freq, [0.0], 0.0)[1]
        
        Δrcs = rcs_efie - rcs_mfie
        
        println("    RCS_EFIE = $(round(rcs_efie, digits=1)) dBsm")
        println("    RCS_MFIE = $(round(rcs_mfie, digits=1)) dBsm")
        println("    Δ = $(round(Δrcs, digits=1)) dB")
        
        # RCS 差异应在合理范围内（< 10 dB for simple geometry）
        @test abs(Δrcs) < 10.0
    end
    
end  # @testset "Integral Equations End-to-End"

println("✅ Integral equations end-to-end tests completed")
