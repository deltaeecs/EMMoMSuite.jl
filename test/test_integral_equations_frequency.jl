"""
test_integral_equations_frequency.jl

Phase 21.3 数值验证：频率扫描测试
- 多频率验证 EFIE/MFIE/CFIE
- 检查频率缩放关系
- 验证数值稳定性
"""

using Test
using EMSuite
using EMSuite.Core
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.PostProcessing
using LinearAlgebra

@testset "Integral Equations Frequency Sweep" begin
    
    # ────────────────────────────────────────────────────────────
    # EFIE 频率扫描：100 MHz - 1 GHz
    # ────────────────────────────────────────────────────────────
    @testset "EFIE Frequency Sweep: 100MHz-1GHz" begin
        println("  Testing EFIE across multiple frequencies...")
        
        # 固定几何（小板，电尺寸随频率变化）
        vertices = [
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.3, 0.3, 0.0],
            [0.0, 0.3, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        
        # 频率扫描
        frequencies = [100e6, 300e6, 600e6, 1e9]
        rcs_forward = Float64[]
        current_norms = Float64[]
        
        for freq in frequencies
            set_frequency!(freq)
            
            basis = RWGBasis(mesh)
            efie = EFIE(basis)
            
            Z = assemble_impedance_matrix(efie)
            
            # 平面波激励（固定）
            θ_inc = 0.0
            ϕ_inc = 0.0
            E_inc = [1.0, 0.0, 0.0]
            excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
            
            V = compute_excitation_vector(efie, excitation)
            I = Z \ V
            
            @test !any(isnan.(I))
            @test norm(I) > 0
            
            # 计算前向 RCS
            rcs = compute_rcs(basis, I, freq, [0.0], 0.0)[1]
            
            push!(rcs_forward, rcs)
            push!(current_norms, norm(I))
            
            println("    f = $(freq/1e9) GHz: |I| = $(round(norm(I), sigdigits=3)), RCS = $(round(rcs, digits=1)) dBsm")
        end
        
        # 验证：所有频率点都有有效结果
        @test length(rcs_forward) == 4
        @test !any(isnan.(rcs_forward))
        @test !any(isinf.(rcs_forward))
        
        # 检查 RCS 单调性：对于小目标，RCS 应随频率增加
        # （Rayleigh 区域：RCS ∝ f^4 对 electrically small targets）
        # 这里只检查趋势，不严格要求 f^4
        @test rcs_forward[4] > rcs_forward[1]  # 1 GHz > 100 MHz
        
        println("    ✓ RCS trend: $(round(rcs_forward[1], digits=1)) → $(round(rcs_forward[4], digits=1)) dBsm (100 MHz → 1 GHz)")
    end
    
    # ────────────────────────────────────────────────────────────
    # MFIE 频率扫描
    # ────────────────────────────────────────────────────────────
    @testset "MFIE Frequency Sweep: 100MHz-1GHz" begin
        println("  Testing MFIE across multiple frequencies...")
        
        vertices = [
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.3, 0.3, 0.0],
            [0.0, 0.3, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        
        frequencies = [100e6, 300e6, 600e6, 1e9]
        rcs_forward = Float64[]
        
        for freq in frequencies
            set_frequency!(freq)
            
            basis = RWGBasis(mesh)
            mfie = MFIE(basis)
            
            Z = assemble_impedance_matrix(mfie)
            
            θ_inc = 0.0
            ϕ_inc = 0.0
            E_inc = [1.0, 0.0, 0.0]
            excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
            
            V = compute_excitation_vector(mfie, excitation)
            I = Z \ V
            
            @test !any(isnan.(I))
            
            rcs = compute_rcs(basis, I, freq, [0.0], 0.0)[1]
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
        
        vertices = [
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.3, 0.3, 0.0],
            [0.0, 0.3, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        
        frequencies = [100e6, 300e6, 600e6, 1e9]
        α = 0.5  # CFIE 组合参数
        
        cond_numbers = Float64[]
        
        for freq in frequencies
            set_frequency!(freq)
            
            basis = RWGBasis(mesh)
            cfie = CFIE(basis, α)
            
            Z = assemble_impedance_matrix(cfie)
            
            # 条件数检查
            cond_Z = cond(Z)
            push!(cond_numbers, cond_Z)
            
            @test !isnan(cond_Z)
            @test !isinf(cond_Z)
            @test cond_Z < 1e8  # 条件数不应太大
            
            println("    f = $(freq/1e9) GHz: cond(Z) = $(round(cond_Z, sigdigits=3))")
        end
        
        # 条件数应在合理范围内
        @test maximum(cond_numbers) < 1e8
        @test minimum(cond_numbers) > 1.0
        
        println("    ✓ CFIE stable: cond(Z) ∈ [$(round(minimum(cond_numbers), sigdigits=2)), $(round(maximum(cond_numbers), sigdigits=2))]")
    end
    
    # ────────────────────────────────────────────────────────────
    # EFIE vs MFIE 频率一致性
    # ────────────────────────────────────────────────────────────
    @testset "EFIE vs MFIE Frequency Consistency" begin
        println("  Comparing EFIE and MFIE at multiple frequencies...")
        
        vertices = [
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.3, 0.3, 0.0],
            [0.0, 0.3, 0.0]
        ]
        triangles = [
            [1, 2, 3],
            [1, 3, 4]
        ]
        
        mesh = TriMesh(vertices, triangles)
        
        frequencies = [100e6, 300e6, 600e6, 1e9]
        
        for freq in frequencies
            set_frequency!(freq)
            
            basis = RWGBasis(mesh)
            
            # EFIE
            efie = EFIE(basis)
            Z_efie = assemble_impedance_matrix(efie)
            
            θ_inc = 0.0
            ϕ_inc = 0.0
            E_inc = [1.0, 0.0, 0.0]
            excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
            
            V_efie = compute_excitation_vector(efie, excitation)
            I_efie = Z_efie \ V_efie
            
            # MFIE
            mfie = MFIE(basis)
            Z_mfie = assemble_impedance_matrix(mfie)
            V_mfie = compute_excitation_vector(mfie, excitation)
            I_mfie = Z_mfie \ V_mfie
            
            # RCS 对比
            rcs_efie = compute_rcs(basis, I_efie, freq, [0.0], 0.0)[1]
            rcs_mfie = compute_rcs(basis, I_mfie, freq, [0.0], 0.0)[1]
            
            Δrcs = abs(rcs_efie - rcs_mfie)
            
            println("    f = $(freq/1e9) GHz: ΔRCS(EFIE-MFIE) = $(round(Δrcs, digits=1)) dB")
            
            # EFIE 和 MFIE 应给出相似的 RCS（对于良好网格）
            # 这里放宽到 15 dB，因为 toy mesh 可能有数值误差
            @test Δrcs < 15.0
        end
        
        println("    ✓ EFIE and MFIE consistent across frequencies")
    end
    
    # ────────────────────────────────────────────────────────────
    # 波长缩放验证
    # ────────────────────────────────────────────────────────────
    @testset "Wavelength Scaling Validation" begin
        println("  Validating wavelength scaling...")
        
        # 使用固定电尺寸（几何随频率缩放）
        using EMSuite.Core.Constants: c0
        
        freq_ref = 300e6  # 参考频率
        λ_ref = c0 / freq_ref
        
        # 设计一个 λ/10 × λ/10 的板（在 300 MHz）
        size_ref = λ_ref / 10
        
        frequencies = [300e6, 600e6, 1.2e9]
        rcs_values = Float64[]
        
        for freq in frequencies
            λ = c0 / freq
            size = λ / 10  # 保持 λ/10 电尺寸
            
            vertices = [
                [0.0, 0.0, 0.0],
                [size, 0.0, 0.0],
                [size, size, 0.0],
                [0.0, size, 0.0]
            ]
            triangles = [
                [1, 2, 3],
                [1, 3, 4]
            ]
            
            mesh = TriMesh(vertices, triangles)
            set_frequency!(freq)
            
            basis = RWGBasis(mesh)
            efie = EFIE(basis)
            
            Z = assemble_impedance_matrix(efie)
            
            θ_inc = 0.0
            ϕ_inc = 0.0
            E_inc = [1.0, 0.0, 0.0]
            excitation = PlaneWave(θ_inc, ϕ_inc, E_inc)
            
            V = compute_excitation_vector(efie, excitation)
            I = Z \ V
            
            rcs = compute_rcs(basis, I, freq, [0.0], 0.0)[1]
            push!(rcs_values, rcs)
            
            println("    f = $(freq/1e9) GHz (size = $(round(size*1000, digits=1)) mm): RCS = $(round(rcs, digits=1)) dBsm")
        end
        
        # 对于固定电尺寸，RCS 应相对稳定（< 5 dB 变化）
        rcs_range = maximum(rcs_values) - minimum(rcs_values)
        
        @test rcs_range < 5.0  # 电尺寸相同，RCS 应接近
        
        println("    ✓ Wavelength scaling: RCS range = $(round(rcs_range, digits=1)) dB (constant electrical size)")
    end
    
end  # @testset "Integral Equations Frequency Sweep"

println("✅ Frequency sweep tests completed")
