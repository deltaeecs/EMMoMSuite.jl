"""
test_coverage_gaps.jl

Round 6-12 修复后的回归测试套件，专注于覆盖缺口。
测试聚焦边界条件、参数验证、数值稳定性。

创建日期：2026-04-01
目的：防止 Phase 19 修复回归
"""

using Test
using EMSuite
using LinearAlgebra
using StaticArrays

@testset "Coverage Gaps - Round 6-12" begin

    # ───────────────────────────────────────────────────────────────────
    # P0.1: Singularities.jl 退化三角形
    # ───────────────────────────────────────────────────────────────────
    
    @testset "Singularities singularF1 - Degenerate Triangles" begin
        using EMSuite.IntegralEquations: singularF1
        
        # Case 1: Nearly degenerate (s ≈ a, makes 1-a/s → 0⁺)
        a, b, c = 1.0, 1.0, 1.999  # s=1.9995, 1-a/s ≈ 0.0005
        result = singularF1(a, b, c)
        @test isfinite(result) "F1 must be finite for near-degenerate triangle"
        @test !isnan(result)   "F1 must not be NaN"
        
        # Case 2: Equilateral (should be well-behaved)
        a = b = c = 1.0
        result = singularF1(a, b, c)
        @test isfinite(result)
        
        # Case 3: Highly degenerate (s barely > edge)
        a = 1.0; b = 1.0; c = 1.0 + 1e-10
        result = singularF1(a, b, c)
        @test isfinite(result)
        
        # Case 4: Long thin triangle (b ≈ c, a very small)
        a = 0.01; b = 1.0; c = 1.005
        result = singularF1(a, b, c)
        @test isfinite(result)
    end
    
    @testset "Singularities singularF21 - Degenerate Triangles" begin
        using EMSuite.IntegralEquations: singularF21
        
        # Degenerate: area2 → 0
        a, b, c = 1.0, 1.0, 1.999
        area2 = 1e-15  # Nearly zero area
        result = singularF21(a, b, c, area2)
        @test isfinite(result) "F21 must be finite despite near-zero area"
        @test !isnan(result)
        
        # Equilateral with normal area
        a = b = c = 1.0
        area2 = sqrt(3) / 4  # Area of equilateral triangle
        result = singularF21(a, b, c, area2)
        @test isfinite(result)
        
        # Long thin triangle
        a = 0.01; b = 1.0; c = 1.005
        area2 = 0.0001  # Very small area
        result = singularF21(a, b, c, area2)
        @test isfinite(result)
    end
    
    @testset "Singularities singularF22 - Degenerate Triangles" begin
        using EMSuite.IntegralEquations: singularF22
        
        # Case 1: Near degenerate
        a, b, c = 1.0, 1.0, 1.999
        area2 = 1e-15
        result = singularF22(a, b, c, area2)
        @test isfinite(result) "F22 must be finite despite near-zero area"
        @test !isnan(result)
        
        # Case 2: Equilateral
        a = b = c = 1.0
        area2 = sqrt(3) / 4
        result = singularF22(a, b, c, area2)
        @test isfinite(result)
        
        # Case 3: Long thin
        a = 0.01; b = 1.0; c = 1.005
        area2 = 0.0001
        result = singularF22(a, b, c, area2)
        @test isfinite(result)
    end
    
    # ───────────────────────────────────────────────────────────────────
    # P0.2: FastExp - Boundary at R_max
    # ───────────────────────────────────────────────────────────────────
    
    @testset "FastExp - Boundary Conditions R≈R_max" begin
        using EMSuite.IntegralEquations: FastExpTable, fast_green_func, fast_exp_ikr
        using EMSuite.Core.Constants: c0
        
        freq = 1e8
        k = 2π * freq / c0
        λ = 2π / k
        table = FastExpTable(k; R_max=20*λ)
        
        # Test R just before R_max
        R_before = 0.9999 * table.R_max
        G1 = fast_green_func(table, R_before)
        @test isfinite(real(G1)) && isfinite(imag(G1))
        
        # Test R at R_max
        G2 = fast_green_func(table, table.R_max)
        @test isfinite(real(G2)) && isfinite(imag(G2))
        
        # Test R just after R_max (should use fallback)
        R_after = 1.0001 * table.R_max
        G3 = fast_green_func(table, R_after)
        @test isfinite(real(G3)) && isfinite(imag(G3))
        
        # Verify accuracy at boundary
        G_exact_before = exp(-im * k * R_before) / (4π * R_before)
        rel_err = abs(G1 - G_exact_before) / abs(G_exact_before)
        @test rel_err < 1e-3 "Relative error at R≈R_max should be < 0.1%"
    end
    
    @testset "FastExp - fast_exp_ikr Boundary" begin
        using EMSuite.IntegralEquations: FastExpTable, fast_exp_ikr
        using EMSuite.Core.Constants: c0
        
        freq = 1e8; k = 2π * freq / c0; λ = 2π / k
        table = FastExpTable(k; R_max=20*λ)
        
        # Test at boundaries
        R_at_max = table.R_max
        val1 = fast_exp_ikr(table, R_at_max)
        @test isfinite(real(val1)) && isfinite(imag(val1))
        
        R_beyond = 1.001 * table.R_max
        val2 = fast_exp_ikr(table, R_beyond)
        @test isfinite(real(val2)) && isfinite(imag(val2))
        
        # Test R=0 edge case (should return 0 from fast_green_func, but test separately)
        # Note: fast_exp_ikr doesn't have R<threshold check, so we test fallback path
        R_tiny = 1e-20
        val_tiny = fast_exp_ikr(table, R_tiny)
        @test isfinite(real(val_tiny)) && isfinite(imag(val_tiny))
    end
    
    # ───────────────────────────────────────────────────────────────────
    # P0.3: CoordinateTransforms - acos clamp boundary
    # ───────────────────────────────────────────────────────────────────
    
    @testset "CoordinateTransforms - acos Boundary Values" begin
        using EMSuite.Geometry: r̂θϕInfo
        
        # Test at boundaries and beyond
        test_cases = [
            ("Below -1", -1.0 - 1e-12),
            ("Exactly -1", -1.0),
            ("Just above -1", -1.0 + 1e-12),
            ("Zero (equator)", 0.0),
            ("Just below +1", 1.0 - 1e-12),
            ("Exactly +1", 1.0),
            ("Above +1", 1.0 + 1e-12),
        ]
        
        for (desc, r3_val) in test_cases
            r_vec = [0.0, 0.0, r3_val]
            # Normalize
            r_norm = norm(r_vec)
            if r_norm > 0
                r_vec = r_vec / r_norm
            end
            
            info = r̂θϕInfo(r_vec)
            @test !any(isnan, info.r̂) "r̂ must not have NaN for case: $desc"
            @test !any(isnan, info.θhat) "θhat must not have NaN for case: $desc"
            @test !any(isnan, info.ϕhat) "ϕhat must not have NaN for case: $desc"
            @test all(isfinite, info.r̂)
            @test all(isfinite, info.θhat)
            @test all(isfinite, info.ϕhat)
        end
        
        # Test with realistic vectors that might have rounding errors
        r_vec_numerical = [0.0, 0.0, 1.0 + eps(Float64)]
        info_num = r̂θϕInfo(r_vec_numerical)
        @test all(isfinite, info_num.r̂)
    end
    
    # ───────────────────────────────────────────────────────────────────
    # P1.1: WavePort Parameter Validation
    # ───────────────────────────────────────────────────────────────────
    
    @testset "WavePort - Parameter Validation @test_throws" begin
        using EMSuite.Ports: WavePort, compute_port_modes
        
        # Test invalid a (≤0)
        @test_throws AssertionError WavePort(1; mode=:TE10, a=0.0, b=0.01)
        @test_throws AssertionError WavePort(1; mode=:TE10, a=-0.01, b=0.01)
        
        # Test invalid b (≤0)
        @test_throws AssertionError WavePort(1; mode=:TE10, a=0.01, b=0.0)
        @test_throws AssertionError WavePort(1; mode=:TE10, a=0.01, b=-0.01)
        
        # Test invalid frequency (≤0) - requires compute_port_modes
        port = WavePort(1; mode=:TE10, a=0.02, b=0.01)
        @test_throws AssertionError compute_port_modes(port, 0.0)
        @test_throws AssertionError compute_port_modes(port, -1e9)
        
        # Verify valid parameters still work
        @test_nowarn WavePort(1; mode=:TE10, a=0.02, b=0.01)
        @test_nowarn compute_port_modes(port, 1e9)
    end
    
    # ───────────────────────────────────────────────────────────────────
    # P1.2: CoaxPort Parameter Validation
    # ───────────────────────────────────────────────────────────────────
    
    @testset "CoaxPort - Parameter Validation @test_throws" begin
        using EMSuite.Ports: CoaxPort, coax_impedance
        
        # Test invalid inner_radius (≤0)
        @test_throws AssertionError CoaxPort(1, 0.0, 2e-3, 1.0)
        @test_throws AssertionError CoaxPort(1, -1e-3, 2e-3, 1.0)
        
        # Test invalid outer_radius (≤inner)
        @test_throws AssertionError CoaxPort(1, 2e-3, 1e-3, 1.0)  # outer < inner
        @test_throws AssertionError CoaxPort(1, 1e-3, 1e-3, 1.0)  # outer == inner
        
        # Test invalid eps_r (≤0)
        @test_throws AssertionError CoaxPort(1, 1e-3, 2e-3, 0.0)
        @test_throws AssertionError CoaxPort(1, 1e-3, 2e-3, -1.0)
        
        # Verify valid parameters still work
        port_valid = @test_nowarn CoaxPort(1, 1e-3, 2e-3, 1.0)
        @test_nowarn coax_impedance(port_valid)
    end
    
    # ───────────────────────────────────────────────────────────────────
    # P2: BasisUtilities - Degenerate Geometry Detection
    # ───────────────────────────────────────────────────────────────────
    
    @testset "BasisUtilities - Degenerate Triangle Detection" begin
        using EMSuite: TriangleMesh, RWGBasis
        using EMSuite.BasisFunctions: get_triangle_info
        
        # Case 1: Collinear vertices (zero area)
        # Three points on a line: (0,0,0), (1,0,0), (2,0,0)
        nodes_collinear = [
            0.0 1.0 2.0
            0.0 0.0 0.0
            0.0 0.0 0.0
        ]
        tris_collinear = [1 2 3]'  # Single triangle with collinear vertices
        
        # This should error due to area ≈ 0
        @test_throws Exception begin
            mesh_col = TriangleMesh(1, nodes_collinear, tris_collinear, [1])
            basis_col = RWGBasis(mesh_col)
            get_triangle_info(Float64, nodes_collinear[:, tris_collinear[1, 1]], 
                                      nodes_collinear[:, tris_collinear[2, 1]], 
                                      nodes_collinear[:, tris_collinear[3, 1]])
        end
        
        # Case 2: Duplicate vertices (zero edge length)
        # Two identical points
        nodes_dup = [
            0.0 0.0 1.0
            0.0 0.0 1.0
            0.0 0.0 0.0
        ]
        tris_dup = [1 2 3]'
        
        @test_throws Exception begin
            get_triangle_info(Float64, nodes_dup[:, 1], nodes_dup[:, 2], nodes_dup[:, 3])
        end
        
        # Case 3: Valid triangle (should not error)
        nodes_valid = [
            0.0 1.0 0.0
            0.0 0.0 1.0
            0.0 0.0 0.0
        ]
        @test_nowarn get_triangle_info(Float64, nodes_valid[:, 1], nodes_valid[:, 2], nodes_valid[:, 3])
    end
    
    # ───────────────────────────────────────────────────────────────────
    # 额外：FastExp R=0 验证（确认 Round 6 修复）
    # ───────────────────────────────────────────────────────────────────
    
    @testset "FastExp - R=0 Handling" begin
        using EMSuite.IntegralEquations: FastExpTable, fast_green_func
        using EMSuite.Core.Constants: c0
        
        freq = 1e8; k = 2π * freq / c0
        table = FastExpTable(k)
        
        # Test R=0 (should return 0, not NaN/Inf)
        G_zero = fast_green_func(table, 0.0)
        @test G_zero == 0.0 "FastExp should return 0 at R=0"
        
        # Test very small R (< threshold)
        G_tiny = fast_green_func(table, 1e-12)
        @test G_tiny == 0.0 "FastExp should return 0 for R < threshold"
    end

end

# 测试统计摘要
println("=" ^ 60)
println("Coverage Gaps Test Suite - Summary")
println("=" ^ 60)
println("P0 Tests (Critical Numerical Stability): 5 test sets")
println("  - Singularities F1/F21/F22: 3 test sets")
println("  - FastExp boundary: 2 test sets")
println("  - CoordinateTransforms acos clamp: 1 test set")
println("P1 Tests (Parameter Validation): 2 test sets")
println("  - WavePort: 1 test set")
println("  - CoaxPort: 1 test set")
println("P2 Tests (Degenerate Geometry): 1 test set")
println("Additional Tests (Round 6 verification): 1 test set")
println("Total: 9 test sets covering Round 6-12 fixes")
println("=" ^ 60)
