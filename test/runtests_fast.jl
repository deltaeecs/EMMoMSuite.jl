"""
runtests_fast.jl

快速测试套件 - 用于 CI/CD 和本地快速验证
运行时间目标：< 5 分钟

包含：
- 核心单元测试
- 边界条件测试
- 参数验证测试
- 不包含大规模数值验证
"""

using Test

@testset "EMSuite Fast Tests" begin
    @testset "Core Modules" begin
        include("test_materials.jl")
        include("test_geometry.jl")
        include("test_basis_functions.jl")
    end
    
    @testset "Integral Equations" begin
        include("test_integral_equations.jl")
        include("test_fastexp.jl")
    end
    
    @testset "Coverage Gaps (Phase 21)" begin
        include("test_coverage_gaps.jl")
    end
    
    @testset "Utilities" begin
        include("test_io.jl")
        include("test_lebedev.jl")
    end
    
    @testset "Ports" begin
        include("test_ports.jl")
    end
end
