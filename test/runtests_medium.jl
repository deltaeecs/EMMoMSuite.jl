"""
runtests_medium.jl

中等规模测试套件
运行时间目标：< 30 分钟

包含：
- 集成测试
- MLFMA 测试
- 并行测试
- 中等规模数值验证
"""

using Test

@testset "EMSuite Medium Tests" begin
    @testset "Solvers" begin
        include("test_solvers.jl")
        include("test_solvers_verification.jl")
        include("test_preconditioners.jl")
    end
    
    @testset "MLFMA" begin
        include("test_mlfma.jl")
    end
    
    @testset "Parallel" begin
        include("test_parallel.jl")
        include("test_mpiarray.jl")
        include("test_parallel_mfie_cfie.jl")
    end
    
    @testset "Integral Equations - Advanced" begin
        include("test_mfie_decomposition.jl")
        include("test_mfie_near_quadrature.jl")
    end
    
    @testset "Volume Methods" begin
        include("test_scfie.jl")
        include("test_pwc.jl")
    end
    
    @testset "Post Processing" begin
        include("test_postprocessing.jl")
        include("test_antenna_feeds.jl")
    end
end
