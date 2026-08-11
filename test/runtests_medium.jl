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

@testset "EMMoMSuite Medium Tests" begin
    @testset "Solvers" begin
        include("test_solvers.jl")
        include("test_solvers_verification.jl")
        include("test_preconditioners.jl")
    end
    
    @testset "MLFMA" begin
        include("test_mlfma.jl")
        include("test_mlfma_params.jl")
    end

    @testset "Fast Solvers (ACA/MLACA)" begin
        include("test_aca.jl")
        include("test_block_evaluator.jl")
        include("test_aca_operator.jl")
        include("test_mlaca_operator.jl")
        include("test_lowfreq_aca.jl")
        include("test_block_lu.jl")
        include("test_block_evaluator_pmchw.jl")
        include("test_aca_mlaca_pmchw.jl")
    end
    
    @testset "Parallel" begin
        include("test_parallel.jl")
        include("test_mpiarray.jl")
        include("test_parallel_mfie_cfie.jl")
    end
    
    @testset "Integral Equations - Advanced" begin
        include("test_mfie_decomposition.jl")
        include("test_mfie_near_quadrature.jl")
        include("test_integral_equations_endtoend.jl")  # Phase 21.2
        include("test_integral_equations_frequency.jl")  # Phase 21.3
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
