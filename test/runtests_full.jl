"""
runtests_full.jl

完整测试套件 - 包含所有测试
运行时间目标：< 2 小时

包含：
- 所有单元测试
- 所有集成测试
- 大规模数值验证
- Legacy 对齐测试
- 长时间运行的精度测试
"""

using Test

@testset "EMMoMSuite Full Test Suite" begin
    # Fast tests
    println("\n" * "="^60)
    println("Running FAST tests...")
    println("="^60)
    include("runtests_fast.jl")
    
    # Medium tests
    println("\n" * "="^60)
    println("Running MEDIUM tests...")
    println("="^60)
    include("runtests_medium.jl")
    
    # Slow tests - large scale verification
    println("\n" * "="^60)
    println("Running SLOW tests...")
    println("="^60)
    
    @testset "Integration & Workflow" begin
        include("test_integration.jl")
        include("test_workflow.jl")
    end
    
    @testset "PMCHW - Comprehensive" begin
        include("test_pmchw_excitation.jl")
        include("test_pmchw.jl")
        include("test_pmchw_operator_shell.jl")
        include("test_pmchw_operator_shell_mlfma.jl")
        include("test_pmchw_gate_s_dense.jl")
    end
    
    @testset "N-Müller - Comprehensive" begin
        include("test_nmuller.jl")
        include("test_nmuller_excitation.jl")
        include("test_nmuller_comparison.jl")
    end
    
    @testset "Special Cases" begin
        include("test_scfie_delta_gap.jl")
        include("test_hex_rbf.jl")
        include("test_feko_reader.jl")
    end
    
    @testset "Advanced Post Processing" begin
        include("test_postprocessing_advanced.jl")
        include("test_accuracy_metrics.jl")
        include("test_benchmark_report_data.jl")
    end
    
    @testset "Release & Legacy" begin
        include("test_release_workflow.jl")
        include("test_legacy_parity.jl")
    end
    
    @testset "MPI Tests (if available)" begin
        if haskey(ENV, "JULIA_MPI_TEST_ENABLE")
            include("test_distributed_gmres.jl")
            include("test_volume_assembly_mpi.jl")
            include("test_mpi_array_utils.jl")
        else
            @info "Skipping MPI tests (set JULIA_MPI_TEST_ENABLE=1 to enable)"
        end
    end
end

# Test summary
println("\n" * "="^60)
println("FULL TEST SUITE COMPLETED")
println("="^60)
