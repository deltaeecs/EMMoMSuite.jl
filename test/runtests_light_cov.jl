# runtests_light_cov.jl
# 轻量覆盖率测试套件，跳过 MLFMA 和大型装配（覆盖率 instrumentation 下过慢）。
# 可通过以下方式运行:
#   julia --project=EMMoMSuite --code-coverage=user --startup-file=no EMMoMSuite/test/runtests_light_cov.jl
#
# 跳过的模块:
#   - test_mlfma.jl              (MLFMA近场装配, ~5M核函数调用, coverage下极慢)
#   - test_parallel.jl           (N=932 SWG MPI装配)
#   - test_parallel_mfie_cfie.jl (CFIE MPI装配)
#   - test_distributed_gmres.jl  (DistGMRES大系统)
#   - test_volume_assembly_mpi.jl (内含N=932 SWG装配)
#   - test_scfie.jl              (SCFIE非并行装配)
#   - test_pwc.jl                (PWC非并行装配, 含大型四面体网格)
#   - test_hex_rbf.jl            (Hex/RBF非并行装配)
#   - test_workflow.jl           (全链路仿真驱动器)
# 以上均在 run_tests.jl 无覆盖率模式下已验证通过（484 tests / 3 分钟）。

using Test

@testset "EMMoMSuite (light coverage)" begin
    include("test_materials.jl")
    include("test_geometry.jl")
    include("test_basis_functions.jl")
    include("test_integral_equations.jl")
    include("test_fastexp.jl")
    # SKIP: test_mfie_decomposition.jl (密集 CFIE 分解验证, ~1330 RWG 装配, 覆盖率插桩下过慢)
    include("test_solvers.jl")
    include("test_solvers_verification.jl")
    # SKIP: test_preconditioners.jl (BlockJacobi 立方块 LU + MLFMA 栈加载, 覆盖率插桩下过慢)
    # SKIP: test_mlfma.jl
    # SKIP: test_parallel.jl
    include("test_mpiarray.jl")
    # SKIP: test_parallel_mfie_cfie.jl
    # SKIP: test_distributed_gmres.jl
    # SKIP: test_volume_assembly_mpi.jl
    include("test_mpi_array_utils.jl")
    include("test_lebedev.jl")
    include("test_lebedev_interp.jl")
    include("test_postprocessing.jl")
    include("test_io.jl")
    include("test_integration.jl")
    # SKIP: test_workflow.jl
    # SKIP: test_scfie.jl
    # SKIP: test_pwc.jl
    # SKIP: test_hex_rbf.jl
end
