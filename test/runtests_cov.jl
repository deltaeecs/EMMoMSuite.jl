# runtests_cov.jl — 覆盖率采集套件
#
# = light_cov 全部用例 + MPI/Parallel 补测（test_mpi_coverage.jl）
# 运行（需采集 src 源码覆盖率）：
#   julia --project=. --code-coverage=user test/runtests_cov.jl
# 注意：Windows/Julia 1.12 若复用了未插桩的编译缓存，将不产生 .cov；
#   可靠做法是先移开 ~/.julia/compiled/v1.12/EMMoMSuite 再运行。
#
# 时间门（实测，Windows Julia 1.12）：
#   - 常规运行（无覆盖率插桩，单线程）：≈ 2m10s ≤ 3 分钟
#   - 覆盖率采集（Julia 1.12 需从源码带插桩重编译）：≈ 10 分钟
using Test

@testset "EMMoMSuite coverage suite" begin
    include("runtests_light_cov.jl")
    include("test_mpi_coverage.jl")
    include("test_cov_small.jl")
    include("test_cov_ports.jl")
    include("test_cov_geometry.jl")
    include("test_cov_post.jl")
    include("test_pmchw.jl")
    include("test_nmuller.jl")
    include("test_pmchw_operator_shell.jl")
    include("test_mlfma.jl")
    include("test_aim_operator.jl")
    include("test_cov_basis.jl")
    include("test_preconditioners.jl")
    include("test_cov_solvers.jl")
    include("test_cov_core.jl")
    include("test_cov_scfie.jl")
    include("test_cov_vefie.jl")
    include("test_cov_io.jl")
    include("test_pmchw_excitation.jl")
    include("test_scfie_delta_gap.jl")
    include("test_cov_ie.jl")
    include("test_cov_abs.jl")
    include("test_cov_fa.jl")
    include("test_cov_post2.jl")
    include("test_cov_volume_assembly.jl")
    include("test_cov_accuracy.jl")
    include("test_cov_misc.jl")
end
