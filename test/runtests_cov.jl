# runtests_cov.jl — 覆盖率采集套件（单线程 ≤3 分钟目标）
#
# = light_cov 全部用例 + MPI/Parallel 补测（test_mpi_coverage.jl）
# 运行（需采集 src 源码覆盖率）：
#   julia --project=. --code-coverage=user --pkgimages=no test/runtests_cov.jl
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
end
