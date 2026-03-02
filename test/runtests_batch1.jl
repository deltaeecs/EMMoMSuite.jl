# 分批覆盖率测试脚本 - 第1批: 基础测试（快速）
# 从 EMSuite/test 目录运行: julia --project=../../EMSuite --code-coverage=user runtests_batch1.jl
using Test
using EMSuite

# 批1: 基础/几何/积分方程
include("test_materials.jl")
include("test_geometry.jl")
include("test_basis_functions.jl")
include("test_integral_equations.jl")
include("test_fastexp.jl")
include("test_mfie_decomposition.jl")
include("test_solvers.jl")
include("test_solvers_verification.jl")
include("test_preconditioners.jl")
