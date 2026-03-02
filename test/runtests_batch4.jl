# 分批覆盖率测试脚本 - 第4批: Lebedev + 后处理 + IO + 集成
using Test
using EMSuite

include("test_lebedev.jl")
include("test_postprocessing.jl")
include("test_io.jl")
include("test_integration.jl")
