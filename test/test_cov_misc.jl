# test_cov_misc.jl — 杂项模块轻量覆盖
#
# 覆盖：
#   Utilities/Logging.jl（init_logging）
#   PostProcessing/RCSBatch.jl（RCSResult 频率响应与角度切面）
#   FastAlgorithms/Lebedev/pinv2interpW.jl（插值权重计算）
#   Geometry/HexMeshIO.jl（read_hex_mesh 错误路径 + validate_mesh）
using Test
using EMMoMSuite
using EMMoMSuite.Utilities
using EMMoMSuite.Geometry
using EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW: interpWeightsInitial, pinv2W!
using LinearAlgebra, Logging, Random, SparseArrays

@testset "Logging: init_logging 文件输出" begin
    logpath = tempname() * ".log"
    init_logging(logpath)
    @info "coverage logging smoke test"
    @test isfile(logpath)
    global_logger(ConsoleLogger(stderr, Logging.Info))
end

@testset "RCSBatch: 频率响应与角度切面" begin
    freqs = [1.0e9, 2.0e9]
    ths = [0.0, 90.0, 180.0]
    phs = [0.0, 90.0]
    base = reshape(collect(1.0:12.0), 2, 3, 2)
    result = RCSResult(freqs, ths, phs, base, base .+ 1.0, base .+ 2.0, base .+ 3.0)

    f, r = rcs_frequency_response(result; theta_idx = 2, phi_idx = 1, polarization = :VV)
    @test f == freqs
    @test length(r) == 2
    @test all(isfinite, r)
    @test r[1] ≈ 10 * log10(3.0)  # base[1,2,1] = 3.0

    a, rcs = rcs_angular_pattern(result; freq_idx = 2, polarization = :HH, cut_plane = :E)
    @test a == ths
    @test length(rcs) == 3
    a2, rcs2 = rcs_angular_pattern(result; cut_plane = :H)
    @test a2 == phs
    @test length(rcs2) == 2

    @test_throws ArgumentError rcs_frequency_response(result; theta_idx = 99)
    @test_throws ArgumentError rcs_angular_pattern(result; cut_plane = :X)
    @test_throws ArgumentError rcs_frequency_response(result; polarization = :ZZ)
end

@testset "pinv2interpW: 插值权重" begin
    Random.seed!(11)
    tN = normalize(randn(3, 4))
    pN = normalize(randn(3, 5))
    w0 = interpWeightsInitial(tN, pN; nInterp = 2)
    @test size(w0) == (10, 8)
    @test all(isfinite, nonzeros(w0))
    @test all(vec(sum(abs2, w0; dims = 2)) .> 0)
    @test all(vec(sum(abs2, w0; dims = 1)) .> 0)

    # pinv2W!：逐行分支（nInterp < size(w, 2) ÷ 2）
    xx = randn(8, 20)
    yy = randn(10, 20)
    w1 = pinv2W!(copy(w0), 2, xx, yy)
    @test size(w1) == (10, 8)
    @test all(isfinite, nonzeros(w1))
    # pinv2W!：满矩阵分支（nInterp >= size(w, 2) ÷ 2）
    w2 = pinv2W!(copy(w0), 5, xx, yy)
    @test size(w2) == (10, 8)
    @test all(isfinite, nonzeros(w2))
end

@testset "HexMeshIO: read_hex_mesh 与 validate_mesh" begin
    # 错误路径：不支持的格式
    @test_throws ErrorException read_hex_mesh("whatever.msh"; format = :vtu)
    # 错误路径：非六面体网格（三角 .msh 走 read_msh_mesh 后类型检查失败）
    msh = tempname() * ".msh"
    open(msh, "w") do io
        println(io, "\$MeshFormat")
        println(io, "4.1 0 8")
        println(io, "\$EndMeshFormat")
        println(io, "\$Nodes")
        println(io, "1 3 1 3")
        println(io, "2 1 0 3")
        println(io, "1")
        println(io, "2")
        println(io, "3")
        println(io, "0.0 0.0 0.0")
        println(io, "1.0 0.0 0.0")
        println(io, "0.0 1.0 0.0")
        println(io, "\$EndNodes")
        println(io, "\$Elements")
        println(io, "1 1 1 1")
        println(io, "2 1 2 1")
        println(io, "1 1 2 3")
        println(io, "\$EndElements")
    end
    try
        @test_throws ErrorException read_hex_mesh(msh)
    finally
        rm(msh; force = true)
    end

    # validate_mesh：合法六面体
    nodes = [
        0.0 1.0 1.0 0.0 0.0 1.0 1.0 0.0;
        0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
        0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
    ]
    hm = HexahedraMesh(1, nodes, reshape(collect(1:8), 8, 1), [1])
    @test validate_mesh(hm)
    # 空网格（0 单元）应失败
    hm0 = HexahedraMesh(1, zeros(3, 8), zeros(Int, 8, 0), Int[])
    @test !validate_mesh(hm0)
    # 倒置四面体（[A,C,B,D] 使体积为负）经 mesh_quality 应失败
    nodes4 = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    tetbad = TetrahedraMesh(1, nodes4, reshape([1, 3, 2, 4], 4, 1), [1])
    @test !validate_mesh(tetbad)
end
