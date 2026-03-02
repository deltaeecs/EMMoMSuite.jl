# test_mpi_array_utils.jl
#
# 测试 MPIArray 工具函数（单进程可执行部分）：
# - arraychunk.jl: ArrayChunk 结构及索引操作
# - indices.jl: slicedim2bounds, slicedim2partition, expandslice,
#               sizeChunks2cuts
#
# 这些函数不依赖多进程 MPI，单进程即可测试。
using Test
using EMSuite
using EMSuite.Parallel: MPIArrays

@testset "MPIArray Utils" begin

    # ──────────────────────────────────────────────────────────────────────────
    # arraychunk.jl
    # ──────────────────────────────────────────────────────────────────────────
    @testset "ArrayChunk - Construction and Interface" begin
        # 1D ArrayChunk
        ac1d = MPIArrays.ArrayChunk(Float64, 2:5)
        @test size(ac1d) == (4,)
        @test length(ac1d) == 4
        @test eltype(ac1d) == Float64

        # Fill and sum
        fill!(ac1d, 3.0)
        @test sum(ac1d) ≈ 12.0

        # getindex: index within range → value
        ac1d[3] = 99.0
        @test ac1d[3] ≈ 99.0

        # getindex: index outside range → zero
        @test ac1d[10] == 0.0

        # 2D ArrayChunk from data array
        data2d = rand(3, 4)
        ac2d = MPIArrays.ArrayChunk(data2d, 1:3, 5:8)
        @test size(ac2d) == (3, 4)
        @test length(ac2d) == 12

        # setindex! within range
        ac2d[2, 6] = 99.0
        @test ac2d[2, 6] ≈ 99.0
    end

    # ──────────────────────────────────────────────────────────────────────────
    # indices.jl: slicedim2bounds
    # ──────────────────────────────────────────────────────────────────────────
    @testset "slicedim2bounds" begin
        # sz=10, nc=3 → chunks of sizes ~4,3,3
        g = MPIArrays.slicedim2bounds(10, 3)
        @test length(g) == 4
        # Grid must start at 1 and end at 11 (exclusive upper)
        @test g[1] == 1
        @test g[end] == 11
        # Consecutive differences give chunk sizes
        sizes = diff(g)
        @test sum(sizes) == 10
        @test all(s -> s > 0, sizes)

        # sz < nc case: should return [[1:sz+1...]; zeros(nc-sz)]
        g2 = MPIArrays.slicedim2bounds(2, 5)
        @test length(g2) == 6
    end

    @testset "slicedim2partition" begin
        # 2D: nc=4 → factors=(2,2) for balanced dims
        parts = MPIArrays.slicedim2partition((10, 10), 4)
        @test prod(parts) == 4

        # 1D: nc=8 → only 1 chunk in dim
        parts1d = MPIArrays.slicedim2partition([16], 8)
        @test prod(parts1d) == 8

        # nc=1 → all ones
        parts_1 = MPIArrays.slicedim2partition((5, 5), 1)
        @test all(==(1), parts_1)
    end

    @testset "sizeChunks2cuts" begin
        # (10,) split into 2 chunks
        cuts = MPIArrays.sizeChunks2cuts(10, [2])
        @test length(cuts) == 1
        @test length(cuts[1]) == 3  # 3 boundary points for 2 chunks

        # (10, 8) split into (2, 2) chunks
        cuts2 = MPIArrays.sizeChunks2cuts((10, 8), (2, 2))
        @test length(cuts2) == 2
    end

    @testset "expandslice" begin
        # Expand 3:7 by 2 within 1:10
        expanded = MPIArrays.expandslice(3:7, 2, 1:10)
        @test expanded == 1:9

        # Expansion clamped by bounds
        expanded2 = MPIArrays.expandslice(1:3, 5, 1:10)
        @test first(expanded2) == 1
        @test last(expanded2) == 8
    end

end
