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

    # ──────────────────────────────────────────────────────────────────────────
    # indices.jl: sizeChunksCuts2indices
    # ──────────────────────────────────────────────────────────────────────────
    @testset "sizeChunksCuts2indices" begin
        # 1D: 10 elements split into 2 chunks
        cuts = MPIArrays.sizeChunks2cuts(10, [2])
        idxs = MPIArrays.sizeChunksCuts2indices((10,), (2,), cuts)
        @test size(idxs) == (2,)
        @test idxs[1] == (1:5,)
        @test idxs[2] == (6:10,)

        # 2D: (6, 4) split into (2, 2) chunks
        cuts2 = MPIArrays.sizeChunks2cuts((6, 4), (2, 2))
        idxs2 = MPIArrays.sizeChunksCuts2indices((6, 4), (2, 2), cuts2)
        @test size(idxs2) == (2, 2)

        # Vector cuts (1D special case)
        cuts_vec = MPIArrays.slicedim2bounds(6, 3)
        idxs3 = MPIArrays.sizeChunksCuts2indices((6,), (3,), cuts_vec)
        @test length(idxs3) == 3
        @test idxs3[1] == (1:2,)
        @test idxs3[3] == (5:6,)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # indices.jl: sizeChunks2idxs
    # ──────────────────────────────────────────────────────────────────────────
    @testset "sizeChunks2idxs" begin
        idxs = MPIArrays.sizeChunks2idxs((10,), (2,))
        @test length(idxs) == 2
        @test idxs[1] == (1:5,)
        @test idxs[2] == (6:10,)

        # 1D integer (non-tuple) overloads
        cuts_a = MPIArrays.sizeChunks2cuts(8, (4,))
        @test length(cuts_a) == 1

        cuts_b = MPIArrays.sizeChunks2cuts((8,), 4)
        @test length(cuts_b) == 1

        cuts_c = MPIArrays.sizeChunks2cuts(8, 4)
        @test length(cuts_c) == 1
    end

    # ──────────────────────────────────────────────────────────────────────────
    # indices.jl: indice2rank
    # ──────────────────────────────────────────────────────────────────────────
    @testset "indice2rank" begin
        # Integer rank2indices (Dict{Integer, NTuple{1}})
        r2i = Dict{Integer,NTuple{1,UnitRange{Int}}}(
            0 => (1:5,),
            1 => (6:10,),
        )
        @test MPIArrays.indice2rank(3, r2i) == 0
        @test MPIArrays.indice2rank(7, r2i) == 1

        # NTuple rank2indices (Dict{Int, Tuple{...}})
        r2i2 = Dict{Int,Tuple{UnitRange{Int},UnitRange{Int}}}(
            0 => (1:4, 1:3),
            1 => (5:8, 1:3),
        )
        @test MPIArrays.indice2rank((2, 2), r2i2) == 0
        @test MPIArrays.indice2rank((6, 2), r2i2) == 1

        # indice2ranks: multiple ranks for range queries
        r2i3 = Dict{Int,Tuple{UnitRange{Int}}}(
            0 => (1:5,),
            1 => (6:10,),
        )
        ranks = MPIArrays.indice2ranks(((1:10,),), r2i3)
        @test 0 in ranks
        @test 1 in ranks
    end

    # ──────────────────────────────────────────────────────────────────────────
    # indices.jl: Base.searchsortedfirst for Vector{UnitRange}
    # ──────────────────────────────────────────────────────────────────────────
    @testset "searchsortedfirst - Vector{UnitRange}" begin
        a = [1:3, 4:7, 8:10]
        # Search for element within first range
        @test Base.searchsortedfirst(a, 2) == 0   # 0 elements fully before range containing 2
        # Search for element within third range
        idx = Base.searchsortedfirst(a, 9)
        @test idx >= 0  # Should find some position
    end

    # ──────────────────────────────────────────────────────────────────────────
    # linearalgebra.jl: axpy! and mul! for MPIArray shims
    # ──────────────────────────────────────────────────────────────────────────
    @testset "MPIArray LinearAlgebra - mul! error dispatch" begin
        # MPIMatrix is not supported for basic mul! with AbstractVector
        # Any MPIMatrix creation requires MPI. We can only test the defined
        # non-MPI helpers here. Test axpy! on ArrayChunk (SubOrMPIArray).
        ac1 = MPIArrays.ArrayChunk(Float64, 1:4)
        fill!(ac1, 2.0)
        ac2 = MPIArrays.ArrayChunk(Float64, 1:4)
        fill!(ac2, 1.0)

        # axpy!: Y += a * X → [1+2*2, 1+2*2, 1+2*2, 1+2*2] = [5,5,5,5]
        LinearAlgebra.axpy!(2.0, ac1, ac2)
        @test all(==(5.0), ac2.data)
    end

end
