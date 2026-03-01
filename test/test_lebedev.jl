using EMSuite
using EMSuite.FastAlgorithms.Lebedev
using Test
using StaticArrays

@testset "Lebedev" begin
    # Debug info
    dict = EMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints.lbnP2FILEDict
    println("Test Dict keys: ", sort(collect(keys(dict))))
    println("Test Dict[3]: ", get(dict, 3, "missing"))

    # Test getlbSortedData
    # p=3 corresponds to 6 points usually, let's check
    # 3.6.txt exists in nodesSorted
    nodes, weights = getlbSortedData(3)
    println("Nodes size: ", size(nodes))

    @test size(nodes, 1) == 3
    @test length(weights) == size(nodes, 2)
    @test size(nodes, 2) == 6 # 3.6.txt

    # Test get_t_nodes
    # t=1 -> p=3
    nodes_t = get_t_nodes(1)
    println("Nodes_t size: ", size(nodes_t))
    @test size(nodes_t, 1) == 3
    @test size(nodes_t, 2) == 6

    # Test nodes2Poles
    poles = nodes2Poles(nodes)
    @test length(poles) == size(nodes, 2)
    @test poles[1].r̂ isa StaticArrays.SVector

    # Test LVI
    @testset "LVI" begin
        # levelIntegralInfoCal
        # We use a large cube size to trigger Lebedev? Or small?
        # truncationLCal(1.0) -> L.
        # If 2L+1 < max(keys(p2nDict)), it uses Lebedev.
        # max p is 131. L ~ 65.
        # L ~ k*d + ...
        # If d=1.0 (wavelengths?), L might be small enough.

        # We need to import levelIntegralInfoCal
        using EMSuite.FastAlgorithms.MLFMA.Interpolation: levelIntegralInfoCal

        # Call it
        # Note: It might fail if dependencies are missing or paths are wrong.
        # But let's try.
        truncL, poles = levelIntegralInfoCal(0.5, Val(:LbTrained1Step))
        println("TruncL: ", truncL)
        @test truncL isa Integer
        @test poles isa LbPolesInfo
    end
end
