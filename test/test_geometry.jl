using Test
using EMSuite
using EMSuite.Geometry

@testset "Geometry" begin
    @testset "TriangleMesh" begin
        # Create a simple mesh manually
        # 2 triangles forming a square
        # (0,1) -- (1,1)
        #   | \      |
        #   |  \     |
        # (0,0) -- (1,0)

        nodes = [
            0.0 1.0 0.0 1.0
            0.0 0.0 1.0 1.0
            0.0 0.0 0.0 0.0
        ]

        # Triangle 1: 1-2-3
        # Triangle 2: 2-4-3
        tris = [
            1 2
            2 4
            3 3
        ]

        tags = [1, 2] # Different tags for testing

        mesh = TriangleMesh(2, nodes, tris, tags)

        @test num_vertices(mesh) == 4
        @test num_elements(mesh) == 2
        @test dimension(mesh) == 2
        @test mesh.tags == [1, 2]
        @test vertices(mesh) == nodes
        @test elements(mesh) == tris
    end

    @testset "MeshIO (Mock NAS)" begin
        # Create a temporary NAS file
        nas_content = """
\$ Test NAS file
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 0.0, 1.0, 0.0
CTRIA3, 1, 100, 1, 2, 3
"""
        path = "test_mesh.nas"
        open(path, "w") do f
            write(f, nas_content)
        end

        try
            mesh = read_nas_mesh(path)
            @test num_vertices(mesh) == 3
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 100
            @test mesh.node[:, 1] ≈ [0.0, 0.0, 0.0]
        finally
            rm(path, force = true)
        end
    end

    @testset "GmshIO (Mock MSH)" begin
        # Create a temporary MSH file (Gmsh 4.1)
        msh_content = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Nodes
1 3 1 3
2 1 0 3
1
2
3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
\$EndNodes
\$Elements
1 1 1 1
2 1 2 1
1 1 2 3
\$EndElements
"""
        # Explanation of Mock Data:
        # Nodes: 1 block, 3 nodes.
        # Block 1: dim=2, tag=1, parametric=0, num=3
        # Node tags: 1, 2, 3
        # Coords: (0,0,0), (1,0,0), (0,1,0)

        # Elements: 1 block, 1 element.
        # Block 1: dim=2, tag=1 (Physical Tag), type=2 (Triangle), num=1
        # Element 1: tag=1, nodes=1, 2, 3

        path = "test_mesh.msh"
        open(path, "w") do f
            write(f, msh_content)
        end

        try
            mesh = read_msh_mesh(path)
            @test num_vertices(mesh) == 3
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 1
            @test mesh.node[:, 1] ≈ [0.0, 0.0, 0.0]
        finally
            rm(path, force = true)
        end
    end

    @testset "GaussQuadrature - Hexahedron and Quadrangle" begin
        using EMSuite.Geometry: GaussQuadratureInfo, gaussQuadratureHexa1D

        # Test Hexahedron GQ (8 = 2³ points)
        gq_hexa = GaussQuadratureInfo(:Hexahedron, 8)
        @test length(gq_hexa.weight) == 8
        @test isapprox(sum(gq_hexa.weight), 1.0; rtol=1e-10)
        @test size(gq_hexa.coordinate, 1) == 8  # 8 shape functions
        @test size(gq_hexa.coordinate, 2) == 8  # 8 points

        # Test Hexahedron GQ (1 = 1³ point)
        gq_hexa1 = GaussQuadratureInfo(:Hexahedron, 1)
        @test length(gq_hexa1.weight) == 1
        @test isapprox(sum(gq_hexa1.weight), 1.0; rtol=1e-10)

        # Test Quadrangle GQ (4 = 2² points)
        gq_quad = GaussQuadratureInfo(:Quadrangle, 4)
        @test length(gq_quad.weight) == 4
        @test isapprox(sum(gq_quad.weight), 1.0; rtol=1e-10)
        @test size(gq_quad.coordinate, 1) == 4  # 4 shape functions

        # Test Quadrangle GQ (1 point)
        gq_quad1 = GaussQuadratureInfo(:Quadrangle, 1)
        @test length(gq_quad1.weight) == 1
        @test isapprox(sum(gq_quad1.weight), 1.0; rtol=1e-10)

        # Test gaussQuadratureHexa1D
        x, w = gaussQuadratureHexa1D(2)
        @test length(x) == 2
        @test length(w) == 2
        @test all(0.0 .<= x .<= 1.0)
        @test isapprox(sum(w), 1.0; rtol=1e-10)

        # Test error path
        @test_throws ErrorException GaussQuadratureInfo(:InvalidGeometry, 3)
    end

    @testset "GaussQuadrature - MeshTypes" begin
        using EMSuite.Geometry: TriangleMesh, TetrahedraMesh
        # TetrahedraMesh construction
        nodes_tet = [0.0 1.0 0.0 0.0;
                     0.0 0.0 1.0 0.0;
                     0.0 0.0 0.0 1.0]
        elements_tet = reshape([1; 2; 3; 4], 4, 1)
        tags_tet = [1]
        mesh_tet = TetrahedraMesh(1, nodes_tet, elements_tet, tags_tet)
        @test num_elements(mesh_tet) == 1
        @test num_vertices(mesh_tet) == 4
        @test dimension(mesh_tet) == 3
    end
end
