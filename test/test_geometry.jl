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

        # TetrahedraMesh with default tags
        mesh_tet2 = TetrahedraMesh(1, nodes_tet, elements_tet)
        @test mesh_tet2.tags == [0]
    end

    @testset "NAS - CTETRA and CHEXA reading" begin
        # Test read_nas_mesh with CTETRA (fixed width)
        nas_ctetra = """
\$ Test NAS file with tetrahedra (fixed width)
GRID    1               0.0     0.0     0.0
GRID    2               1.0     0.0     0.0
GRID    3               0.0     1.0     0.0
GRID    4               0.0     0.0     1.0
CTETRA  1       100     1       2       3       4
"""
        path = "test_ctetra.nas"
        open(path, "w") do f
            write(f, nas_ctetra)
        end
        try
            mesh = read_nas_mesh(path)
            @test mesh isa TetrahedraMesh
            @test num_vertices(mesh) == 4
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 100
        finally
            rm(path, force=true)
        end

        # Test read_nas_mesh with CHEXA (free-field format)
        nas_chexa_ff = """
\$ Test NAS file with hexahedra (free field)
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 1.0, 1.0, 0.0
GRID, 4, , 0.0, 1.0, 0.0
GRID, 5, , 0.0, 0.0, 1.0
GRID, 6, , 1.0, 0.0, 1.0
GRID, 7, , 1.0, 1.0, 1.0
GRID, 8, , 0.0, 1.0, 1.0
CHEXA, 1, 200, 1, 2, 3, 4, 5, 6,
+, 7, 8
"""
        path2 = "test_chexa.nas"
        open(path2, "w") do f
            write(f, nas_chexa_ff)
        end
        try
            mesh = read_nas_mesh(path2)
            @test mesh isa HexahedraMesh
            @test num_vertices(mesh) == 8
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 200
        finally
            rm(path2, force=true)
        end

        # Test read_mixed_nas_mesh
        nas_mixed = """
\$ Mixed mesh: triangles + tetrahedra
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 0.0, 1.0, 0.0
GRID, 4, , 0.0, 0.0, 1.0
CTRIA3, 1, 100, 1, 2, 3
CTETRA, 1, 200, 1, 2, 3, 4
"""
        path3 = "test_mixed.nas"
        open(path3, "w") do f
            write(f, nas_mixed)
        end
        try
            surf_mesh, vol_mesh, hexa_mesh = read_mixed_nas_mesh(path3)
            @test surf_mesh isa TriangleMesh
            @test num_elements(surf_mesh) == 1
            @test vol_mesh isa TetrahedraMesh
            @test num_elements(vol_mesh) == 1
            @test hexa_mesh isa HexahedraMesh
            @test num_elements(hexa_mesh) == 0
        finally
            rm(path3, force=true)
        end
    end

    @testset "NAS - parse_nastran_float and GRID* format" begin
        using EMSuite.Geometry: parse_nastran_float

        # Normal floats
        @test parse_nastran_float("1.5", Float64) ≈ 1.5
        @test parse_nastran_float("0", Float64) ≈ 0.0
        @test parse_nastran_float("", Float64) ≈ 0.0

        # Nastran compressed scientific notation
        @test parse_nastran_float("1.78-15", Float64) ≈ 1.78e-15
        @test parse_nastran_float("3.2+5", Float64) ≈ 3.2e5
        @test parse_nastran_float("-1.5-3", Float64) ≈ -1.5e-3
    end

    @testset "NAS - write_nas_mesh" begin
        # Create a mesh, write it, read it back
        nodes = [
            0.0 1.0 0.0 1.0;
            0.0 0.0 1.0 1.0;
            0.0 0.0 0.0 0.0
        ]
        tris = [1 2; 2 4; 3 3]
        tags = [100, 200]
        mesh_orig = TriangleMesh(2, nodes, tris, tags)

        path = "test_write.nas"
        write_nas_mesh(path, mesh_orig)
        try
            @test isfile(path)
            mesh_read = read_nas_mesh(path)
            @test num_elements(mesh_read) == 2
            @test num_vertices(mesh_read) == 4
        finally
            rm(path, force=true)
        end
    end

    @testset "MSH - Gmsh tetrahedra reading" begin
        # Create Gmsh .msh file with tetrahedra (element type 4)
        # Gmsh tetrahedral type = 4
        msh_tet = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
\$EndNodes
\$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
\$EndElements
"""
        path = "test_tet.msh"
        open(path, "w") do f
            write(f, msh_tet)
        end
        try
            mesh = read_msh_mesh(path)
            @test mesh isa TetrahedraMesh
            @test num_vertices(mesh) == 4
            @test num_elements(mesh) == 1
        finally
            rm(path, force=true)
        end
    end
end
