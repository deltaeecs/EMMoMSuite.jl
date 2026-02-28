using Test
using EMSuite
using StaticArrays
using LinearAlgebra

@testset "RWG Basis Functions" begin
    # Create a simple mesh with 2 triangles sharing an edge
    # Nodes:
    # 1: (0,0,0)
    # 2: (1,0,0)
    # 3: (0,1,0)
    # 4: (1,1,0)
    
    nodes = [
        0.0 1.0 0.0 1.0;
        0.0 0.0 1.0 1.0;
        0.0 0.0 0.0 0.0
    ]
    
    # Elements (Triangles)
    # T1: 1-2-3
    # T2: 2-4-3 (Note: 2-3 is the shared edge)
    elements = [
        1 2;
        2 4;
        3 3
    ]
    
    # Tags (Region 1)
    tags = [1, 1]
    
    mesh = TriangleMesh(1, nodes, elements, tags)
    
    @test num_vertices(mesh) == 4
    @test num_elements(mesh) == 2
    
    # Create RWG basis
    basis = RWGBasis(mesh)
    
    # We expect 1 internal edge (2-3)
    # Edges of T1: 1-2, 2-3, 3-1
    # Edges of T2: 2-4, 4-3, 3-2
    # Shared: 2-3
    
    @test num_basis(basis) == 1
    
    rwg = basis.functions[1]
    @test rwg.is_boundary == false
    
    # Check support
    # Should be T1 and T2 (indices 1 and 2)
    # Order depends on sorting, but it should contain 1 and 2
    @test (rwg.support[1] == 1 && rwg.support[2] == 2) || (rwg.support[1] == 2 && rwg.support[2] == 1)
    
    # Check edge length
    # Edge 2-3: (1,0,0) to (0,1,0) -> length sqrt(2)
    @test isapprox(rwg.edge_length, sqrt(2.0))
    
    # Check center
    # Midpoint of (1,0,0) and (0,1,0) -> (0.5, 0.5, 0)
    @test isapprox(rwg.center, SVector(0.5, 0.5, 0.0))
end

@testset "SWG Basis Functions" begin
    # Create a simple mesh with 2 tetrahedra sharing a face
    # Nodes:
    # 1: (0,0,0)
    # 2: (1,0,0)
    # 3: (0,1,0)
    # 4: (0,0,1)
    # 5: (1,1,1)
    
    nodes = [
        0.0 1.0 0.0 0.0 1.0;
        0.0 0.0 1.0 0.0 1.0;
        0.0 0.0 0.0 1.0 1.0
    ]
    
    # Elements (Tetrahedra)
    # T1: 1-2-3-4
    # T2: 2-3-4-5 (Shared face: 2-3-4)
    elements = [
        1 2;
        2 3;
        3 4;
        4 5
    ]
    
    tags = [1, 1]
    
    mesh = TetrahedraMesh(2, nodes, elements, tags)
    
    @test num_vertices(mesh) == 5
    @test num_elements(mesh) == 2
    
    # Create SWG basis
    basis = SWGBasis(mesh)
    
    # SWG creates basis functions for all faces (boundary + internal)
    # Faces of T1: 1-2-3, 1-2-4, 1-3-4, 2-3-4  (3 boundary + 1 shared)
    # Faces of T2: 2-3-4, 2-3-5, 2-4-5, 3-4-5  (3 boundary + 1 shared)
    # Total: 6 boundary + 1 internal = 7
    
    @test num_basis(basis) == 7
    
    # Find internal (non-boundary) basis function
    internal_idx = findfirst(f -> !f.is_boundary, basis.functions)
    @test internal_idx !== nothing
    
    swg = basis.functions[internal_idx]
    @test swg.is_boundary == false
    
    # Check support
    @test (swg.support[1] == 1 && swg.support[2] == 2) || (swg.support[1] == 2 && swg.support[2] == 1)
    
    # Check area of face 2-3-4
    # v2=(1,0,0), v3=(0,1,0), v4=(0,0,1)
    # Side lengths: sqrt(2), sqrt(2), sqrt(2) -> Equilateral triangle
    # Area = sqrt(3)/4 * a^2 = sqrt(3)/4 * 2 = sqrt(3)/2
    @test isapprox(swg.area, sqrt(3)/2)
end

@testset "PWC Basis Functions" begin
    # Single tetrahedron
    nodes = [
        0.0 1.0 0.0 0.0;
        0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 1.0
    ]
    elements = [1; 2; 3; 4] # 4x1 matrix
    elements = reshape(elements, 4, 1)
    
    tags = [1]
    
    mesh = TetrahedraMesh(1, nodes, elements, tags)
    
    basis = PWCBasis(mesh)
    
    # 3 DOFs per tetrahedron (x, y, z components)
    @test num_basis(basis) == 3
    
    pwc = basis.functions[1]
    @test pwc.support == 1
    
    # Volume of unit tetrahedron = 1/6
    @test isapprox(pwc.volume, 1.0/6.0)
    
    # Check inBfsID mapping: tet 1 → global IDs 1,2,3
    @test pwc.inBfsID == SVector(1, 2, 3)
end

@testset "RBF Basis Functions" begin
    # Two hexahedra sharing a face
    # Hex 1: Unit cube at origin
    # Hex 2: Unit cube shifted by (1,0,0)
    
    # Nodes
    # 1-8: Hex 1
    # 1:(0,0,0), 2:(1,0,0), 3:(1,1,0), 4:(0,1,0) (Bottom z=0)
    # 5:(0,0,1), 6:(1,0,1), 7:(1,1,1), 8:(0,1,1) (Top z=1)
    
    # 9-12: Extra nodes for Hex 2 (sharing 2,3,6,7)
    # 9:(2,0,0), 10:(2,1,0)
    # 11:(2,0,1), 12:(2,1,1)
    
    nodes = [
        0.0 1.0 1.0 0.0 0.0 1.0 1.0 0.0 2.0 2.0 2.0 2.0;
        0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 0.0 1.0 0.0 1.0;
        0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 0.0 0.0 1.0 1.0
    ]
    
    # Elements
    # Hex 1: 1,2,3,4, 5,6,7,8
    # Hex 2: 2,9,10,3, 6,11,12,7 (Careful with ordering)
    # Standard Hex ordering: 1-2-3-4 (bottom), 5-6-7-8 (top)
    # Face 2 (Right) of Hex 1 is 2-3-7-6.
    # For Hex 2, this should be Face 4 (Left) or similar.
    # Hex 2 nodes:
    # Bottom: 2, 9, 10, 3
    # Top: 6, 11, 12, 7
    
    elements = [
        1 2;
        2 9;
        3 10;
        4 3;
        5 6;
        6 11;
        7 12;
        8 7
    ]
    
    tags = [1, 1]
    
    mesh = HexahedraMesh(2, nodes, elements, tags)
    
    basis = RBFBasis(mesh)
    
    # Shared face: 2-3-7-6
    # Hex 1 Face 1: 2-3-7-6
    # Hex 2 Face 4: 2-9-12-6? No.
    # Hex 2 Face 4 (Left): 2-9-11-6? No.
    # Hex 2 Face 4 is (1,2,6,5) local indices.
    # Local 1 is 2. Local 2 is 9. Local 6 is 11. Local 5 is 6.
    # So Face 4 of Hex 2 is 2-9-11-6.
    # Wait, shared face is 2-3-7-6.
    # In Hex 2: 2 is node 1, 3 is node 4, 7 is node 8, 6 is node 5.
    # So shared face is 1-4-8-5 in local indices of Hex 2.
    # That corresponds to Face 2 in my list? No.
    # My list: 2. (1,4,8,5).
    # So Hex 1 Face 1 matches Hex 2 Face 2.
    
    # We expect 11 total basis functions:
    # 2 hexes × 6 faces = 12 faces; 1 shared internal face → 1 internal + 10 boundary = 11
    
    @test num_basis(basis) == 11
    
    # Find the internal (non-boundary) basis function
    internal_idx = findfirst(f -> !f.is_boundary, basis.functions)
    @test internal_idx !== nothing
    
    rbf = basis.functions[internal_idx]
    @test rbf.is_boundary == false
    
    # Area of face (unit square) = 1
    @test isapprox(rbf.area, 1.0)
end
