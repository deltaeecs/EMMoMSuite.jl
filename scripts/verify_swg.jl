using EMSuite.CoreModule
using EMSuite.Geometry
using EMSuite.BasisFunctions
using LinearAlgebra
using StaticArrays

function verify_swg()
    println("Verifying SWG Basis...")
    
    # Create a simple mesh: Two tetrahedra sharing a face
    # Nodes
    # 1: (0,0,0)
    # 2: (1,0,0)
    # 3: (0,1,0)
    # 4: (0,0,1)
    # 5: (0,0,-1)
    
    nodes = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0;
        0.0 0.0 -1.0
    ]
    
    # Tet 1: 1, 2, 3, 4
    # Tet 2: 1, 2, 3, 5
    # Shared face: 1, 2, 3
    
    tets = [
        1 2 3 4;
        1 2 3 5
    ]
    
    mesh = TetrahedraMesh(2, collect(nodes'), collect(tets'))
    
    basis = SWGBasis(mesh)
    
    println("Number of basis functions: ", num_basis(basis))
    
    # Expected:
    # Tet 1 has 4 faces.
    # Tet 2 has 4 faces.
    # Shared face (1,2,3) -> 1 internal basis function.
    # Boundary faces:
    # Tet 1: (1,2,4), (2,3,4), (3,1,4) -> 3 boundary
    # Tet 2: (1,2,5), (2,3,5), (3,1,5) -> 3 boundary
    # Total faces = 7.
    # Internal = 1.
    # If boundary are skipped, num_basis = 1.
    
    for (i, swg) in enumerate(basis.functions)
        println("Basis $i: ID=$(swg.id), Boundary=$(swg.is_boundary), Support=$(swg.support)")
    end
end

verify_swg()
