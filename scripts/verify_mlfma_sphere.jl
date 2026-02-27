using EMSuite.CoreModule
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using StaticArrays

function verify_mlfma()
    # 1. Create Mesh (Sphere)
    radius = 1.0
    # mesh = load_mesh("sphere_r1.nas") 
    # If not, generate a simple mesh
    if !isfile("sphere_r1.nas")
        println("Generating sphere mesh...")
        # Simple octahedron subdivision or similar?
        # For now, let's assume we can use a built-in or just create a small plate
        # Actually, let's use a plate for simplicity if sphere is not available
        # But sphere is better for MLFMA testing (closed surface)
        
        # Let's try to load the one from previous tests if available
        # Or just define a small mesh manually
    end
    
    # Better: Use a small mesh defined here
    # 2 triangles
    # vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    # elements = [1 2 3; 1 2 4; 1 3 4; 2 3 4] # Tetrahedron
    # mesh = TriangleMesh(vertices, elements)
    
    # Let's use a slightly larger mesh for MLFMA to have levels
    # Load sphere mesh from previous context if possible
    # Or use `EMSuite.Geometry.read_nastran` if I have a file.
    
    # I'll assume "sphere_r1.nas" exists in current dir or I'll create a dummy one.
    # Actually, I'll just create a simple mesh in code.
    
    # Tetrahedron
    nodes = [
        1.0 0.0 -1.0/sqrt(2);
        -1.0 0.0 -1.0/sqrt(2);
        0.0 1.0 1.0/sqrt(2);
        0.0 -1.0 1.0/sqrt(2)
    ]
    # Scale to radius 1
    nodes = nodes ./ norm(nodes[1,:])
    
    # Faces
    faces = [
        1 3 2;
        1 2 4;
        1 4 3;
        2 3 4
    ]
    
    faces_mat = collect(faces')
    nt = size(faces_mat, 2)
    mesh = TriangleMesh(nt, collect(nodes'), faces_mat)
    
    # 2. Basis
    basis = RWGBasis(mesh)
    num_edges = length(basis.functions)
    println("Number of unknowns: ", num_edges)
    
    # 3. Operator
    freq = 300e6 # 300 MHz -> lambda = 1m
    # k = 2π * freq / 299792458.0
    operator = EFIE(freq)
    
    # 4. Dense Matrix
    println("Assembling dense matrix...")
    Z_dense = assemble_impedance_matrix(operator, basis)
    
    # 5. MLFMA
    println("Assembling MLFMA...")
    # Force at least 2 levels
    # Box size ~ 2.0. Leaf size ~ 0.5 -> 2 levels.
    mlfma = MLFMAOperator(operator, basis, 0.5)
    
    # 6. Compare
    x = ones(ComplexF64, num_edges)
    
    y_dense = Z_dense * x
    y_mlfma = zeros(ComplexF64, num_edges)
    mul!(y_mlfma, mlfma, x)
    
    diff = norm(y_dense - y_mlfma)
    rel_err = diff / norm(y_dense)
    
    println("Dense Norm: ", norm(y_dense))
    println("MLFMA Norm: ", norm(y_mlfma))
    println("Difference: ", diff)
    println("Relative Error: ", rel_err)
    println("Ratio: ", norm(y_mlfma) / norm(y_dense))
    
end

verify_mlfma()
