using EMSuite
using JLD2
using LinearAlgebra
using Test
using StaticArrays

function verify_efie()
    println("Loading reference data...")
    data = load("reference_data.jld2")
    Z_ref = data["Z_EFIE"]
    freq = data["freq"]
    
    println("Frequency: ", freq)
    
    # Setup EMSuite
    println("Reading mesh from JLD2...")
    nodes = data["nodes"]
    triangles = data["triangles"]
    
    # Convert to TriangleMesh
    # nodes might be 3xN or Nx3.
    # triangles might be 3xM or Mx3.
    
    # Check dimensions
    println("Nodes size: ", size(nodes))
    println("Triangles size: ", size(triangles))
    
    # Assuming standard Julia column-major (3xN)
    if size(nodes, 1) != 3
        nodes = nodes'
    end
    if size(triangles, 1) != 3
        triangles = triangles'
    end
    
    # Create mesh manually
    # TriangleMesh(trinum, node, triangles)
    
    n_elems = size(triangles, 2)
    
    # Ensure types
    nodes = convert(Matrix{Float64}, nodes)
    triangles = convert(Matrix{Int}, triangles)
    
    mesh = EMSuite.Geometry.TriangleMesh(n_elems, nodes, triangles)
    
    println("Setting up basis...")
    basis = RWGBasis(mesh)
    
    println("Setting up EFIE...")
    efie = EFIE(freq)
    
    println("Assembling matrix...")
    Z_calc = assemble_impedance_matrix(efie, basis)
    
    println("Comparing matrices...")
    diff = norm(Z_calc - Z_ref) / norm(Z_ref)
    println("Relative Error: ", diff * 100, "%")
    
    # Check self-term specifically
    println("Z_ref[1,1]: ", Z_ref[1,1])
    println("Z_calc[1,1]: ", Z_calc[1,1])
end

verify_efie()
