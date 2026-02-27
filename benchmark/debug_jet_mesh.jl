using EMSuite
using LinearAlgebra
using Printf

function debug_jet_mesh()
    println("=== Debugging Jet Mesh and Basis ===")
    
    # 1. Load Mesh
    mesh_file = joinpath("f:/OneDrive/MoM/MoM_AllinOne/meshfiles/jet_100MHz.nas")
    println("Loading mesh: $mesh_file")
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    println("Number of elements: ", num_elements(mesh))
    println("Number of vertices: ", num_vertices(mesh))
    
    # Check for degenerate triangles
    min_area = Inf
    max_area = -Inf
    verts = vertices(mesh)
    elems = elements(mesh)
    
    for i in 1:num_elements(mesh)
        # Assuming 3-node triangles
        idx = elems[:, i]
        v1 = verts[:, idx[1]]
        v2 = verts[:, idx[2]]
        v3 = verts[:, idx[3]]
        
        area = 0.5 * norm(cross(v2-v1, v3-v1))
        min_area = min(min_area, area)
        max_area = max(max_area, area)
        if area < 1e-10
            println("Warning: Degenerate triangle $i, Area: $area")
        end
    end
    println("Min Area: $min_area")
    println("Max Area: $max_area")
    
    # 2. Basis
    println("Setting up RWG basis...")
    basis = RWGBasis(mesh)
    println("Number of unknowns (edges): ", num_basis(basis))
    
    # Check for non-manifold edges (shared by > 2 triangles)
    # We need to reconstruct the edge list logic to check this
    nt = num_elements(mesh)
    tris = elements(mesh)
    EdgeInfo = Tuple{Int, Int}
    all_edges = Vector{EdgeInfo}(undef, nt * 3)
    idx = 1
    for t in 1:nt
        v1, v2, v3 = tris[:, t]
        all_edges[idx] = (min(v2, v3), max(v2, v3))
        idx += 1
        all_edges[idx] = (min(v3, v1), max(v3, v1))
        idx += 1
        all_edges[idx] = (min(v1, v2), max(v1, v2))
        idx += 1
    end
    sort!(all_edges)
    
    non_manifold_count = 0
    boundary_count = 0
    internal_count = 0
    
    i = 1
    while i <= length(all_edges)
        count = 1
        while i + count <= length(all_edges) && all_edges[i+count] == all_edges[i]
            count += 1
        end
        
        if count == 1
            boundary_count += 1
        elseif count == 2
            internal_count += 1
        else
            non_manifold_count += 1
            println("Warning: Edge $(all_edges[i]) shared by $count triangles!")
        end
        i += count
    end
    
    println("Internal Edges (Manifold): $internal_count")
    println("Boundary Edges: $boundary_count")
    println("Non-Manifold Edges: $non_manifold_count")
    
    if num_basis(basis) != internal_count
        println("WARNING: Basis count mismatch! Basis: $(num_basis(basis)), Internal Edges: $internal_count")
    else
        println("Basis count matches internal edges.")
    end
    
end

debug_jet_mesh()
