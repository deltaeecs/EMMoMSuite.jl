using EMMoMSuite
using EMMoMSuite.Geometry
using Printf

function check_mesh()
    mesh_file = joinpath(@__DIR__, "../../../deps/fixtures/AllinOne/meshfiles/sphere_600MHz.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
    mesh = read_nas_mesh(mesh_file)
    println("Mesh Type: $(typeof(mesh))")
    println("Vertices: $(num_vertices(mesh))")
    println("Elements: $(num_elements(mesh))")
    
    # Calculate Bounding Box
    min_coords = [Inf, Inf, Inf]
    max_coords = [-Inf, -Inf, -Inf]
    
    nodes = mesh.node
    for i in 1:size(nodes, 2)
        v = nodes[:, i]
        min_coords = min.(min_coords, v)
        max_coords = max.(max_coords, v)
    end
    
    println("Bounding Box:")
    println("  Min: $min_coords")
    println("  Max: $max_coords")
    
    dims = max_coords - min_coords
    println("  Dimensions: $dims")
    
    # Check if it looks like a sphere (approx equal dimensions)
    if abs(dims[1] - dims[2]) < 0.1 && abs(dims[2] - dims[3]) < 0.1
        println("Geometry appears to be a Sphere or Cube.")
        radius = dims[1] / 2
        println("Approx Radius: $radius")
    else
        println("Geometry is likely NOT a sphere.")
    end
end

check_mesh()
