using EMSuite.Geometry

mesh_file = joinpath("c:/Users/12253/OneDrive/MoM/MoM_AllinOne/meshfiles/plate_and_metal_1dot2GHz.nas")
if isfile(mesh_file)
    println("File exists: ", mesh_file)
    try
        mesh = read_nas_mesh(mesh_file)
        println("Mesh read successfully.")
        println("Type: ", typeof(mesh))
        println("Unique Tags: ", unique(mesh.tags))
        if mesh isa TetrahedraMesh
            println("Tetrahedra: ", length(mesh.tetras))
        elseif mesh isa TriangleMesh
            println("Triangles: ", length(mesh.triangles))
        end
    catch e
        println("Error reading mesh: ", e)
    end
else
    println("File not found: ", mesh_file)
end
