using EMSuite.Geometry

mesh_file = joinpath("c:/Users/12253/OneDrive/MoM/MoM_AllinOne/meshfiles/sphere_600MHz.nas")
if isfile(mesh_file)
    println("File exists: ", mesh_file)
    # Try to read it
    try
        mesh = read_nas_mesh(mesh_file)
        println("Mesh read successfully.")
        println("Type: ", typeof(mesh))
        if mesh isa TetrahedraMesh
            println("Tetrahedra: ", length(mesh.tetrahedra))
        elseif mesh isa TriangleMesh
            println("Triangles: ", length(mesh.triangles))
        end
    catch e
        println("Error reading mesh: ", e)
    end
else
    println("File not found: ", mesh_file)
end
