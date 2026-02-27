using JLD2
using StaticArrays

data = load("c:/Users/12253/OneDrive/MoM/EMSuite/sphere_mesh_data.jld2")
println("Keys: ", keys(data))
mesh = data["mesh"]
println("Mesh Type: ", typeof(mesh))
println("Num Vertices: ", length(mesh.vertices))
if hasproperty(mesh, :tetrahedra)
    println("Num Tetrahedra: ", length(mesh.tetrahedra))
end
if hasproperty(mesh, :triangles)
    println("Num Triangles: ", length(mesh.triangles))
end
