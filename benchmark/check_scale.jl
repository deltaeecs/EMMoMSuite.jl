using EMSuite
using EMSuite.Geometry

mesh_file = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
println("Reading with scale=1.0:")
s1, v1 = read_mixed_nas_mesh(mesh_file; scale=1.0)
println("  vol x range: ", extrema(v1.node[1,:]))
println("  vol y range: ", extrema(v1.node[2,:]))

println("Reading with scale=0.001:")
s2, v2 = read_mixed_nas_mesh(mesh_file; scale=0.001)
println("  vol x range: ", extrema(v2.node[1,:]))
println("  vol y range: ", extrema(v2.node[2,:]))
