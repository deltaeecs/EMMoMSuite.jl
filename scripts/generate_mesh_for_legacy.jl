using EMSuite
using EMSuite.Geometry

# Parameters (Must match verify_SEFIE_mlfma_rcs.jl)
radius = 1.0
n_theta = 24
n_phi = 48

println("Generating sphere mesh (radius=$radius, n_theta=$n_theta, n_phi=$n_phi)...")
mesh = generate_sphere_mesh(radius, n_theta, n_phi)

output_file = joinpath(@__DIR__, "..", "temp_sphere.nas")
println("Saving to $output_file...")
write_nas_mesh(output_file, mesh)
println("Done.")
