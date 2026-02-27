using EMSuite
using JLD2

# Generate a small sphere mesh
radius = 0.5
n_theta = 6
n_phi = 12
mesh = generate_sphere_mesh(radius, n_theta, n_phi)

# Extract raw data
nodes = mesh.node
triangles = mesh.triangles

# Save to file
@save "sphere_mesh_data.jld2" nodes triangles radius
println("Sphere mesh saved to sphere_mesh_data.jld2")
