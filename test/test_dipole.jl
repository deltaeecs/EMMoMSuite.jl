using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.CoreModule
using EMSuite.PostProcessing
using LinearAlgebra
using StaticArrays
using Test

# 1. Generate Mesh (Dipole Antenna)
# Dipole length L = lambda/2
freq = 300e6
lambda = 3e8 / freq
L = lambda / 2
radius = L / 200

# Cylinder mesh
# We need a thin cylinder.
# Let's use the new generator
mesh = generate_cylinder_mesh(radius, L, 12, 20; closed=true)

println("Mesh generated: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) elements")

# 2. Setup Basis
basis = RWGBasis(mesh)
println("Basis functions: $(num_basis(basis))")

# 3. Find the center edge for excitation
# We look for edges at z=0

    
# Find all edges at the center (z=0)
center_edges = Int[]
tol = L / 100 # Tolerance for z=0

for n in 1:num_basis(basis)
    bf = basis.functions[n]
    # Get edge center
    # We can use bf.center if available, or calculate it
    # bf struct has center field? Let's check RWG.jl
    # Yes, it has center::SVector{3, FT}
    
    if abs(bf.center[3]) < tol
        push!(center_edges, n)
    end
end

println("Excitation edges: $center_edges")

# 4. Define Source (Delta Gap)
source = DeltaGapSource(freq, center_edges, 1.0 + 0im)

# 5. Assemble and Solve
efie = EFIE(freq)
Z = assemble_impedance_matrix(efie, basis)
V = excitation_vector(source, basis)

solver = GMRESSolver(tol=1e-4, maxiter=500)
I_sol = solve!(solver, Z, V)

# 6. Calculate Input Impedance
# Zin = V_gap / I_total
# V_gap = 1.0 (since we set V[k] = 1.0 * l_k)
# I_total = sum(I_k * l_k) for k in center_edges

I_total = zero(ComplexF64)
for idx in center_edges
    l_k = basis.functions[idx].edge_length
    global I_total += I_sol[idx] * l_k
end

Zin = 1.0 / I_total
println("Input Impedance: $Zin Ohms")

# Theoretical approx for half-wave dipole: 73 + j42.5
# Note: Thin wire approximation vs surface mesh might differ.

# 7. Far Field Pattern
theta = collect(range(0, pi, length=37))
phi = [0.0]

# Need to convert mesh to TriangleInfo for farField function
# Or update farField to accept mesh
# The current farField signature is:
# farField(θs_obs, ϕs_obs, ICoeff, trianglesInfo, source, BFT)
# It seems it expects trianglesInfo.

# Let's create trianglesInfo
# Use the updated get_triangle_info from BasisFunctions which now handles signed IDs
trianglesInfo = [get_triangle_info(mesh, basis, i) for i in 1:num_elements(mesh)]

# Set frequency in parameters for k0
set_frequency!(freq)

ff = farField(theta, phi, I_sol, trianglesInfo, nothing, RWG)

# Check if pattern is omnidirectional in H-plane (xy plane, theta=pi/2)
# and figure-8 in E-plane (xz plane, phi=0)

# E-plane (phi=0, varying theta)
# Max at theta=pi/2, Null at theta=0, pi
max_ff = maximum(norm.(ff))
null_ff_1 = norm(ff[1]) # theta=0
null_ff_2 = norm(ff[end]) # theta=pi

println("Max Far Field: $max_ff")
println("Null 1: $null_ff_1")
println("Null 2: $null_ff_2")

@test null_ff_1 < max_ff * 0.1
@test null_ff_2 < max_ff * 0.1

println("Dipole test completed.")

# 8. Near Field Calculation
println("Calculating Near Field...")
points = [SVector(100.0, 0.0, 0.0), SVector(0.0, 0.0, 100.0)] # Broadside, Axis (Far Field region)
E_near = calculate_near_field(points, basis, I_sol)

println("E-field at (100,0,0) (Broadside): $(norm(E_near[1]))")
println("E-field at (0,0,100) (Axis): $(norm(E_near[2]))")

# Broadside should be much larger than Axis in Far Field
@test norm(E_near[1]) > norm(E_near[2]) * 10
println("Near Field test completed.")
