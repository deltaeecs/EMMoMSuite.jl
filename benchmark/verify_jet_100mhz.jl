using EMSuite
using EMSuite.Core
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.IntegralEquations.Impedance
using EMSuite.Solvers
using EMSuite.PostProcessing
using LinearAlgebra
using CSV
using DataFrames
using Printf
using StaticArrays

# --- Configuration ---
freq = 100e6
c0 = 299792458.0
k = 2 * π * freq / c0

# Paths
mesh_file = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/jet_100MHz.nas")
feko_file = joinpath(@__DIR__, "../../MoM_AllinOne/deps/compare_feko/jet_100MHzRCS.csv")

if !isfile(mesh_file)
    error("Mesh file not found: $mesh_file")
end

# --- 1. Load Mesh ---
println("Loading mesh from: $mesh_file")
nodes, triangles, tags = read_nas_mesh(mesh_file)
trinum = size(triangles, 2)
# Use default constructor for TriangleMesh
mesh = TriangleMesh(trinum, nodes, triangles, tags)
println("Mesh loaded: $(size(nodes, 2)) nodes, $trinum triangles.")

# --- 2. Basis Functions ---
println("Initializing RWG Basis...")
basis = RWGBasis(mesh)
println("Number of edges (unknowns): $(num_basis(basis))")

# --- 3. Integral Equation (EFIE) ---
println("Initializing EFIE Operator...")
efie = EFIE(freq)

println("Assembling Impedance Matrix...")
Z = assemble_impedance_matrix(efie, basis)
println("Matrix size: $(size(Z))")

# --- 4. Excitation ---
println("Setting up Excitation...")
# Incident from (90, 0) -> Propagating towards (90, 180)
# Polarization along z (theta-hat at 90,0 is -z, so alpha=0 means -z? 
# MoM_Basics: E = -theta_hat. At (90,0), theta_hat = (0,0,-1). So E = (0,0,1).
inc_theta = π/2
inc_phi = π # Propagation direction
pol = [0.0, 0.0, 1.0]
source = PlaneWave(freq, inc_theta, inc_phi, pol)

V = excitation_vector(efie, source, basis)

# --- 5. Solve ---
println("Solving System (Direct)...")
I_coeff = Z \ V
println("Solved.")

# --- 6. RCS Calculation ---
println("Calculating RCS...")
theta_obs = collect(range(-π, π, length=721))
phi_obs = [0.0, π/2]

# Precompute triangle info for RCS
tris_info = [get_triangle_info(mesh, basis, t) for t in 1:trinum]

# Calculate RCS for phi=0 cut
rcs_phi0 = radarCrossSection(theta_obs, [0.0], I_coeff, tris_info, RWG)
# Calculate RCS for phi=90 cut
rcs_phi90 = radarCrossSection(theta_obs, [π/2], I_coeff, tris_info, RWG)

# Convert to dBsm
rcs_phi0_db = 10 .* log10.(vec(rcs_phi0))
rcs_phi90_db = 10 .* log10.(vec(rcs_phi90))

# --- 7. Compare with Feko ---
println("Loading Feko Reference...")
if isfile(feko_file)
    df_feko = CSV.read(feko_file, DataFrame; delim=' ', ignorerepeated=true)
    # Feko file columns: THETA PHI magn. phase ...
    # We need to match theta. Feko theta is in degrees.
    # Our theta_obs is in radians.
    theta_deg = rad2deg.(theta_obs)
    
    # Filter Feko data for phi=0
    # Assuming Feko file structure: THETA PHI ...
    # We might need to inspect the file structure more closely if this fails.
    
    # For now, just save our results
    output_file = joinpath(@__DIR__, "jet_100MHz_results.csv")
    println("Saving results to $output_file")
    
    df_out = DataFrame(
        Theta_Deg = theta_deg,
        RCS_Phi0_dB = rcs_phi0_db,
        RCS_Phi90_dB = rcs_phi90_db
    )
    CSV.write(output_file, df_out)
    println("Results saved.")
else
    println("Feko reference file not found: $feko_file")
end
