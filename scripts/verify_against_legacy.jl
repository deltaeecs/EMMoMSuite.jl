using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) # Activate EMSuite
using EMSuite
using Serialization
using LinearAlgebra
using Printf
using SparseArrays

# Load Legacy Data
data_file = joinpath(@__DIR__, "legacy_data.jls")
if !isfile(data_file)
    error("Legacy data file not found. Run generate_legacy_data.jl first.")
end

println("Loading legacy data...")
data = deserialize(data_file)
Z_EFIE_old = data["Z_EFIE"]
V_EFIE_old = data["V_EFIE"]
Z_MFIE_old = data["Z_MFIE"]
V_MFIE_old = data["V_MFIE"]

# Setup EMSuite
freq = 300e6
mesh_file = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/sphere_600MHz.nas")
println("Loading mesh from $mesh_file...")
mesh = read_nas_mesh(mesh_file)
println("Mesh loaded. Nodes: $(size(mesh.node, 2)), Triangles: $(size(mesh.triangles, 2))")
basis = RWGBasis(mesh)
println("Basis functions generated. Count: $(length(basis.functions))")

# PlaneWave(freq, theta, phi, polarization_vector)
# Legacy: PlaneWave(π/2, 0.0, 0.0, 1.0) -> theta=90, phi=0, pol_theta=0, pol_phi=1 (y-polarized)
# At theta=90, phi=0: r_hat=(1,0,0), theta_hat=(0,0,-1), phi_hat=(0,1,0)
# So pol_phi=1 means E is along y-axis.
source = PlaneWave(freq, π/2, 0.0, [0.0, 1.0, 0.0])

# --- EFIE ---
println("\n--- EFIE Verification ---")
efie = EFIE(freq)
println("Assembling EFIE Matrix...")
Z_EFIE_new = assemble_impedance_matrix(efie, basis)
println("Calculating EFIE Excitation...")
V_EFIE_new = excitation_vector(efie, source, basis)

# Convert to dense if sparse for comparison
Z_EFIE_new_dense = Matrix(Z_EFIE_new)
Z_EFIE_old_dense = Matrix(Z_EFIE_old)

diff_Z_EFIE = norm(Z_EFIE_new_dense - Z_EFIE_old_dense) / norm(Z_EFIE_old_dense)
diff_V_EFIE = norm(V_EFIE_new - V_EFIE_old) / norm(V_EFIE_old)
@printf "Z Matrix Rel Error: %.4e\n" diff_Z_EFIE
@printf "V Vector Rel Error: %.4e\n" diff_V_EFIE

# --- MFIE ---
println("\n--- MFIE Verification ---")
mfie = MFIE(freq)
println("Assembling MFIE Matrix...")
Z_MFIE_new = assemble_impedance_matrix(mfie, basis)
println("Calculating MFIE Excitation...")
V_MFIE_new = excitation_vector(mfie, source, basis)

Z_MFIE_new_dense = Matrix(Z_MFIE_new)
Z_MFIE_old_dense = Matrix(Z_MFIE_old)

diff_Z_MFIE = norm(Z_MFIE_new_dense - Z_MFIE_old_dense) / norm(Z_MFIE_old_dense)
diff_V_MFIE = norm(V_MFIE_new - V_MFIE_old) / norm(V_MFIE_old)
@printf "Z Matrix Rel Error: %.4e\n" diff_Z_MFIE
@printf "V Vector Rel Error: %.4e\n" diff_V_MFIE

if diff_Z_MFIE > 1e-4
    println("\nDetailed MFIE Analysis:")
    diag_old = diag(Z_MFIE_old_dense)
    diag_new = diag(Z_MFIE_new_dense)
    diff_diag = norm(diag_new - diag_old)/norm(diag_old)
    @printf "Diagonal Rel Error: %.4e\n" diff_diag
    
    # Check a few random elements
    println("Sample Elements (Old vs New):")
    for i in 1:5
        r, c = rand(1:size(Z_MFIE_old, 1), 2)
        @printf "(%d, %d): %s vs %s\n" r c string(Z_MFIE_old_dense[r,c]) string(Z_MFIE_new_dense[r,c])
    end
end
