using EMSuite
using JLD2
using LinearAlgebra
using Printf
using SparseArrays

# Load reference data
data = load("../reference_data.jld2")
nodes = data["nodes"]
triangles = data["triangles"]
ref_Z_EFIE = data["Z_EFIE"]
ref_Z_MFIE = data["Z_MFIE"]
freq = data["freq"]

println("Loaded reference data.")

# Create EMSuite objects
mesh = TriangleMesh(size(triangles, 2), nodes, triangles)
basis = RWGBasis(mesh)
efie = EFIE(freq)
mfie = MFIE(freq)

# --- Verify EFIE ---
println("\n--- Verifying EFIE Matrix ---")
Z_EFIE = assemble_impedance_matrix(efie, basis)

# Compare
diff_EFIE = Z_EFIE - ref_Z_EFIE
norm_diff_EFIE = norm(diff_EFIE) / norm(ref_Z_EFIE)
max_diff_EFIE = maximum(abs.(diff_EFIE))

println("EFIE Relative Norm Error: ", norm_diff_EFIE)
println("EFIE Max Element Error: ", max_diff_EFIE)

if norm_diff_EFIE < 1e-10
    println("EFIE VERIFIED: SUCCESS")
else
    println("EFIE VERIFIED: FAILURE")
    # Print some details
    rows, cols = size(Z_EFIE)
    for i in 1:min(5, rows)
        for j in 1:min(5, cols)
            println("Z[$i,$j]: EMSuite=$(Z_EFIE[i,j]), Ref=$(ref_Z_EFIE[i,j])")
        end
    end
end

# --- Verify MFIE ---
println("\n--- Verifying MFIE Matrix ---")
Z_MFIE = assemble_impedance_matrix(mfie, basis)

# Compare
diff_MFIE = Z_MFIE - ref_Z_MFIE
norm_diff_MFIE = norm(diff_MFIE) / norm(ref_Z_MFIE)
max_diff_MFIE = maximum(abs.(diff_MFIE))

println("MFIE Relative Norm Error: ", norm_diff_MFIE)
println("MFIE Max Element Error: ", max_diff_MFIE)

if norm_diff_MFIE < 1e-10
    println("MFIE VERIFIED: SUCCESS")
else
    println("MFIE VERIFIED: FAILURE")
    # Print some details
    rows, cols = size(Z_MFIE)
    for i in 1:min(5, rows)
        for j in 1:min(5, cols)
            println("Z[$i,$j]: EMSuite=$(Z_MFIE[i,j]), Ref=$(ref_Z_MFIE[i,j])")
        end
    end
end
