using EMSuite
using JLD2
using LinearAlgebra
using Test
using StaticArrays
using Printf
using Statistics

# Load Golden Data
golden_file = joinpath(@__DIR__, "..", "legacy_golden_data.jld2")
if !isfile(golden_file)
    error("Golden data file not found: $golden_file")
end

println("Loading Golden Data...")
data = load(golden_file)
Z_legacy = data["Z"]
tri_data = data["tri_data"]
basis_data = data["bf_data"]

# Extract Triangle 1 Data
tri_idx = 1
tri_info = tri_data[tri_idx]
l1 = abs(tri_info["edge_lengths"][1])
l2 = abs(tri_info["edge_lengths"][2])
l3 = abs(tri_info["edge_lengths"][3])
area = tri_info["area"]
area2 = area^2

println("Triangle 1:")
println("  Edges: ", [l1, l2, l3])
println("  Area: ", area)

# Compute Singular Integrals using EMSuite
using EMSuite.IntegralEquations.EFIEModule.Singularities
sF1 = Singularities.singularF1(l1, l2, l3)
sF21 = Singularities.singularF21(l1, l2, l3, area2)

println("EMSuite Singular Integrals:")
println("  sF1: ", sF1)
println("  sF21: ", sF21)

# Try to load Legacy Singular Integrals if possible (requires MoM_Kernels)
# Since we can't easily load MoM_Kernels here without full env, we rely on the grep results.
# The grep results showed the implementation is IDENTICAL.

# Check Quadrature Weights
# EMSuite
gq = EMSuite.Geometry.GaussQuadratureInfo(:Triangle, 7, Float64)
println("EMSuite Quadrature (N=7):")
println("  Sum Weights: ", sum(gq.weight))
println("  Weights: ", gq.weight)

# Legacy Parity Check
# If Legacy weights sum to 0.5, then EMSuite (sum=1.0) is 2x larger.
# If Legacy weights sum to 1.0, then they match.

# Check Z[1,1]
Z11_leg = Z_legacy[1, 1]
println("Legacy Z[1,1]: ", Z11_leg)

# Calculate Z[1,1] manually using EMSuite logic
k = 2π * 300e6 / 299792458.0
eta = 376.73031346177
factor = im * k * eta / (16 * π)
C4divk2 = 4 / k^2

F1 = C4divk2 * sF1
val_singular = sF21 - F1

println("Manual Calculation:")
println("  k: ", k)
println("  eta: ", eta)
println("  factor: ", factor)
println("  C4divk2: ", C4divk2)
println("  F1: ", F1)
println("  val_singular: ", val_singular)

# Extract Basis Function 1 Data
bf_info = basis_data[1]
edge_len_bf = bf_info["edge_length"]
println("Basis Function 1:")
println("  Edge Length: ", edge_len_bf)

# Check if BF1 is on Triangle 1
tri_ids = bf_info["triangles"]
println("  Triangles: ", tri_ids)

if 1 in tri_ids
    println("  BF1 is on Triangle 1.")
else
    println("  BF1 is NOT on Triangle 1. This script assumes Z[1,1] is self-term of Tri 1.")
end

# Re-calculate Z using correct edge length
lm = edge_len_bf
ln = edge_len_bf
inv_areas = 1.0 / area2
Z_singular = val_singular * area2 * inv_areas * lm * ln * factor
println("  Z_singular (Corrected Edge): ", Z_singular)

ratio = Z_singular / Z11_leg
println("  Ratio (Corrected / Legacy): ", ratio)
println("  Inverse Ratio: ", 1.0/ratio)


