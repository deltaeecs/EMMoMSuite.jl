#!/usr/bin/env julia
"""
Phase 8.5b: Type stability audit using @code_warntype

Checks key hot-path functions for type instabilities (Any types, unions).
"""

using EMSuite
using StaticArrays
using LinearAlgebra
using InteractiveUtils

println("="^60)
println("  Phase 8.5b: Type Stability Audit")
println("="^60)

# Setup: create a small mesh and operators
mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
freq = 3e8

mesh = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis = RWGBasis(mesh)

efie = EFIE(freq)
mfie = MFIE(freq)
cfie = CFIE(freq, 0.5)

# Get TriangleInfo for testing
nt = EMSuite.CoreModule.num_elements(mesh)
tri1 = EMSuite.BasisFunctions.get_triangle_info(mesh, basis, 1)
tri2 = EMSuite.BasisFunctions.get_triangle_info(mesh, basis, 2)
# Find a non-adjacent triangle
tri_far = EMSuite.BasisFunctions.get_triangle_info(mesh, basis, min(100, nt))

N = num_basis(basis)

# Test 1: EFIE interaction (far field)
println("\n" * "-"^60)
println("  1. efie_interaction! (far-field, pre-computed quad points)")
println("-"^60)

gq = efie.gq_far
N_points = length(gq.weight)
quad_points = Vector{SVector{N_points, SVector{3, Float64}}}(undef, nt)
for t in 1:nt
    v_indices = mesh.triangles[:, t]
    v1 = SVector{3, Float64}(mesh.node[:, v_indices[1]])
    v2 = SVector{3, Float64}(mesh.node[:, v_indices[2]])
    v3 = SVector{3, Float64}(mesh.node[:, v_indices[3]])
    quad_points[t] = SVector{N_points, SVector{3, Float64}}(
        v1 * gq.coordinate[1, i] + v2 * gq.coordinate[2, i] + v3 * gq.coordinate[3, i]
        for i in 1:N_points
    )
end

Z = zeros(ComplexF64, 3, 3)
try
    @code_warntype EMSuite.IntegralEquations.EFIEModule.efie_interaction!(Z, efie, tri1, tri_far, quad_points)
catch e
    println("  SKIPPED: $e")
end

# Test 2: EFIE interaction (self-term)
println("\n" * "-"^60)
println("  2. efie_interaction! (self-term)")
println("-"^60)
fill!(Z, 0)
try
    @code_warntype EMSuite.IntegralEquations.EFIEModule.efie_interaction!(Z, efie, tri1, tri1)
catch e
    println("  SKIPPED: $e")
end

# Test 3: MFIE interaction (K-term fast)
println("\n" * "-"^60)
println("  3. mfie_interaction! (K-term)")
println("-"^60)
fill!(Z, 0)
try
    @code_warntype EMSuite.IntegralEquations.MFIEModule.mfie_interaction!(Z, mfie, tri1, tri_far)
catch e
    println("  SKIPPED: $e")
end

# Test 4: MFIE calc_self_term!
println("\n" * "-"^60)
println("  4. MFIE calc_self_term!")
println("-"^60)
fill!(Z, 0)
try
    @code_warntype EMSuite.IntegralEquations.MFIEModule.calc_self_term!(Z, mfie, tri1)
catch e
    println("  SKIPPED: $e")
end

# Test 5: Green's function
println("\n" * "-"^60)
println("  5. green_function_free_space")
println("-"^60)
r1 = SVector(0.0, 0.0, 0.0)
r2 = SVector(1.0, 1.0, 1.0)
try
    @code_warntype EMSuite.IntegralEquations.Kernels.green_function_free_space(r1, r2, 1.0)
catch e
    println("  SKIPPED: $e")
end

# Test 6: assemble_generic
println("\n" * "-"^60)
println("  6. assemble_generic (via assemble_impedance_matrix EFIE)")
println("-"^60)
try
    @code_warntype EMSuite.IntegralEquations.Impedance.assemble_generic(efie, basis, (Z, op, t1, t2) -> EMSuite.IntegralEquations.EFIEModule.efie_interaction!(Z, op, t1, t2); symmetric=true)
catch e
    println("  SKIPPED: $e")
end

# Summary: check for type instabilities by looking at output
println("\n" * "="^60)
println("  Type Stability Audit Complete")
println("  (Check for 'Any' or 'Union' type annotations above)")
println("="^60)
