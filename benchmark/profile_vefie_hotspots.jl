#!/usr/bin/env julia
"""
profile_vefie_hotspots.jl

快速性能分析：定位 V-EFIE 组装的热点函数。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using Profile
using Printf

const MOM_DIR = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles")

println("="^70)
println("  V-EFIE Performance Hotspot Analysis")
println("="^70)

# ── Load small mesh for quick profiling ─────────────────────────────────────
mesh_file = joinpath(MOM_DIR, "plate_and_metal_1dot2GHz.nas")
isfile(mesh_file) || error("Mesh not found: $mesh_file")

println("\nLoading mesh...")
surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=1.0)
vol_basis = SWGBasis(vol_mesh)
n_vol = num_basis(vol_basis)
println("  SWG basis functions: $n_vol")

freq = 1.2e9
eps_r = 2.0 * (1.0 - 0.0002im)
permittivities = fill(eps_r, num_elements(vol_mesh))

# ── Warm-up run ──────────────────────────────────────────────────────────────
println("\n[1/2] Warm-up compilation...")
vefie = VEFIE(freq, permittivities)
GC.gc()
@time assemble_impedance_matrix(vefie, vol_basis, permittivities)

# ── Profile run ──────────────────────────────────────────────────────────────
println("\n[2/2] Profiling assembly...")
vefie = VEFIE(freq, permittivities)
GC.gc()

Profile.clear()
@profile assemble_impedance_matrix(vefie, vol_basis, permittivities)

println("\n" * "="^70)
println("  Top 20 Hotspot Functions (by self time)")
println("="^70)

# Print flat profile (top 20 by self time)
Profile.print(maxdepth=15, mincount=100, sortedby=:count, format=:flat)

println("\n" * "="^70)
println("  Analysis complete. Key observations:")
println("="^70)
println("  1. Check which functions dominate CPU time")
println("  2. Look for unnecessary allocations (GC time)")
println("  3. Identify if bottleneck is:")
println("     - Green's function computation (fast_green_func)")
println("     - Basis function evaluation")
println("     - Matrix assembly / locking")
println("     - Memory access patterns")
println("="^70)
