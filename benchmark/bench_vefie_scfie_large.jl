# bench_vefie_scfie_large.jl
# Quick timing for VEFIE + SCFIE on the Phase 8 benchmark mesh (plate_and_metal_1dot2GHz.nas).
# Expected baselines (4 threads, Phase 8.9):
#   VEFIE: 66.24 s  (Legacy: 46.13 s)
#   SCFIE: 96.94 s  (Legacy: 66.68 s)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using Printf

const MOM_DIR = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles")

println("Threads: $(Threads.nthreads())")
println()

# ── Load mesh ────────────────────────────────────────────────────────────────
mesh_file = joinpath(MOM_DIR, "plate_and_metal_1dot2GHz.nas")
isfile(mesh_file) || error("Mesh not found: $mesh_file")

println("Loading plate_and_metal_1dot2GHz.nas ...")
surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=1.0)
println("  Surface: $(num_elements(surf_mesh)) triangles")
println("  Volume:  $(num_elements(vol_mesh)) tetrahedra")

surf_basis = RWGBasis(surf_mesh)
vol_basis  = SWGBasis(vol_mesh)
n_surf = num_basis(surf_basis)
n_vol  = num_basis(vol_basis)
println("  RWG: $n_surf,  SWG: $n_vol,  Total: $(n_surf+n_vol)")

freq  = 1.2e9
eps_r = 2.0 * (1.0 - 0.0002im)
permittivities = fill(eps_r, num_elements(vol_mesh))

# ── VEFIE (SWG only) ─────────────────────────────────────────────────────────
println("\n--- VEFIE (SWG, N=$n_vol) ---")
vefie = VEFIE(freq, permittivities)
GC.gc()
t_vefie = @elapsed begin
    Z_vv = assemble_impedance_matrix(vefie, vol_basis, permittivities)
end
@printf "  Assembly time: %.2f s\n" t_vefie
@printf "  Phase 8.9 baseline: 66.24 s  →  speedup %.2f×\n" (66.24 / t_vefie)
@printf "  Legacy: 46.13 s  →  ratio %.3f×\n"               (t_vefie / 46.13)

# ── SCFIE (RWG+SWG) ──────────────────────────────────────────────────────────
println("\n--- SCFIE (RWG+SWG, N=$(n_surf+n_vol)) ---")
scfie = SCFIE(freq, permittivities; alpha=0.5)
GC.gc()
t_scfie = @elapsed begin
    Z = assemble_impedance_matrix(scfie, surf_basis, vol_basis)
end
@printf "  Assembly time: %.2f s\n" t_scfie
@printf "  Phase 8.9 baseline: 96.94 s  →  speedup %.2f×\n" (96.94 / t_scfie)
@printf "  Legacy: 66.68 s  →  ratio %.3f×\n"               (t_scfie / 66.68)

println("\nDone.")
