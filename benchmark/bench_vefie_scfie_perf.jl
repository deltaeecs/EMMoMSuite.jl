"""
bench_vefie_scfie_perf.jl

Quick assembly-time benchmark for VEFIE (SWG) and SCFIE.
Compares against Phase 8.9 baseline (VEFIE=66.24s, SCFIE=96.94s).
"""
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using Printf

const MOM_DIR = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles")

function bench_vefie()
    mesh_file = joinpath(MOM_DIR, "Tetra.nas")
    if !isfile(mesh_file)
        @warn "Tetra.nas not found at $mesh_file, skipping VEFIE bench."
        return
    end
    println("\n" * "="^60)
    println("  VEFIE (SWG) Assembly Benchmark")
    println("="^60)
    println("Threads: $(Threads.nthreads())")

    mesh = read_nas_mesh(mesh_file)
    max_coord = maximum(abs.(mesh.node))
    if max_coord > 10.0
        mesh.node .*= 0.001
    end
    center = (minimum(mesh.node, dims=2) + maximum(mesh.node, dims=2)) / 2
    mesh.node .-= center

    freq = 300e6
    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))

    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("N (SWG unknowns) = $N")

    vefie = VEFIE(freq, permittivities)

    # Warm-up (small mesh would be ideal but just measure directly)
    println("Assembling Z (VEFIE)...")
    t = @elapsed Z = assemble_impedance_matrix(vefie, basis, permittivities)
    @printf "  Assembly time: %.2f s\n" t
    @printf "  Phase 8.9 baseline: 66.24 s (Legacy: 46.13 s)\n"
    @printf "  Speedup vs Phase 8.9 baseline: %.2f×\n" (66.24 / t)
    @printf "  vs Legacy: %.2f×\n" (46.13 / t)
    return t
end

function bench_scfie()
    surf_file = joinpath(MOM_DIR, "plate_and_metal_1dot2GHz.nas")
    if !isfile(surf_file)
        @warn "plate_and_metal_1dot2GHz.nas not found, skipping SCFIE bench."
        return
    end
    println("\n" * "="^60)
    println("  SCFIE (RWG+SWG) Assembly Benchmark")
    println("="^60)
    println("Threads: $(Threads.nthreads())")

    freq = 1.2e9
    eps_r = 2.0 * (1.0 - 0.0002im)

    mesh = read_nas_mesh(surf_file)

    # This mesh has both surface tris and volume tets
    surf_mesh = extract_surface_mesh(mesh)
    vol_mesh  = extract_volume_mesh(mesh)

    if isnothing(surf_mesh) || isnothing(vol_mesh)
        @warn "Could not split surf/vol mesh, skipping."
        return
    end

    permittivities = fill(eps_r, num_elements(vol_mesh))
    surf_basis = RWGBasis(surf_mesh)
    vol_basis  = SWGBasis(vol_mesh)
    N_total    = num_basis(surf_basis) + num_basis(vol_basis)
    @printf "N_surf=%d  N_vol=%d  N_total=%d\n" num_basis(surf_basis) num_basis(vol_basis) N_total

    scfie = SCFIE(freq, permittivities)

    println("Assembling Z (SCFIE)...")
    t = @elapsed Z = assemble_impedance_matrix(scfie, surf_basis, vol_basis)
    @printf "  Assembly time: %.2f s\n" t
    @printf "  Phase 8.9 baseline: 96.94 s (Legacy: 66.68 s)\n"
    @printf "  Speedup vs Phase 8.9 baseline: %.2f×\n" (96.94 / t)
    @printf "  vs Legacy: %.2f×\n" (66.68 / t)
    return t
end

println("EMSuite VEFIE/SCFIE Performance Benchmark")
println("Phase 9 optimizations: row-parallel cyclic + O(4) direct scatter")

t_vefie = bench_vefie()
t_scfie = bench_scfie()

println("\nSummary:")
@printf "  VEFIE: %.2f s\n" (t_vefie !== nothing ? t_vefie : NaN)
@printf "  SCFIE: %.2f s\n" (t_scfie !== nothing ? t_scfie : NaN)
