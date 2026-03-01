#!/usr/bin/env julia
"""
Performance profiling for the 3 underperforming cases.
Measures each component separately to identify exact bottlenecks.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using Printf
using Statistics

const MESH_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles")

function profile_efie_jet()
    println("\n" * "="^60)
    println("Case 1: Jet EFIE Direct (N=14559)")
    println("="^60)
    
    freq = 1e8
    set_frequency!(freq)
    
    # Mesh + basis
    t0 = time()
    mesh = read_nas_mesh(joinpath(MESH_DIR, "jet_100MHz.nas"), scale=1.0)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    t_mesh = time() - t0
    @printf("  Mesh+Basis: %.3f s  (N=%d)\n", t_mesh, N)
    
    # EFIE operator
    efie = EFIE(freq)
    
    # Assembly (3 trials)
    times = Float64[]
    local Z
    for trial in 1:3
        t = @elapsed begin
            Z = assemble_impedance_matrix(efie, basis)
        end
        push!(times, t)
        @printf("  EFIE Assembly trial %d: %.3f s\n", trial, t)
    end
    t_assembly = median(times)
    @printf("  EFIE Assembly (median): %.3f s\n", t_assembly)
    
    # Excitation
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    t0 = time()
    V = excitation_vector(efie, source, basis)
    t_exc = time() - t0
    @printf("  Excitation: %.3f s\n", t_exc)
    
    # LU solve
    t0 = time()
    I_coeff = Z \ V
    t_solve = time() - t0
    @printf("  LU Solve: %.3f s\n", t_solve)
    
    # RCS
    theta = range(-π, π, length=73)
    phi = [0.0]
    t0 = time()
    rcs = radarCrossSection(theta, phi, I_coeff, basis)
    t_rcs = time() - t0
    @printf("  RCS: %.3f s\n", t_rcs)
    
    total = t_mesh + t_assembly + t_exc + t_solve + t_rcs
    @printf("  TOTAL: %.3f s  (Legacy: 46.43s, Ratio: %.2f×)\n", total, total/46.43)
end

function profile_cfie_jet()
    println("\n" * "="^60)
    println("Case 2: Jet CFIE Direct (N=14559)")
    println("="^60)
    
    freq = 1e8
    set_frequency!(freq)
    
    # Mesh + basis
    t0 = time()
    mesh = read_nas_mesh(joinpath(MESH_DIR, "jet_100MHz.nas"), scale=1.0)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    t_mesh = time() - t0
    @printf("  Mesh+Basis: %.3f s  (N=%d)\n", t_mesh, N)
    
    cfie = CFIE(freq, 0.5)
    
    # CFIE Assembly (3 trials)
    times_cfie = Float64[]
    local Z_cfie
    for trial in 1:3
        t = @elapsed begin
            Z_cfie = assemble_impedance_matrix(cfie, basis)
        end
        push!(times_cfie, t)
        @printf("  CFIE Assembly trial %d: %.3f s\n", trial, t)
    end
    t_cfie = median(times_cfie)
    @printf("  CFIE Assembly (median): %.3f s\n", t_cfie)
    
    # Breakdown: measure EFIE and MFIE separately (1 trial each)
    t_efie = @elapsed assemble_impedance_matrix(cfie.efie, basis)
    t_mfie = @elapsed assemble_impedance_matrix(cfie.mfie, basis)
    @printf("    .EFIE sub-assembly: %.3f s\n", t_efie)
    @printf("    .MFIE sub-assembly: %.3f s\n", t_mfie)
    @printf("    .Overhead (matrix add etc): %.3f s\n", t_cfie - t_efie - t_mfie)
    
    # LU solve
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    t0 = time()
    I_coeff = Z_cfie \ V
    t_solve = time() - t0
    @printf("  LU Solve: %.3f s\n", t_solve)
    
    total = t_mesh + t_cfie + t_solve
    @printf("  TOTAL: %.3f s  (Legacy: 64.21s, Ratio: %.2f×)\n", total, total/64.21)
end

function profile_scfie_plate()
    println("\n" * "="^60)
    println("Case 3: SCFIE Direct (Mixed Plate)")
    println("="^60)
    
    freq = 12e8
    set_frequency!(freq)
    
    # Mesh + basis
    t0 = time()
    mesh_file = joinpath(MESH_DIR, "plate_and_metal_1dot2GHz.nas")
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=1.0)
    rwg_basis = RWGBasis(surf_mesh)
    swg_basis = SWGBasis(vol_mesh)
    n_surf = num_basis(rwg_basis)
    n_vol = num_basis(swg_basis)
    N = n_surf + n_vol
    
    n_tet = EMSuite.CoreModule.num_elements(vol_mesh)
    eps_r = ComplexF64(2.0 * (1 - 0.0002im))
    perms = fill(eps_r, n_tet)
    t_mesh = time() - t0
    @printf("  Mesh+Basis: %.3f s  (N_surf=%d, N_vol=%d, N=%d, n_tet=%d)\n", t_mesh, n_surf, n_vol, N, n_tet)
    
    scfie = SCFIE(freq, perms; alpha=0.5)
    
    # Full SCFIE Assembly (3 trials)
    times = Float64[]
    local Z
    for trial in 1:3
        t = @elapsed begin
            Z = assemble_impedance_matrix(scfie, rwg_basis, swg_basis)
        end
        push!(times, t)
        @printf("  SCFIE Assembly trial %d: %.3f s\n", trial, t)
    end
    t_assembly = median(times)
    @printf("  SCFIE Assembly (median): %.3f s\n", t_assembly)
    
    # Breakdown: measure sub-blocks separately (1 trial)
    cfie_sub = CFIE(freq, 0.5)
    t_ss = @elapsed assemble_impedance_matrix(cfie_sub, rwg_basis)
    @printf("    .Z_SS (CFIE %dx%d): %.3f s\n", n_surf, n_surf, t_ss)
    
    vefie = VEFIE(freq, perms)
    t_vv = @elapsed assemble_impedance_matrix(vefie, swg_basis)
    @printf("    .Z_VV (VEFIE %dx%d): %.3f s\n", n_vol, n_vol, t_vv)
    
    coupling_and_fss = t_assembly - t_ss - t_vv
    @printf("    .Coupling + Fss: %.3f s (estimated)\n", max(0.0, coupling_and_fss))
    
    # Excitation + Solve
    source = PlaneWave(freq, π/4, π, [0.0, 0.0, 1.0])
    V = excitation_vector(source, rwg_basis, swg_basis)
    t0 = time()
    I_coeff = Z \ V
    t_solve = time() - t0
    @printf("  LU Solve: %.3f s\n", t_solve)
    
    total = t_mesh + t_assembly + t_solve
    @printf("  TOTAL: %.3f s  (Legacy: 66.52s, Ratio: %.2f×)\n", total, total/66.52)
end

println("EMSuite Performance Profiling v2")
println("Julia: ", VERSION)
println("Threads: ", Threads.nthreads())
println()

# Warmup
println("Warming up...")
set_frequency!(3e8)
mesh_w = read_nas_mesh(joinpath(@__DIR__, "plate_benchmark.nas"), scale=1.0)
basis_w = RWGBasis(mesh_w)
Z_w = assemble_impedance_matrix(EFIE(3e8), basis_w)
println("Warmup complete. N=", num_basis(basis_w))

# Profile each case
profile_efie_jet()
profile_cfie_jet()
profile_scfie_plate()

println("\n" * "="^60)
println("DONE — All profiling complete")
println("="^60)
