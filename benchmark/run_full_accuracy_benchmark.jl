# Comprehensive Accuracy Benchmark for EMMoMSuite (v2 - memory-aware)
# Runs:
#   1. SEFIE Direct (Jet, EFIE, 100MHz) vs Legacy
#   2. SCFIE Direct (Jet, CFIE alpha=0.5, 100MHz) vs SEFIE Direct
#   3. SEFIE MLFMA (Jet, EFIE+GMRES, 100MHz) vs Legacy & vs Direct
#   4. SCFIE MLFMA (Sphere, CFIE alpha=0.5, 600MHz) vs Legacy

using EMMoMSuite
using LinearAlgebra
using Printf
using Statistics
using CSV, DataFrames
using SparseArrays

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne")
const LEGACY_BASELINE_DIR = joinpath(@__DIR__, "..", "test_results", "legacy_baseline")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "test_results", "emsuite_verification")

mkpath(OUTPUT_DIR)

# LU Preconditioner wrapper for GMRES
struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

# Utility: compare RCS against Legacy baseline CSV
function compare_rcs(RCS_dBsm, θs, legacy_file; label="")
    if !isfile(legacy_file)
        println("  [WARN] Legacy baseline not found: $legacy_file")
        return nothing
    end
    df = CSV.read(legacy_file, DataFrame)
    
    d0  = RCS_dBsm[:, 1] .- df.RCS_Phi0_dBsm
    d90 = RCS_dBsm[:, 2] .- df.RCS_Phi90_dBsm
    
    println("\n  --- $label Phi=0 vs Legacy ---")
    @printf("  Mean Diff : %+.4f dB\n", mean(d0))
    @printf("  RMSE      : %.4f dB\n",  sqrt(mean(d0.^2)))
    @printf("  Max |Diff|: %.4f dB\n",  maximum(abs.(d0)))
    
    println("  --- $label Phi=90 vs Legacy ---")
    @printf("  Mean Diff : %+.4f dB\n", mean(d90))
    @printf("  RMSE      : %.4f dB\n",  sqrt(mean(d90.^2)))
    @printf("  Max |Diff|: %.4f dB\n",  maximum(abs.(d90)))
    
    # Typical angles
    println("  --- Typical Angles (Phi=0) ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("  θ=%+7.1f° | EMMoMSuite=%+8.2f | Legacy=%+8.2f | Diff=%+6.3f dB\n",
                td, RCS_dBsm[i,1], df.RCS_Phi0_dBsm[i], d0[i])
        end
    end
    println("  --- Typical Angles (Phi=90) ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("  θ=%+7.1f° | EMMoMSuite=%+8.2f | Legacy=%+8.2f | Diff=%+6.3f dB\n",
                td, RCS_dBsm[i,2], df.RCS_Phi90_dBsm[i], d90[i])
        end
    end
    
    return (
        mean_phi0=mean(d0), rmse_phi0=sqrt(mean(d0.^2)), max_phi0=maximum(abs.(d0)),
        mean_phi90=mean(d90), rmse_phi90=sqrt(mean(d90.^2)), max_phi90=maximum(abs.(d90))
    )
end

# Utility: compare two RCS arrays
function compare_rcs_arrays(RCS_A, RCS_B, θs; labelA="A", labelB="B")
    d0  = RCS_A[:, 1] .- RCS_B[:, 1]
    d90 = RCS_A[:, 2] .- RCS_B[:, 2]
    
    println("\n  --- $labelA vs $labelB (Phi=0) ---")
    @printf("  Mean Diff : %+.4f dB\n", mean(d0))
    @printf("  RMSE      : %.4f dB\n",  sqrt(mean(d0.^2)))
    @printf("  Max |Diff|: %.4f dB\n",  maximum(abs.(d0)))
    
    println("  --- $labelA vs $labelB (Phi=90) ---")
    @printf("  Mean Diff : %+.4f dB\n", mean(d90))
    @printf("  RMSE      : %.4f dB\n",  sqrt(mean(d90.^2)))
    @printf("  Max |Diff|: %.4f dB\n",  maximum(abs.(d90)))
    
    # Typical angles
    println("  --- Typical Angles (Phi=0) ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("  θ=%+7.1f° | %s=%+8.2f | %s=%+8.2f | Diff=%+6.3f dB\n",
                td, labelA, RCS_A[i,1], labelB, RCS_B[i,1], d0[i])
        end
    end
    
    return (
        mean_phi0=mean(d0), rmse_phi0=sqrt(mean(d0.^2)), max_phi0=maximum(abs.(d0)),
        mean_phi90=mean(d90), rmse_phi90=sqrt(mean(d90.^2)), max_phi90=maximum(abs.(d90))
    )
end

println("=" ^ 70)
println("EMMoMSuite Comprehensive Accuracy Benchmark (v2)")
println("=" ^ 70)
println("Threads: $(Threads.nthreads())")
println()

θs_721 = collect(LinRange(-π, π, 721))
ϕs_2 = [0.0, π/2]

# ============================================================================
# 1. SEFIE Direct (Jet, 100 MHz, EFIE)
# ============================================================================
println("\n" * "=" ^ 70)
println("1. SEFIE Direct (Jet, 100 MHz, EFIE)")
println("=" ^ 70)

sefie_direct_I = nothing
sefie_direct_RCS = nothing
sefie_direct_stats = nothing
sefie_direct_timing = nothing

let
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")
    
    efie = EFIE(freq)
    
    t_asm = @elapsed Z = assemble_impedance_matrix(efie, basis)
    println("Assembly: $(round(t_asm, digits=2))s")
    
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    t_slv = @elapsed I_coeff = solve!(LUSolver(), Z, V)
    println("LU Solve: $(round(t_slv, digits=2))s")
    println("Total: $(round(t_asm + t_slv, digits=2))s")
    
    RCS_res = radarCrossSection(θs_721, ϕs_2, I_coeff, basis)
    RCS_dBsm = 10 * log10.(RCS_res[2])
    
    stats = compare_rcs(RCS_dBsm, θs_721, joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv"); label="SEFIE Direct")
    
    global sefie_direct_I = I_coeff
    global sefie_direct_RCS = RCS_dBsm
    global sefie_direct_stats = stats
    global sefie_direct_timing = (N=N, t_asm=t_asm, t_slv=t_slv)
end

GC.gc()  # Free Z matrix memory

# ============================================================================
# 2. SCFIE Direct (Jet, 100 MHz, CFIE alpha=0.5) — vs SEFIE Direct
# ============================================================================
println("\n" * "=" ^ 70)
println("2. SCFIE Direct (Jet, 100 MHz, CFIE alpha=0.5)")
println("=" ^ 70)

scfie_direct_jet_I = nothing
scfie_direct_jet_RCS = nothing
scfie_direct_jet_timing = nothing
scfie_vs_sefie_stats = nothing

let
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")
    
    cfie = CFIE(freq, 0.5)
    
    t_asm = @elapsed Z = assemble_impedance_matrix(cfie, basis)
    println("Assembly: $(round(t_asm, digits=2))s")
    
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    
    t_slv = @elapsed I_coeff = Z \ V
    println("LU Solve: $(round(t_slv, digits=2))s")
    println("Total: $(round(t_asm + t_slv, digits=2))s")
    
    RCS_res = radarCrossSection(θs_721, ϕs_2, I_coeff, basis)
    RCS_dBsm = 10 * log10.(abs.(RCS_res[2]))
    
    # Compare with SEFIE Direct
    stats_vs_sefie = compare_rcs_arrays(RCS_dBsm, sefie_direct_RCS, θs_721;
        labelA="CFIE", labelB="EFIE")
    
    # Also compare with Legacy SEFIE baseline
    stats_vs_legacy = compare_rcs(RCS_dBsm, θs_721, joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv"); label="CFIE Direct vs Legacy EFIE")
    
    # Coefficient comparison
    coeff_diff = norm(I_coeff - sefie_direct_I) / norm(sefie_direct_I) * 100
    @printf("  Coeff Relative Diff (CFIE vs EFIE): %.2f%%\n", coeff_diff)
    
    global scfie_direct_jet_I = I_coeff
    global scfie_direct_jet_RCS = RCS_dBsm
    global scfie_direct_jet_timing = (N=N, t_asm=t_asm, t_slv=t_slv)
    global scfie_vs_sefie_stats = stats_vs_sefie
end

GC.gc()

# ============================================================================
# 3. SEFIE MLFMA (Jet, 100 MHz, EFIE + GMRES)
# ============================================================================
println("\n" * "=" ^ 70)
println("3. SEFIE MLFMA (Jet, 100 MHz, EFIE + GMRES)")
println("=" ^ 70)

sefie_mlfma_stats = nothing
sefie_mlfma_timing = nothing
sefie_mlfma_vs_direct = nothing

let
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    lambda = 299792458.0 / freq
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")
    
    efie = EFIE(freq)
    leaf_size = 0.35 * lambda
    
    t_mlfma = @elapsed Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
    println("MLFMA Setup: $(round(t_mlfma, digits=2))s")
    
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    P_near = lu(Z_mlfma.Z_near)
    P = LUPreconditioner(P_near)
    solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
    
    t_gmres = @elapsed I_mlfma = solve!(solver, Z_mlfma, V, Pl=P)
    println("GMRES Solve: $(round(t_gmres, digits=2))s")
    println("Total: $(round(t_mlfma + t_gmres, digits=2))s")
    
    RCS_res = radarCrossSection(θs_721, ϕs_2, I_mlfma, basis)
    RCS_mlfma_dBsm = 10 * log10.(RCS_res[2])
    
    stats_legacy = compare_rcs(RCS_mlfma_dBsm, θs_721, joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv"); label="SEFIE MLFMA vs Legacy")
    
    # vs Direct
    d0_vs_direct = RCS_mlfma_dBsm[:, 1] .- sefie_direct_RCS[:, 1]
    d90_vs_direct = RCS_mlfma_dBsm[:, 2] .- sefie_direct_RCS[:, 2]
    coeff_diff = norm(I_mlfma - sefie_direct_I) / norm(sefie_direct_I) * 100
    
    println("\n  --- SEFIE MLFMA vs SEFIE Direct ---")
    @printf("  Coeff Relative Diff: %.2f%%\n", coeff_diff)
    @printf("  Phi=0  Mean Diff: %+.4f dB, RMSE: %.4f dB\n", mean(d0_vs_direct), sqrt(mean(d0_vs_direct.^2)))
    @printf("  Phi=90 Mean Diff: %+.4f dB, RMSE: %.4f dB\n", mean(d90_vs_direct), sqrt(mean(d90_vs_direct.^2)))
    
    global sefie_mlfma_stats = stats_legacy
    global sefie_mlfma_timing = (N=N, t_mlfma=t_mlfma, t_gmres=t_gmres, leaf=leaf_size/lambda)
    global sefie_mlfma_vs_direct = (
        coeff_diff_pct=coeff_diff,
        rmse_phi0=sqrt(mean(d0_vs_direct.^2)),
        rmse_phi90=sqrt(mean(d90_vs_direct.^2))
    )
end

GC.gc()

# ============================================================================
# 4. SCFIE MLFMA (Sphere, 600 MHz, CFIE alpha=0.5 + GMRES) vs Legacy
# ============================================================================
println("\n" * "=" ^ 70)
println("4. SCFIE MLFMA (Sphere, 600 MHz, CFIE alpha=0.5 + GMRES)")
println("=" ^ 70)

scfie_mlfma_stats = nothing
scfie_mlfma_timing = nothing

let
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "sphere_600MHz.nas")
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 6e8
    set_frequency!(freq)
    lambda = 299792458.0 / freq
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")
    
    cfie = CFIE(freq, 0.5)
    leaf_size = 0.25 * lambda
    
    t_mlfma = @elapsed Z_mlfma = MLFMAOperator(cfie, basis, leaf_size)
    println("MLFMA Setup: $(round(t_mlfma, digits=2))s")
    
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    
    P_near = lu(Z_mlfma.Z_near)
    P = LUPreconditioner(P_near)
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    
    t_gmres = @elapsed I_mlfma = solve!(solver, Z_mlfma, V, Pl=P)
    println("GMRES Solve: $(round(t_gmres, digits=2))s")
    println("Total: $(round(t_mlfma + t_gmres, digits=2))s")
    
    RCS_res = radarCrossSection(θs_721, ϕs_2, I_mlfma, basis)
    RCS_mlfma_dBsm = 10 * log10.(abs.(RCS_res[2]))
    
    stats_legacy = compare_rcs(RCS_mlfma_dBsm, θs_721, joinpath(LEGACY_BASELINE_DIR, "SCFIE_Direct_Sphere.csv"); label="SCFIE MLFMA vs Legacy")
    
    global scfie_mlfma_stats = stats_legacy
    global scfie_mlfma_timing = (N=N, t_mlfma=t_mlfma, t_gmres=t_gmres, leaf=leaf_size/lambda)
end

GC.gc()

# ============================================================================
# SUMMARY
# ============================================================================
println("\n" * "=" ^ 70)
println("BENCHMARK SUMMARY")
println("=" ^ 70)

println("\n--- 1. SEFIE Direct (Jet 100MHz, N=$(sefie_direct_timing.N)) ---")
@printf("  Timing: Assembly=%.2fs, Solve=%.2fs, Total=%.2fs\n", 
    sefie_direct_timing.t_asm, sefie_direct_timing.t_slv,
    sefie_direct_timing.t_asm + sefie_direct_timing.t_slv)
if sefie_direct_stats !== nothing
    @printf("  vs Legacy Phi=0:  Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        sefie_direct_stats.mean_phi0, sefie_direct_stats.rmse_phi0, sefie_direct_stats.max_phi0)
    @printf("  vs Legacy Phi=90: Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        sefie_direct_stats.mean_phi90, sefie_direct_stats.rmse_phi90, sefie_direct_stats.max_phi90)
end

println("\n--- 2. SCFIE Direct (Jet 100MHz CFIE, N=$(scfie_direct_jet_timing.N)) ---")
@printf("  Timing: Assembly=%.2fs, Solve=%.2fs, Total=%.2fs\n",
    scfie_direct_jet_timing.t_asm, scfie_direct_jet_timing.t_slv,
    scfie_direct_jet_timing.t_asm + scfie_direct_jet_timing.t_slv)
if scfie_vs_sefie_stats !== nothing
    @printf("  vs EFIE Direct Phi=0:  Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        scfie_vs_sefie_stats.mean_phi0, scfie_vs_sefie_stats.rmse_phi0, scfie_vs_sefie_stats.max_phi0)
    @printf("  vs EFIE Direct Phi=90: Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        scfie_vs_sefie_stats.mean_phi90, scfie_vs_sefie_stats.rmse_phi90, scfie_vs_sefie_stats.max_phi90)
end

println("\n--- 3. SEFIE MLFMA (Jet 100MHz, N=$(sefie_mlfma_timing.N)) ---")
@printf("  Timing: MLFMA=%.2fs, GMRES=%.2fs, Total=%.2fs (leaf=%.2fλ)\n",
    sefie_mlfma_timing.t_mlfma, sefie_mlfma_timing.t_gmres,
    sefie_mlfma_timing.t_mlfma + sefie_mlfma_timing.t_gmres,
    sefie_mlfma_timing.leaf)
if sefie_mlfma_stats !== nothing
    @printf("  vs Legacy Phi=0:  Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        sefie_mlfma_stats.mean_phi0, sefie_mlfma_stats.rmse_phi0, sefie_mlfma_stats.max_phi0)
    @printf("  vs Legacy Phi=90: Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        sefie_mlfma_stats.mean_phi90, sefie_mlfma_stats.rmse_phi90, sefie_mlfma_stats.max_phi90)
end
if sefie_mlfma_vs_direct !== nothing
    @printf("  vs Direct: Coeff=%.2f%%, RMSE_phi0=%.4f dB, RMSE_phi90=%.4f dB\n",
        sefie_mlfma_vs_direct.coeff_diff_pct,
        sefie_mlfma_vs_direct.rmse_phi0, sefie_mlfma_vs_direct.rmse_phi90)
end

println("\n--- 4. SCFIE MLFMA (Sphere 600MHz, N=$(scfie_mlfma_timing.N)) ---")
@printf("  Timing: MLFMA=%.2fs, GMRES=%.2fs, Total=%.2fs (leaf=%.2fλ)\n",
    scfie_mlfma_timing.t_mlfma, scfie_mlfma_timing.t_gmres,
    scfie_mlfma_timing.t_mlfma + scfie_mlfma_timing.t_gmres,
    scfie_mlfma_timing.leaf)
if scfie_mlfma_stats !== nothing
    @printf("  vs Legacy Phi=0:  Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        scfie_mlfma_stats.mean_phi0, scfie_mlfma_stats.rmse_phi0, scfie_mlfma_stats.max_phi0)
    @printf("  vs Legacy Phi=90: Mean=%+.4f RMSE=%.4f Max=%.4f dB\n",
        scfie_mlfma_stats.mean_phi90, scfie_mlfma_stats.rmse_phi90, scfie_mlfma_stats.max_phi90)
end

println("\n\nBenchmark complete.")
