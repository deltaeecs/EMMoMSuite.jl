# Test SEFIE Direct RCS vs Legacy (corrected dB conversion)
using EMSuite
using LinearAlgebra
using Printf
using Statistics
using CSV, DataFrames

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")
const LEGACY_BASELINE_DIR = joinpath(@__DIR__, "..", "test_results", "legacy_baseline")

println("=" ^ 60)
println("SEFIE Direct vs Legacy (Corrected)")
println("=" ^ 60)

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
println("Solve: $(round(t_slv, digits=2))s")
println("Total: $(round(t_asm + t_slv, digits=2))s")

# RCS — use RCS_res[2] (linear) and convert manually
θs = collect(LinRange(-π, π, 721))
ϕs = [0.0, π/2]
RCS_res = radarCrossSection(θs, ϕs, I_coeff, basis)
RCS_linear = RCS_res[2]   # ← INDEX 2: total linear (m²)
RCS_dBsm = 10 * log10.(RCS_linear)

# Also check: does RCS_res[3] == 10*log10(RCS_res[2])?
RCS_from_index3 = RCS_res[3]
println("\nSanity check: RCS_res[2][1,1] = $(RCS_linear[1,1])")
println("  10*log10(RCS_res[2][1,1]) = $(10*log10(RCS_linear[1,1]))")
println("  RCS_res[3][1,1] = $(RCS_from_index3[1,1])")
println("  Difference: $(abs(10*log10(RCS_linear[1,1]) - RCS_from_index3[1,1]))")

# Compare with Legacy
baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
if isfile(baseline_file)
    df = CSV.read(baseline_file, DataFrame)
    d0 = RCS_dBsm[:, 1] .- df.RCS_Phi0_dBsm
    d90 = RCS_dBsm[:, 2] .- df.RCS_Phi90_dBsm

    println("\n--- Phi=0 vs Legacy ---")
    println("Mean Diff: $(round(mean(d0), digits=4)) dB")
    println("RMSE: $(round(sqrt(mean(d0.^2)), digits=4)) dB")
    println("Max |Diff|: $(round(maximum(abs.(d0)), digits=4)) dB")

    println("\n--- Phi=90 vs Legacy ---")
    println("Mean Diff: $(round(mean(d90), digits=4)) dB")
    println("RMSE: $(round(sqrt(mean(d90.^2)), digits=4)) dB")
    println("Max |Diff|: $(round(maximum(abs.(d90)), digits=4)) dB")

    println("\n--- 典型角度 (Phi=0) ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("θ=%+7.1f° | EMSuite=%+8.2f | Legacy=%+8.2f | Diff=%+6.3f dB\n",
                td, RCS_dBsm[i,1], df.RCS_Phi0_dBsm[i], d0[i])
        end
    end

    println("\n--- 典型角度 (Phi=90) ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("θ=%+7.1f° | EMSuite=%+8.2f | Legacy=%+8.2f | Diff=%+6.3f dB\n",
                td, RCS_dBsm[i,2], df.RCS_Phi90_dBsm[i], d90[i])
        end
    end
else
    println("Legacy baseline not found")
end
