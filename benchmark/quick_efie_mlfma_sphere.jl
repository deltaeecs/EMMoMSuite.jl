"""
Quick check: EFIE MLFMA vs Legacy Direct on sphere_600MHz.
If this also has large RMSE, the problem is geometry/octree-specific, not CFIE-specific.
"""

using EMSuite
using LinearAlgebra
using Printf
using Statistics
using CSV, DataFrames

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

println("=" ^ 72)
println("  Quick Check: EFIE MLFMA vs Legacy on Sphere 600MHz")
println("=" ^ 72)

freq = 6e8
lambda = 299792458.0 / freq
set_frequency!(freq)

mesh_file = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "sphere_600MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
basis = RWGBasis(mesh)
N = num_basis(basis)
println("  Unknowns: $N")

efie = EFIE(freq)
source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

# MLFMA
leaf_size = 0.25 * lambda
println("  Leaf size: $(round(leaf_size, digits=4))m = 0.25λ")
t1 = @elapsed mlfma_op = MLFMAOperator(efie, basis, leaf_size)
println("  MLFMA setup: $(round(t1, digits=2))s")

V = excitation_vector(efie, source, basis)
P = LUPreconditioner(lu(mlfma_op.Z_near))
solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
t2 = @elapsed I_mlfma = solve!(solver, mlfma_op, V, Pl=P)
println("  GMRES: $(round(t2, digits=2))s")

# RCS
θs = collect(LinRange(-π, π, 721))
ϕs = [0.0, π/2]
RCS_res = radarCrossSection(θs, ϕs, I_mlfma, basis)
dBsm = 10 * log10.(RCS_res[2])

# Compare with Legacy SCFIE baseline (same mesh) 
# Note: Legacy baseline was from SCFIE Direct, not EFIE Direct. 
# But for relative comparison the RMSE tells us about MLFMA accuracy.
baseline_file = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SCFIE_Direct_Sphere.csv")
if isfile(baseline_file)
    df_base = CSV.read(baseline_file, DataFrame)
    diff = dBsm[:, 1] .- df_base.RCS_Phi0_dBsm

    println("\n--- EFIE MLFMA vs Legacy SCFIE Direct (Phi=0) ---")
    println("  Note: Different equations, expect systematic offset from EFIE vs CFIE difference")
    println("  Mean Diff:   $(round(mean(diff), digits=4)) dB")
    println("  RMSE:        $(round(sqrt(mean(diff.^2)), digits=4)) dB")
    println("  Max |Diff|:  $(round(maximum(abs.(diff)), digits=4)) dB")

    # Typical angles
    println("\n--- 典型角度 ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("  θ=%+7.1f° | Legacy(SCFIE)=%+8.2f | EFIE_MLFMA=%+8.2f | Diff=%+6.3f dB\n",
                td, df_base.RCS_Phi0_dBsm[i], dBsm[i,1], diff[i])
        end
    end
end
