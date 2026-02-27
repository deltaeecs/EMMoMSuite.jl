"""
Diagnostic: CFIE MLFMA with alpha=1.0 (pure EFIE through CFIE code path).
This isolates whether the CFIE code path for EFIE is correct.
If alpha=1.0 matches EFIE MLFMA, the issue is specifically in the MFIE term.
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
println("  Diagnostic: CFIE MLFMA alpha sweep on Sphere 600MHz")
println("=" ^ 72)

freq = 6e8
lambda = 299792458.0 / freq
set_frequency!(freq)

mesh_file = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "sphere_600MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)

source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

θs = collect(LinRange(-π, π, 721))
ϕs = [0.0, π/2]

baseline_file = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SCFIE_Direct_Sphere.csv")
RCS_legacy_dB = CSV.read(baseline_file, DataFrame).RCS_Phi0_dBsm

leaf_size = 0.25 * lambda

for alpha in [1.0, 0.9, 0.6]
    println("\n--- CFIE alpha=$alpha ---")
    
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    
    cfie = CFIE(freq, alpha)
    
    t1 = @elapsed mlfma_op = MLFMAOperator(cfie, basis, leaf_size)
    println("  MLFMA setup: $(round(t1, digits=1))s")
    
    V = excitation_vector(cfie, source, basis)
    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=false)
    t2 = @elapsed I = solve!(solver, mlfma_op, V, Pl=P)
    println("  GMRES: $(round(t2, digits=1))s")
    
    RCS_res = radarCrossSection(θs, ϕs, I, basis)
    dBsm = 10 * log10.(RCS_res[2])
    
    diff = dBsm[:, 1] .- RCS_legacy_dB
    rmse = sqrt(mean(diff.^2))
    maxd = maximum(abs.(diff))
    meand = mean(diff)
    
    @printf("  Mean Diff: %+.4f dB | RMSE: %.4f dB | Max: %.4f dB\n", meand, rmse, maxd)
    
    # Show a few angles
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, 0.0, 90.0]
            @printf("    θ=%+7.1f° | Diff=%+6.3f dB\n", td, diff[i])
        end
    end
end
