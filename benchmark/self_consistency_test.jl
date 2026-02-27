# Self-consistency: SEFIE Direct vs SEFIE MLFMA (no Legacy comparison)
using EMSuite
using LinearAlgebra
using Printf
using Statistics

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

println("=" ^ 60)
println("Self-Consistency: SEFIE Direct vs SEFIE MLFMA")
println("=" ^ 60)

mesh_file = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)
basis = RWGBasis(mesh)
N = num_basis(basis)
println("Unknowns: $N")

efie = EFIE(freq)

# Direct solve
println("\n--- Direct ---")
t1 = @elapsed Z = assemble_impedance_matrix(efie, basis)
println("Assembly: $(round(t1, digits=2))s")

source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V = excitation_vector(efie, source, basis)

t2 = @elapsed I_direct = solve!(LUSolver(), Z, V)
println("Solve: $(round(t2, digits=2))s")

# MLFMA solve
println("\n--- MLFMA ---")
λ = 299792458.0 / freq
basis_mlfma = RWGBasis(mesh)
t3 = @elapsed mlfma_op = MLFMAOperator(efie, basis_mlfma, 0.35λ)
println("MLFMA setup: $(round(t3, digits=2))s")

V_mlfma = excitation_vector(efie, source, basis_mlfma)
P = LUPreconditioner(lu(mlfma_op.Z_near))
solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
t4 = @elapsed I_mlfma = solve!(solver, mlfma_op, V_mlfma, Pl=P)
println("GMRES: $(round(t4, digits=2))s")

# Compare coefficients
coeff_diff = norm(I_direct - I_mlfma) / norm(I_direct)
println("\n--- Coefficient Comparison ---")
println("||I_direct - I_mlfma|| / ||I_direct|| = $(round(coeff_diff*100, digits=4))%")

# Compare RCS
θs = collect(LinRange(-π, π, 721))
ϕs = [0.0, π/2]

RCS_d = radarCrossSection(θs, ϕs, I_direct, basis)
RCS_m = radarCrossSection(θs, ϕs, I_mlfma, basis_mlfma)

dBsm_d = 10 * log10.(RCS_d[2])
dBsm_m = 10 * log10.(RCS_m[2])

diff = dBsm_d[:, 1] .- dBsm_m[:, 1]
println("\n--- RCS Comparison (Phi=0) ---")
println("Mean |Diff|: $(round(mean(abs.(diff)), digits=4)) dB")
println("Max |Diff|: $(round(maximum(abs.(diff)), digits=4)) dB")
println("RMSE: $(round(sqrt(mean(diff.^2)), digits=4)) dB")

println("\n--- 典型角度 (Phi=0) ---")
for (i, th) in enumerate(θs)
    td = round(th * 180 / π, digits=1)
    if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
        @printf("θ=%+7.1f° | Direct=%+8.2f | MLFMA=%+8.2f | Diff=%+6.3f dB\n",
            td, dBsm_d[i,1], dBsm_m[i,1], diff[i])
    end
end

println("\n--- Timing Summary ---")
println("Direct: Assembly=$(round(t1, digits=1))s + Solve=$(round(t2, digits=1))s = $(round(t1+t2, digits=1))s")
println("MLFMA:  Setup=$(round(t3, digits=1))s + GMRES=$(round(t4, digits=1))s = $(round(t3+t4, digits=1))s")
println("Speedup: $(round((t1+t2)/(t3+t4), digits=2))x")
