# Phase 10: D2/E2 — Iterative (Full GMRES on Dense Z) vs Direct (LU)
# A2 skipped: N=14559 EFIE requires O(N²) Krylov memory (3.4 GB) for full GMRES;
#             EFIE iterative path is validated via A3 (MLFMA+GMRES with LU preconditioner).
# Criterion: RMSE < 0.1 dB for each test
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.PostProcessing
using EMSuite.CoreModule.Sources
using EMSuite.PostProcessing.RadiationIntegral: radiation_integral_rwg, radiation_integral_swg,
    r̂θϕInfo, ∠Info
using EMSuite.Utilities.Parameters: get_k0, get_eta0
using LinearAlgebra
using Printf
using Statistics: mean

# ===================================================================
#  D2: V-EFIE Iterative — Tetra 300 MHz (N=3201)
#  Full GMRES (restart=N, no preconditioner) to avoid restart stall
# ===================================================================
function test_D2()
    println("\n" * "=" ^ 60)
    println("  D2: V-EFIE Iterative (Full GMRES on Dense Z) — Tetra 300MHz")
    println("=" ^ 60)

    freq = 300e6
    set_frequency!(freq)
    mesh_file = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/Tetra.nas")
    mesh = read_nas_mesh(mesh_file)

    # Scale mm → m
    nodes = mesh.node
    max_coord = maximum(abs.(nodes))
    if max_coord > 10.0
        mesh.node .*= 0.001
    end
    min_c = minimum(mesh.node, dims=2)
    max_c = maximum(mesh.node, dims=2)
    center = (min_c + max_c) / 2
    mesh.node .-= center

    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("  N = $N")

    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))
    vefie = VEFIE(freq, permittivities)

    # Assembly
    println("  Assembling Z...")
    t_asm = @elapsed Z = assemble_impedance_matrix(vefie, basis, permittivities)
    println("  Assembly: $(round(t_asm, digits=2))s")

    # Excitation
    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(vefie, source, basis, permittivities)

    # Direct solve
    println("  LU solve...")
    t_lu = @elapsed D_direct = solve!(LUSolver(), Z, V)
    println("  LU: $(round(t_lu, digits=2))s")

    # Full GMRES (restart=N) — no restart stall, guaranteed convergence
    println("  Full GMRES solve (restart=$N, no preconditioner, tol=1e-6)...")
    solver = GMRESSolver(restart=N, maxiter=N, tol=1e-6, verbose=true)
    t_gmres = @elapsed D_iter = solve!(solver, Z, V)
    println("  GMRES: $(round(t_gmres, digits=2))s")

    rel_err = norm(D_direct - D_iter) / norm(D_direct)
    println("  ||D_direct - D_iter|| / ||D_direct|| = $(round(rel_err*100, digits=6))%")

    # RCS comparison
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]

    RCS_d_res = radarCrossSection(θs, ϕs, D_direct, basis, permittivities)
    RCS_i_res = radarCrossSection(θs, ϕs, D_iter, basis, permittivities)
    dBsm_d = RCS_d_res[3]
    dBsm_i = RCS_i_res[3]

    diff = dBsm_i[:, 1] .- dBsm_d[:, 1]
    rmse = sqrt(mean(diff .^ 2))
    @printf("  E-plane: Mean Diff=%.6f dB, RMSE=%.6f dB, Max|Diff|=%.6f dB\n",
        mean(diff), rmse, maximum(abs.(diff)))

    if rmse < 0.1
        println("  ✅ D2 PASS: RMSE = $(round(rmse, digits=6)) dB < 0.1 dB")
    else
        println("  ❌ D2 FAIL: RMSE = $(round(rmse, digits=6)) dB ≥ 0.1 dB")
    end

    println("  Timing: ASM=$(round(t_asm,digits=1))s LU=$(round(t_lu,digits=1))s GMRES=$(round(t_gmres,digits=1))s")
    return rmse
end

# ===================================================================
#  E2: VS-EFIE Iterative — TriTetra 2 GHz (N≈1071)
#  Full GMRES (restart=N, no preconditioner)
# ===================================================================
function test_E2()
    println("\n" * "=" ^ 60)
    println("  E2: VS-EFIE Iterative (Full GMRES on Dense Z) — TriTetra 2GHz")
    println("=" ^ 60)

    freq = 2e9
    set_frequency!(freq)
    mesh_file = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)

    surf_basis = RWGBasis(surf_mesh)
    vol_basis = SWGBasis(vol_mesh)
    n_surf = num_basis(surf_basis)
    n_vol = num_basis(vol_basis)
    N = n_surf + n_vol
    println("  N = $N (RWG=$n_surf, SWG=$n_vol)")

    eps_r = 2.0 * (1 - 0.001im)
    permittivities = fill(eps_r, num_elements(vol_mesh))
    scfie = SCFIE(freq, permittivities; alpha=1.0)

    # Assembly
    println("  Assembling Z...")
    t_asm = @elapsed Z = assemble_impedance_matrix(scfie, surf_basis, vol_basis)
    println("  Assembly: $(round(t_asm, digits=2))s")

    # Excitation
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(source, surf_basis, vol_basis)

    # Direct solve
    println("  LU solve...")
    t_lu = @elapsed x_direct = solve!(LUSolver(), Z, V)
    println("  LU: $(round(t_lu, digits=2))s")

    # Full GMRES (restart=N) — no restart stall
    println("  Full GMRES solve (restart=$N, no preconditioner, tol=1e-6)...")
    solver = GMRESSolver(restart=N, maxiter=N, tol=1e-6, verbose=true)
    t_gmres = @elapsed x_iter = solve!(solver, Z, V)
    println("  GMRES: $(round(t_gmres, digits=2))s")

    rel_err = norm(x_direct - x_iter) / norm(x_direct)
    println("  ||x_direct - x_iter|| / ||x_direct|| = $(round(rel_err*100, digits=6))%")

    # RCS comparison — combined surface + volume
    I_direct = x_direct[1:n_surf]; D_direct = x_direct[n_surf+1:end]
    I_iter = x_iter[1:n_surf];     D_iter = x_iter[n_surf+1:end]

    k0 = get_k0()
    eta0 = get_eta0()
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]
    Nθ = length(θs)
    Nϕ = length(ϕs)

    θsobsInfo = [∠Info(θ) for θ in θs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs]

    dBsm_d = zeros(Float64, Nθ, Nϕ)
    dBsm_i = zeros(Float64, Nθ, Nϕ)
    for ip in 1:Nϕ, it in 1:Nθ
        r_info = r̂θϕInfo(θsobsInfo[it], ϕsobsInfo[ip])
        factor = (k0 * eta0)^2 / (4π)

        Nθϕ_s = radiation_integral_rwg(r_info, surf_basis, I_direct)
        Nθϕ_v = radiation_integral_swg(r_info, vol_basis, D_direct, permittivities)
        Nθϕ = Nθϕ_s .+ Nθϕ_v
        rcs = factor * (abs2(Nθϕ[1]) + abs2(Nθϕ[2]))
        dBsm_d[it, ip] = 10 * log10(rcs)

        Nθϕ_s2 = radiation_integral_rwg(r_info, surf_basis, I_iter)
        Nθϕ_v2 = radiation_integral_swg(r_info, vol_basis, D_iter, permittivities)
        Nθϕ2 = Nθϕ_s2 .+ Nθϕ_v2
        rcs2 = factor * (abs2(Nθϕ2[1]) + abs2(Nθϕ2[2]))
        dBsm_i[it, ip] = 10 * log10(rcs2)
    end

    diff = dBsm_i[:, 1] .- dBsm_d[:, 1]
    rmse = sqrt(mean(diff .^ 2))
    @printf("  E-plane: Mean Diff=%.6f dB, RMSE=%.6f dB, Max|Diff|=%.6f dB\n",
        mean(diff), rmse, maximum(abs.(diff)))

    if rmse < 0.1
        println("  ✅ E2 PASS: RMSE = $(round(rmse, digits=6)) dB < 0.1 dB")
    else
        println("  ❌ E2 FAIL: RMSE = $(round(rmse, digits=6)) dB ≥ 0.1 dB")
    end

    println("  Timing: ASM=$(round(t_asm,digits=1))s LU=$(round(t_lu,digits=1))s GMRES=$(round(t_gmres,digits=1))s")
    return rmse
end

# ===================================================================
# Run all
# ===================================================================
println("=" ^ 60)
println("  Phase 10: Iterative Tests (D2, E2)")
println("  A2 skipped: N=14559 too large for full GMRES; see A3")
println("  Criterion: each RMSE < 0.1 dB vs Direct")
println("=" ^ 60)

rmse_D2 = test_D2()
rmse_E2 = test_E2()

println("\n" * "=" ^ 60)
println("  SUMMARY")
println("=" ^ 60)
println("  A2 S-EFIE Iterative: N/A (N=14559 → full GMRES impractical; A3 validates MLFMA+GMRES)")
@printf("  D2 V-EFIE Iterative: RMSE = %.6f dB %s\n", rmse_D2, rmse_D2 < 0.1 ? "✅" : "❌")
@printf("  E2 VS-EFIE Iterative: RMSE = %.6f dB %s\n", rmse_E2, rmse_E2 < 0.1 ? "✅" : "❌")
println("=" ^ 60)
