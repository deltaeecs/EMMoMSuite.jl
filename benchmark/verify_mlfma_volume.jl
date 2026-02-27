# Phase 10: D3 + E3 — MLFMA self-consistency for volume equations
# D3: VEFIE MLFMA vs VEFIE Direct on Tetra.nas
# E3: VSEFIE MLFMA vs VSEFIE Direct on TriTetra.nas
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

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

# ===================================================================
#  D3: V-EFIE MLFMA — Tetra 300 MHz
# ===================================================================
function test_D3()
    println("\n" * "=" ^ 60)
    println("  D3: V-EFIE MLFMA vs Direct — Tetra 300MHz")
    println("=" ^ 60)

    freq = 300e6
    c0 = 299792458.0
    lambda = c0 / freq
    set_frequency!(freq)
    mesh_file = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/Tetra.nas")
    mesh = read_nas_mesh(mesh_file)

    # Scale mm → m
    max_coord = maximum(abs.(mesh.node))
    if max_coord > 10.0
        mesh.node .*= 0.001
    end
    min_c = minimum(mesh.node, dims=2)
    max_c = maximum(mesh.node, dims=2)
    center = (min_c + max_c) / 2
    mesh.node .-= center

    basis_d = SWGBasis(mesh)
    basis_m = SWGBasis(mesh)
    N = num_basis(basis_d)
    println("  N = $N")

    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))
    vefie = VEFIE(freq, permittivities)

    # Direct solve
    println("  --- Direct ---")
    t_asm = @elapsed Z = assemble_impedance_matrix(vefie, basis_d, permittivities)
    println("  Assembly: $(round(t_asm, digits=2))s")

    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0])
    V_d = excitation_vector(vefie, source, basis_d, permittivities)

    t_lu = @elapsed D_direct = solve!(LUSolver(), Z, V_d)
    println("  LU: $(round(t_lu, digits=2))s")

    # MLFMA solve
    println("  --- MLFMA ---")
    leaf_size = 0.35 * lambda
    t_setup = @elapsed mlfma_op = MLFMAOperator(vefie, basis_m, leaf_size)
    println("  MLFMA setup: $(round(t_setup, digits=2))s")

    V_m = excitation_vector(vefie, source, basis_m, permittivities)
    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    t_gmres = @elapsed D_mlfma = solve!(solver, mlfma_op, V_m; Pl=P)
    println("  GMRES: $(round(t_gmres, digits=2))s")

    # Coefficient comparison
    rel_err = norm(D_direct - D_mlfma) / norm(D_direct)
    println("  ||D_direct - D_mlfma|| / ||D_direct|| = $(round(rel_err*100, digits=4))%")

    # RCS comparison
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]

    RCS_d = radarCrossSection(θs, ϕs, D_direct, basis_d, permittivities)
    RCS_m = radarCrossSection(θs, ϕs, D_mlfma, basis_m, permittivities)
    dBsm_d = RCS_d[3]
    dBsm_m = RCS_m[3]

    diff = dBsm_m[:, 1] .- dBsm_d[:, 1]
    rmse = sqrt(mean(diff .^ 2))
    @printf("  E-plane: Mean Diff=%.4f dB, RMSE=%.4f dB, Max|Diff|=%.4f dB\n",
        mean(diff), rmse, maximum(abs.(diff)))

    println("\n  --- Key Angles (E-plane) ---")
    @printf("  %5s  %12s  %12s  %10s\n", "θ(°)", "Direct", "MLFMA", "Diff(dB)")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        if idx <= length(θs)
            @printf("  %5d  %12.4f  %12.4f  %10.4f\n",
                deg, dBsm_d[idx, 1], dBsm_m[idx, 1], diff[idx])
        end
    end

    if rmse < 2.0
        println("  ✅ D3 PASS: RMSE = $(round(rmse, digits=4)) dB < 2.0 dB")
    else
        println("  ❌ D3 FAIL: RMSE = $(round(rmse, digits=4)) dB ≥ 2.0 dB")
    end

    println("  Timing: Direct=$(round(t_asm+t_lu,digits=1))s, MLFMA=$(round(t_setup+t_gmres,digits=1))s")
    return rmse
end

# ===================================================================
#  E3: VS-EFIE MLFMA — TriTetra 2 GHz
# ===================================================================
function test_E3()
    println("\n" * "=" ^ 60)
    println("  E3: VS-EFIE MLFMA vs Direct — TriTetra 2GHz")
    println("=" ^ 60)

    freq = 2e9
    c0 = 299792458.0
    lambda = c0 / freq
    set_frequency!(freq)
    mesh_file = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)

    surf_basis_d = RWGBasis(surf_mesh)
    vol_basis_d = SWGBasis(vol_mesh)
    surf_basis_m = RWGBasis(surf_mesh)
    vol_basis_m = SWGBasis(vol_mesh)
    n_surf = num_basis(surf_basis_d)
    n_vol = num_basis(vol_basis_d)
    N = n_surf + n_vol
    println("  N = $N (RWG=$n_surf, SWG=$n_vol)")

    eps_r = 2.0 * (1 - 0.001im)
    permittivities = fill(eps_r, num_elements(vol_mesh))
    scfie = SCFIE(freq, permittivities; alpha=1.0)

    # Direct solve
    println("  --- Direct ---")
    t_asm = @elapsed Z = assemble_impedance_matrix(scfie, surf_basis_d, vol_basis_d)
    println("  Assembly: $(round(t_asm, digits=2))s")

    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V_d = excitation_vector(source, surf_basis_d, vol_basis_d)

    t_lu = @elapsed x_direct = solve!(LUSolver(), Z, V_d)
    println("  LU: $(round(t_lu, digits=2))s")

    I_d = x_direct[1:n_surf]; D_d = x_direct[n_surf+1:end]

    # MLFMA solve
    println("  --- MLFMA ---")
    leaf_size = 0.25 * lambda
    t_setup = @elapsed mlfma_op = MLFMAOperator(scfie, [surf_basis_m, vol_basis_m], leaf_size)
    println("  MLFMA setup: $(round(t_setup, digits=2))s")

    V_m = excitation_vector(source, surf_basis_m, vol_basis_m)
    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    t_gmres = @elapsed x_mlfma = solve!(solver, mlfma_op, V_m; Pl=P)
    println("  GMRES: $(round(t_gmres, digits=2))s")

    I_m = x_mlfma[1:n_surf]; D_m = x_mlfma[n_surf+1:end]

    # Coefficient comparison
    rel_err = norm(x_direct - x_mlfma) / norm(x_direct)
    println("  ||x_direct - x_mlfma|| / ||x_direct|| = $(round(rel_err*100, digits=4))%")

    # RCS comparison
    k0 = get_k0()
    eta0 = get_eta0()
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]
    Nθ = length(θs)
    Nϕ = length(ϕs)

    θsobsInfo = [∠Info(θ) for θ in θs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs]

    dBsm_d = zeros(Float64, Nθ, Nϕ)
    dBsm_m = zeros(Float64, Nθ, Nϕ)
    for ip in 1:Nϕ, it in 1:Nθ
        r_info = r̂θϕInfo(θsobsInfo[it], ϕsobsInfo[ip])
        factor = (k0 * eta0)^2 / (4π)

        Nθϕ_s = radiation_integral_rwg(r_info, surf_basis_d, I_d)
        Nθϕ_v = radiation_integral_swg(r_info, vol_basis_d, D_d, permittivities)
        Nθϕ = Nθϕ_s .+ Nθϕ_v
        rcs = factor * (abs2(Nθϕ[1]) + abs2(Nθϕ[2]))
        dBsm_d[it, ip] = 10 * log10(rcs)

        Nθϕ_s2 = radiation_integral_rwg(r_info, surf_basis_m, I_m)
        Nθϕ_v2 = radiation_integral_swg(r_info, vol_basis_m, D_m, permittivities)
        Nθϕ2 = Nθϕ_s2 .+ Nθϕ_v2
        rcs2 = factor * (abs2(Nθϕ2[1]) + abs2(Nθϕ2[2]))
        dBsm_m[it, ip] = 10 * log10(rcs2)
    end

    diff = dBsm_m[:, 1] .- dBsm_d[:, 1]
    rmse = sqrt(mean(diff .^ 2))
    @printf("  E-plane: Mean Diff=%.4f dB, RMSE=%.4f dB, Max|Diff|=%.4f dB\n",
        mean(diff), rmse, maximum(abs.(diff)))

    println("\n  --- Key Angles (E-plane) ---")
    @printf("  %5s  %12s  %12s  %10s\n", "θ(°)", "Direct", "MLFMA", "Diff(dB)")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        if idx <= Nθ
            @printf("  %5d  %12.4f  %12.4f  %10.4f\n",
                deg, dBsm_d[idx, 1], dBsm_m[idx, 1], diff[idx])
        end
    end

    if rmse < 2.0
        println("  ✅ E3 PASS: RMSE = $(round(rmse, digits=4)) dB < 2.0 dB")
    else
        println("  ❌ E3 FAIL: RMSE = $(round(rmse, digits=4)) dB ≥ 2.0 dB")
    end

    println("  Timing: Direct=$(round(t_asm+t_lu,digits=1))s, MLFMA=$(round(t_setup+t_gmres,digits=1))s")
    return rmse
end

# ===================================================================
# Run all
# ===================================================================
println("=" ^ 60)
println("  Phase 10: MLFMA Self-Consistency Tests (D3, E3)")
println("  Criterion: each RMSE < 2.0 dB vs Direct")
println("=" ^ 60)

rmse_D3 = test_D3()
rmse_E3 = test_E3()

println("\n" * "=" ^ 60)
println("  SUMMARY")
println("=" ^ 60)
@printf("  D3 V-EFIE MLFMA: RMSE = %.4f dB %s\n", rmse_D3, rmse_D3 < 2.0 ? "✅" : "❌")
@printf("  E3 VS-EFIE MLFMA: RMSE = %.4f dB %s\n", rmse_E3, rmse_E3 < 2.0 ? "✅" : "❌")
println("=" ^ 60)
