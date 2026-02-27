# E1: SCFIE (VSEFIE) Direct — EMSuite vs Legacy on TriTetra.nas
# Note: Legacy uses EFIE for surface → set alpha=1.0 in EMSuite SCFIE
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.PostProcessing
using EMSuite.CoreModule.Sources
using LinearAlgebra
using Printf
using Statistics: mean

using EMSuite.PostProcessing.RadiationIntegral: radiation_integral_rwg, radiation_integral_swg,
    r̂θϕInfo, ∠Info
using EMSuite.Utilities.Parameters: get_k0, get_eta0

function verify_vsefie_direct()
    println("==================================================")
    println("   E1: VSEFIE Direct — EMSuite SCFIE vs Legacy    ")
    println("==================================================")

    # 1. Parameters (match Legacy)
    freq = 2e9
    c0 = 299792458.0
    lambda = c0 / freq
    EMSuite.Utilities.Parameters.set_frequency!(freq)
    println("f = $(freq/1e9) GHz, λ = $(round(lambda*1e3, digits=2)) mm")

    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
    @assert isfile(mesh_file) "Mesh not found: $mesh_file"
    println("Loading mesh: TriTetra.nas")
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001) # mm → m

    println("  Surface: $(num_vertices(surf_mesh)) vertices, $(num_elements(surf_mesh)) triangles")
    println("  Volume:  $(num_vertices(vol_mesh)) vertices, $(num_elements(vol_mesh)) tetrahedra")

    # 3. Basis
    println("Setting up RWG + SWG Basis...")
    surf_basis = RWGBasis(surf_mesh)
    vol_basis = SWGBasis(vol_mesh)
    n_surf = num_basis(surf_basis)
    n_vol = num_basis(vol_basis)
    n_total = n_surf + n_vol
    println("  RWG: $n_surf, SWG: $n_vol, Total: $n_total")

    # 4. Permittivity (match Legacy: 2.0*(1-0.001im))
    eps_r = 2.0 * (1 - 0.001im)
    permittivities = fill(eps_r, num_elements(vol_mesh))

    # 5. SCFIE Operator (alpha=1.0 → pure EFIE on surface, like Legacy)
    println("Setting up SCFIE (alpha=1.0 = pure EFIE on surface)...")
    scfie = SCFIE(freq, permittivities; alpha=1.0)

    # 6. Assembly
    println("Assembling Z matrix...")
    t_asm = @elapsed Z = assemble_impedance_matrix(scfie, surf_basis, vol_basis)
    println("  Done in $(round(t_asm, digits=2)) s, size=$(size(Z))")
    println("  ||Z|| = $(norm(Z))")

    # 7. Excitation (plane wave from -z toward +z, x-pol — match Legacy convention)
    # Legacy PlaneWave(π, 0, 0, 1) has k̂=[0,0,+1] (wave from -z toward +z)
    # EMSuite PlaneWave(freq, 0, 0, pol) has k_dir=[0,0,+1] matching Legacy
    println("Computing excitation...")
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(source, surf_basis, vol_basis)
    println("  ||V|| = $(norm(V))")

    # 8. Solve
    println("Solving (LU)...")
    solver = LUSolver()
    t_solve = @elapsed x = solve!(solver, Z, V)
    println("  Done in $(round(t_solve, digits=2)) s, ||x|| = $(norm(x))")

    # Split solution
    I_surf = x[1:n_surf]
    D_vol = x[n_surf+1:end]
    println("  ||I_surf|| = $(norm(I_surf)), ||D_vol|| = $(norm(D_vol))")

    # Diagnostics: Z block norms
    Z_SS = Z[1:n_surf, 1:n_surf]
    Z_SV = Z[1:n_surf, n_surf+1:end]
    Z_VS = Z[n_surf+1:end, 1:n_surf]
    Z_VV = Z[n_surf+1:end, n_surf+1:end]
    println("\n=== Z Block Norms ===")
    println("  ||Z_SS|| = $(norm(Z_SS))")
    println("  ||Z_SV|| = $(norm(Z_SV))")
    println("  ||Z_VS|| = $(norm(Z_VS))")
    println("  ||Z_VV|| = $(norm(Z_VV))")
    println("  ||Z_SV||/||Z_VS|| = $(norm(Z_SV)/norm(Z_VS))")
    
    # Diagnostics: V block norms
    V_s = V[1:n_surf]
    V_v = V[n_surf+1:end]
    println("  ||V_surf|| = $(norm(V_s))")
    println("  ||V_vol||  = $(norm(V_v))")
    
    # Incident direction check
    println("\n=== Incident Direction ===")
    st, ct = sincos(source.theta)
    sp, cp = sincos(source.phi)
    k_dir = [st*cp, st*sp, ct]
    println("  k_dir = $k_dir")
    println("  E_pol = $(source.polarization)")

    # 9. RCS — combine surface + volume contributions
    println("Computing RCS...")
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]

    k0 = get_k0()
    eta0 = get_eta0()
    Nθ = length(θs)
    Nϕ = length(ϕs)
    
    θsobsInfo = [∠Info(θ) for θ in θs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs]
    
    RCS_dB = zeros(Float64, Nθ, Nϕ)
    RCS_surf_dB = zeros(Float64, Nθ, Nϕ)
    RCS_vol_dB = zeros(Float64, Nθ, Nϕ)
    for ip in 1:Nϕ, it in 1:Nθ
        r_info = r̂θϕInfo(θsobsInfo[it], ϕsobsInfo[ip])
        # Surface radiation
        Nθϕ_s = radiation_integral_rwg(r_info, surf_basis, I_surf)
        # Volume radiation
        Nθϕ_v = radiation_integral_swg(r_info, vol_basis, D_vol, permittivities)
        # Combined
        Nθϕ = Nθϕ_s .+ Nθϕ_v
        factor = (k0 * eta0)^2 / (4π)
        rcs = factor * (abs2(Nθϕ[1]) + abs2(Nθϕ[2]))
        RCS_dB[it, ip] = 10 * log10(rcs)
        
        # Surface-only
        rcs_s = factor * (abs2(Nθϕ_s[1]) + abs2(Nθϕ_s[2]))
        RCS_surf_dB[it, ip] = rcs_s > 0 ? 10*log10(rcs_s) : -999.0
        # Volume-only
        rcs_v = factor * (abs2(Nθϕ_v[1]) + abs2(Nθϕ_v[2]))
        RCS_vol_dB[it, ip] = rcs_v > 0 ? 10*log10(rcs_v) : -999.0
    end
    
    println("\n=== Surface vs Volume RCS at key angles (E-plane) ===")
    @printf("  %5s  %12s  %12s  %12s\n", "θ(°)", "Surf-only", "Vol-only", "Combined")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        if idx <= Nθ
            @printf("  %5d  %12.4f  %12.4f  %12.4f\n",
                deg, RCS_surf_dB[idx, 1], RCS_vol_dB[idx, 1], RCS_dB[idx, 1])
        end
    end

    # Save
    outdir = joinpath(@__DIR__, "../test_results/emsuite_verification")
    mkpath(outdir)
    outfile = joinpath(outdir, "VSEFIE_SWG_TriTetra.csv")
    open(outfile, "w") do io
        println(io, "theta_deg,phi_deg,RCS_dBsm")
        for p in 1:length(ϕs)
            for t in 1:length(θs)
                @printf(io, "%.4f,%.4f,%.6f\n",
                    rad2deg(θs[t]), rad2deg(ϕs[p]), RCS_dB[t, p])
            end
        end
    end
    println("  RCS saved to $outfile")

    # 10. Load Legacy baseline
    legacy_file = joinpath(@__DIR__, "../test_results/legacy_baseline/VSEFIE_SWG_TriTetra.csv")
    if !isfile(legacy_file)
        println("WARNING: Legacy baseline not found!")
        return
    end
    println("\nLoading Legacy baseline...")

    legacy_lines = readlines(legacy_file)
    legacy_theta = Float64[]
    legacy_phi = Float64[]
    legacy_rcs = Float64[]
    for line in legacy_lines[2:end]
        parts = split(line, ",")
        push!(legacy_theta, parse(Float64, parts[1]))
        push!(legacy_phi, parse(Float64, parts[2]))
        push!(legacy_rcs, parse(Float64, parts[3]))
    end

    # Compare E-plane
    mask_e = legacy_phi .== 0.0
    legacy_e_rcs = legacy_rcs[mask_e]
    emsuite_e_rcs = RCS_dB[:, 1]
    n_match = min(length(legacy_e_rcs), length(emsuite_e_rcs))
    diff = emsuite_e_rcs[1:n_match] .- legacy_e_rcs[1:n_match]

    println("\n==================================================")
    println("   Comparison: EMSuite vs Legacy (E-plane)        ")
    println("==================================================")
    println("  Points: $n_match")
    println("  Mean Diff: $(round(mean(diff), digits=4)) dB")
    rmse = sqrt(mean(diff .^ 2))
    println("  RMSE:      $(round(rmse, digits=4)) dB")
    println("  Max Diff:  $(round(maximum(abs.(diff)), digits=4)) dB")

    diff_dm = diff .- mean(diff)
    println("  RMSE (de-meaned): $(round(sqrt(mean(diff_dm.^2)), digits=4)) dB")

    println("\n=== Key Angle Comparison ===")
    @printf("  %5s  %12s  %12s  %10s\n", "θ(°)", "EMSuite", "Legacy", "Diff(dB)")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        if idx <= n_match
            @printf("  %5d  %12.4f  %12.4f  %10.4f\n",
                deg, emsuite_e_rcs[idx], legacy_e_rcs[idx],
                emsuite_e_rcs[idx] - legacy_e_rcs[idx])
        end
    end

    if rmse < 2.0
        println("\n✅ E1-VSEFIE PASS: RMSE = $(round(rmse, digits=4)) dB < 2.0 dB threshold")
    else
        println("\n❌ E1-VSEFIE FAIL: RMSE = $(round(rmse, digits=4)) dB ≥ 2.0 dB threshold")
    end
end

verify_vsefie_direct()
