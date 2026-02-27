# Diagnostic Legacy VSEFIE: get Z block norms, solution norms, and per-component RCS
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../LegacyBenchmark"))

using MoM_Basics
using MoM_Kernels
using LinearAlgebra
using Printf

function diagnose_legacy_vsefie()
    println("==================================================")
    println("   Legacy VSEFIE Diagnostics: TriTetra             ")
    println("==================================================")

    setPrecision!(Float64)
    SimulationParams.SHOWIMAGE = false

    filename = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
    meshUnit = :mm
    frequency = 2e9
    ieT = :EFIE
    sbfT = :RWG
    vbfT = :SWG
    solverT = :direct

    source = PlaneWave(Float64(π), 0.0, 0.0, 1.0)

    inputParameters(; frequency=frequency, ieT=ieT)
    updateVSBFTParams!(; sbfT=sbfT, vbfT=vbfT)

    meshData, εᵣs = getMeshData(filename; meshUnit=meshUnit)
    ngeo, nbf, geosInfo, bfsInfo = getBFsFromMeshData(meshData; sbfT=sbfT, vbfT=vbfT)

    eps_r = 2.0 * (1 - 0.001im)
    setGeosPermittivity!(geosInfo, eps_r)

    # Get counts
    tris = geosInfo[1]
    geosV = geosInfo[2]
    ntri = length(tris)
    ntet = length(geosV)

    # Count basis functions per type
    n_surf = 0
    for tri in tris
        for m in tri.inBfsID
            m > n_surf && (n_surf = m)
        end
    end
    n_total = nbf
    n_vol = n_total - n_surf
    println("  n_surf=$n_surf, n_vol=$n_vol, n_total=$n_total")

    # Assembly
    println("Assembling Z matrix...")
    Zmat = getImpedanceMatrix(geosInfo, nbf)
    println("  ||Z|| = $(norm(Zmat))")

    # Z block norms
    Z_SS = Zmat[1:n_surf, 1:n_surf]
    Z_SV = Zmat[1:n_surf, n_surf+1:end]
    Z_VS = Zmat[n_surf+1:end, 1:n_surf]
    Z_VV = Zmat[n_surf+1:end, n_surf+1:end]
    println("\n=== Z Block Norms ===")
    println("  ||Z_SS|| = $(norm(Z_SS))")
    println("  ||Z_SV|| = $(norm(Z_SV))")
    println("  ||Z_VS|| = $(norm(Z_VS))")
    println("  ||Z_VV|| = $(norm(Z_VV))")
    println("  ||Z_SV||/||Z_VS|| = $(norm(Z_SV)/norm(Z_VS))")

    # Excitation
    V = getExcitationVector(geosInfo, size(Zmat, 1), source)
    V_s = V[1:n_surf]
    V_v = V[n_surf+1:end]
    println("  ||V_surf|| = $(norm(V_s))")
    println("  ||V_vol||  = $(norm(V_v))")

    # Solve
    ICoeff, _ = solve(Zmat, V; solverT=solverT)
    I_surf = ICoeff[1:n_surf]
    D_vol = ICoeff[n_surf+1:end]
    println("  ||I_surf|| = $(norm(I_surf))")
    println("  ||D_vol||  = $(norm(D_vol))")

    # Incident direction
    println("\n=== Incident Direction ===")
    println("  k̂ = $(source.k̂)")
    println("  E_v = $(source.E_v)")

    # RCS
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]

    # Combined RCS
    RCSθsϕs, RCSθsϕsdB, RCS, RCSdB = radarCrossSection(θs, ϕs, ICoeff, geosInfo)
    
    # Surface-only RCS
    println("\nComputing surface-only RCS...")
    RCS_s_θsϕs, RCS_s_dB_θsϕs, RCS_s, RCS_s_dB = radarCrossSection(θs, ϕs, ICoeff, [geosInfo[1]])
    
    # Volume-only RCS: need to zero out surface coefficients
    # Actually, we need to be careful. ICoeff contains both surface and volume unknowns.
    # For surface-only: use ICoeff[1:n_surf] with just surface triangles
    # For volume-only: use ICoeff[n_surf+1:end] with just volume tets
    # But radarCrossSection expects full ICoeff... 
    # Let me try with zeroed coefficients
    I_surf_only = copy(ICoeff)
    I_surf_only[n_surf+1:end] .= 0
    
    I_vol_only = copy(ICoeff) 
    I_vol_only[1:n_surf] .= 0

    RCS_so_θsϕs, _, _, RCS_so_dB = radarCrossSection(θs, ϕs, I_surf_only, geosInfo)
    RCS_vo_θsϕs, _, _, RCS_vo_dB = radarCrossSection(θs, ϕs, I_vol_only, geosInfo)

    println("\n=== Surface vs Volume RCS at key angles (E-plane) ===")
    @printf("  %5s  %12s  %12s  %12s\n", "θ(°)", "Surf-only", "Vol-only", "Combined")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        @printf("  %5d  %12.4f  %12.4f  %12.4f\n",
            deg, RCS_so_dB[idx, 1], RCS_vo_dB[idx, 1], RCSdB[idx, 1])
    end
end

diagnose_legacy_vsefie()
