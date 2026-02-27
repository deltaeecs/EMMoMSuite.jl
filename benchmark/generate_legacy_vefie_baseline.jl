# Generate Legacy VEFIE baseline using MoM_AllinOne.
# Mesh: Tetra.nas (pure dielectric body), 300 MHz, εr=2.0
# Basis: SWG on tetrahedra
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../LegacyBenchmark"))

using MoM_Basics
using MoM_Kernels
using LinearAlgebra
using Printf

function generate_legacy_vefie()
    println("==================================================")
    println("   Legacy VEFIE Baseline: SWG on Tetra            ")
    println("==================================================")

    setPrecision!(Float64)
    SimulationParams.SHOWIMAGE = false

    # Parameters — matches EMSuite verify_VEFIE_direct.jl
    filename = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/Tetra.nas")
    @assert isfile(filename) "Mesh file not found: $filename"
    meshUnit = :mm
    frequency = 300e6
    ieT = :EFIE
    sbfT = :nothing
    vbfT = :SWG
    solverT = :direct

    # Source: plane wave from +z, x-polarized
    # Legacy PlaneWave(theta, phi, alpha, V): theta=π → incidence from +z
    source = PlaneWave(Float64(π), 0.0, 0.0, 1.0)

    # 1. Update Parameters
    inputParameters(; frequency=frequency, ieT=ieT)
    updateVSBFTParams!(; sbfT=sbfT, vbfT=vbfT)

    # 2. Read Mesh
    println("Reading mesh: $filename (unit=$meshUnit)...")
    meshData, εᵣs = getMeshData(filename; meshUnit=meshUnit)

    # 3. Basis Functions
    println("Generating SWG basis...")
    ngeo, nbf, geosInfo, bfsInfo = getBFsFromMeshData(meshData; sbfT=sbfT, vbfT=vbfT)
    println("  ngeo=$ngeo, nbf=$nbf")

    # 4. Permittivity
    println("Setting εr=2.0...")
    setGeosPermittivity!(geosInfo, 2.0 + 0.0im)

    # 5. Impedance Matrix
    println("Assembling Z matrix...")
    t_asm = @elapsed Zmat = getImpedanceMatrix(geosInfo, nbf)
    println("  Assembly: $(round(t_asm, digits=2)) s")
    println("  Z size: $(size(Zmat))")
    println("  Z[1,1] = $(Zmat[1,1])")
    println("  ||Z|| = $(norm(Zmat))")

    # 6. Excitation
    println("Computing excitation vector...")
    V = getExcitationVector(geosInfo, size(Zmat, 1), source)
    println("  ||V|| = $(norm(V))")

    # 7. Solve
    println("Solving with direct solver...")
    t_solve = @elapsed begin
        ICoeff, _ = solve(Zmat, V; solverT=solverT)
    end
    println("  Solve: $(round(t_solve, digits=2)) s")
    println("  ||I|| = $(norm(ICoeff))")

    # 8. RCS
    println("Computing RCS...")
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]  # E-plane and H-plane

    t_rcs = @elapsed begin
        RCSθsϕs, RCSθsϕsdB, RCS, RCSdB = radarCrossSection(θs, ϕs, ICoeff, geosInfo)
    end
    println("  RCS: $(round(t_rcs, digits=2)) s")

    # 9. Save
    outdir = joinpath(@__DIR__, "../test_results/legacy_baseline")
    mkpath(outdir)
    outfile = joinpath(outdir, "VEFIE_SWG_Tetra.csv")
    open(outfile, "w") do io
        println(io, "theta_deg,phi_deg,RCS_dBsm")
        for p in 1:length(ϕs)
            for t in 1:length(θs)
                @printf(io, "%.4f,%.4f,%.6f\n",
                    rad2deg(θs[t]), rad2deg(ϕs[p]), RCSdB[t, p])
            end
        end
    end
    println("Saved Legacy baseline to $outfile")
    println("  Lines: $(length(θs) * length(ϕs))")

    # Print key values for comparison
    println("\n=== Key RCS Values (E-plane, phi=0) ===")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        @printf("  θ=%3d°: %10.4f dBsm\n", deg, RCSdB[idx, 1])
    end

    return Zmat, ICoeff, RCSdB
end

generate_legacy_vefie()
