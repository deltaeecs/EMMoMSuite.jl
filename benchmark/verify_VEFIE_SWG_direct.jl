# D1: VEFIE Direct — EMMoMSuite SWG vs Legacy SWG on Tetra.nas
# Mesh: Tetra.nas (pure dielectric body), 300 MHz, εr=2.0
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using EMMoMSuite.CoreModule.Sources
using LinearAlgebra
using Printf
using DelimitedFiles

function verify_vefie_swg()
    println("==================================================")
    println("   D1: VEFIE Direct — EMMoMSuite SWG vs Legacy       ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    c0 = 299792458.0
    lambda = c0 / freq
    EMMoMSuite.Utilities.Parameters.set_frequency!(freq)
    println("f = $(freq/1e6) MHz, λ = $(round(lambda, digits=4)) m")

    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../deps/fixtures/AllinOne/meshfiles/Tetra.nas")
    @assert isfile(mesh_file) "Mesh not found: $mesh_file"
    println("Loading mesh: $mesh_file")
    mesh = read_nas_mesh(mesh_file)

    # Scale mm → m (same as Legacy meshUnit=:mm)
    nodes = mesh.node
    max_coord = maximum(abs.(nodes))
    if max_coord > 10.0
        println("Scaling mesh by 0.001 (mm → m)...")
        mesh.node .*= 0.001
    end

    # Center the mesh
    min_c = minimum(mesh.node, dims=2)
    max_c = maximum(mesh.node, dims=2)
    center = (min_c + max_c) / 2
    mesh.node .-= center

    println("  $(num_vertices(mesh)) vertices, $(num_elements(mesh)) tetrahedra")

    # 3. SWG Basis
    println("Setting up SWG Basis...")
    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("  N = $N unknowns")

    # 4. Permittivity
    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))

    # 5. VEFIE Operator
    println("Setting up VEFIE...")
    vefie = VEFIE(freq, permittivities)

    # 6. Assembly
    println("Assembling Z matrix...")
    t_asm = @elapsed Z = assemble_impedance_matrix(vefie, basis, permittivities)
    println("  Done in $(round(t_asm, digits=2)) s, size=$(size(Z))")
    println("  Z[1,1] = $(Z[1,1])")
    println("  ||Z|| = $(norm(Z))")

    # 7. Excitation (plane wave from +z, x-pol)
    println("Computing excitation...")
    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(vefie, source, basis, permittivities)
    println("  ||V|| = $(norm(V))")

    # 8. Solve
    println("Solving (LU)...")
    solver = LUSolver()
    t_solve = @elapsed D = solve!(solver, Z, V)
    println("  Done in $(round(t_solve, digits=2)) s, ||D|| = $(norm(D))")

    # 9. RCS
    println("Computing RCS...")
    θs = collect(0:1.0:180.0) .* (π / 180.0)
    ϕs = [0.0, π / 2]
    RCS_components, RCS_total, RCS_dB = radarCrossSection(θs, ϕs, D, basis, permittivities)

    # Save EMMoMSuite result
    outdir = joinpath(@__DIR__, "../test_results/emsuite_verification")
    mkpath(outdir)
    outfile = joinpath(outdir, "VEFIE_SWG_Tetra.csv")
    open(outfile, "w") do io
        println(io, "theta_deg,phi_deg,RCS_dBsm")
        for p in 1:length(ϕs)
            for t in 1:length(θs)
                @printf(io, "%.4f,%.4f,%.6f\n",
                    rad2deg(θs[t]), rad2deg(ϕs[p]), RCS_dB[t, p])
            end
        end
    end
    println("  EMMoMSuite RCS saved to $outfile")

    # 10. Load Legacy baseline
    legacy_file = joinpath(@__DIR__, "../test_results/legacy_baseline/VEFIE_SWG_Tetra.csv")
    if !isfile(legacy_file)
        println("WARNING: Legacy baseline not found: $legacy_file")
        return
    end
    println("\nLoading Legacy baseline: $legacy_file")

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

    # Compare E-plane (phi=0) only
    mask_e = legacy_phi .== 0.0
    legacy_e_theta = legacy_theta[mask_e]
    legacy_e_rcs = legacy_rcs[mask_e]
    emsuite_e_rcs = RCS_dB[:, 1]

    n_match = min(length(legacy_e_rcs), length(emsuite_e_rcs))
    diff = emsuite_e_rcs[1:n_match] .- legacy_e_rcs[1:n_match]

    println("\n==================================================")
    println("   Comparison: EMMoMSuite vs Legacy (E-plane, phi=0) ")
    println("==================================================")
    println("  Points: $n_match")
    println("  Mean Diff: $(round(mean(diff), digits=4)) dB")
    println("  RMSE:      $(round(sqrt(mean(diff.^2)), digits=4)) dB")
    println("  Max Diff:  $(round(maximum(abs.(diff)), digits=4)) dB")

    # De-meaned RMSE
    mean_offset = mean(diff)
    diff_demeaned = diff .- mean_offset
    rmse_dm = sqrt(mean(diff_demeaned .^ 2))
    println("  RMSE (de-meaned): $(round(rmse_dm, digits=4)) dB")

    println("\n=== Key Angle Comparison ===")
    @printf("  %5s  %12s  %12s  %10s\n", "θ(°)", "EMMoMSuite", "Legacy", "Diff(dB)")
    for deg in [0, 30, 60, 90, 120, 150, 180]
        idx = deg + 1
        if idx <= n_match
            @printf("  %5d  %12.4f  %12.4f  %10.4f\n",
                deg, emsuite_e_rcs[idx], legacy_e_rcs[idx],
                emsuite_e_rcs[idx] - legacy_e_rcs[idx])
        end
    end

    # H-plane comparison
    mask_h = legacy_phi .== 90.0
    if any(mask_h)
        legacy_h_rcs = legacy_rcs[mask_h]
        emsuite_h_rcs = RCS_dB[:, 2]
        n_h = min(length(legacy_h_rcs), length(emsuite_h_rcs))
        diff_h = emsuite_h_rcs[1:n_h] .- legacy_h_rcs[1:n_h]
        println("\n=== H-plane (phi=90°) ===")
        println("  RMSE: $(round(sqrt(mean(diff_h.^2)), digits=4)) dB")
        println("  Max Diff: $(round(maximum(abs.(diff_h)), digits=4)) dB")
    end

    # PASS/FAIL
    rmse = sqrt(mean(diff .^ 2))
    if rmse < 1.0
        println("\n✅ D1-SWG PASS: RMSE = $(round(rmse, digits=4)) dB < 1.0 dB threshold")
    else
        println("\n❌ D1-SWG FAIL: RMSE = $(round(rmse, digits=4)) dB ≥ 1.0 dB threshold")
    end
end

using Statistics: mean

verify_vefie_swg()
