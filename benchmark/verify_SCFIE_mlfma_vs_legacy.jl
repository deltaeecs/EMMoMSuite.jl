"""
    verify_SCFIE_mlfma_vs_legacy.jl

C3: S-CFIE MLFMA vs Legacy Direct 基线对比。
在 sphere_600MHz 网格上运行 CFIE MLFMA，与 Legacy SCFIE_Direct_Sphere.csv 对比。
直接跳过 Direct 求解（26k 未知数太慢）。

用法: julia --project=. benchmark/verify_SCFIE_mlfma_vs_legacy.jl
"""

using EMMoMSuite
using LinearAlgebra
using Printf
using Statistics
using CSV, DataFrames

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

function verify_SCFIE_mlfma_vs_legacy()
    println("=" ^ 72)
    println("  C3: S-CFIE MLFMA vs Legacy Direct — Sphere 600MHz")
    println("=" ^ 72)

    # 1. Parameters
    freq = 6e8      # 600 MHz
    alpha = 0.6     # Match Legacy default
    lambda = 299792458.0 / freq

    # 2. Mesh
    mesh_dir = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne", "meshfiles")
    mesh_file = joinpath(mesh_dir, "sphere_600MHz.nas")
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file")
    end
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    set_frequency!(freq)

    # 3. Basis and Equation
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("  Unknowns: $N")

    cfie = CFIE(freq, alpha)
    println("  CFIE alpha=$alpha")

    # 4. Excitation — match Legacy: Direction -x, Pol +z
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

    # 5. MLFMA solve
    println("\n--- MLFMA ---")
    leaf_size = 0.25 * lambda
    t1 = @elapsed mlfma_op = MLFMAOperator(cfie, basis, leaf_size)
    println("  MLFMA setup: $(round(t1, digits=2))s")

    V = excitation_vector(cfie, source, basis)
    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    t2 = @elapsed I_mlfma = solve!(solver, mlfma_op, V, Pl=P)
    println("  GMRES: $(round(t2, digits=2))s")

    # 6. RCS — 721 points matches Legacy baseline 
    θs = collect(LinRange(-π, π, 721))
    ϕs = [0.0, π/2]

    RCS_res = radarCrossSection(θs, ϕs, I_mlfma, basis)
    dBsm = 10 * log10.(RCS_res[2])

    # 7. Compare with Legacy
    baseline_file = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SCFIE_Direct_Sphere.csv")
    if isfile(baseline_file)
        df_base = CSV.read(baseline_file, DataFrame)
        RCS_legacy_dB = df_base.RCS_Phi0_dBsm

        diff = dBsm[:, 1] .- RCS_legacy_dB

        println("\n--- MLFMA vs Legacy Direct (Phi=0) ---")
        println("  Mean Diff:   $(round(mean(diff), digits=4)) dB")
        println("  Mean |Diff|: $(round(mean(abs.(diff)), digits=4)) dB")
        println("  RMSE:        $(round(sqrt(mean(diff.^2)), digits=4)) dB")
        println("  Max |Diff|:  $(round(maximum(abs.(diff)), digits=4)) dB")

        println("\n--- 典型角度 (Phi=0) ---")
        for (i, th) in enumerate(θs)
            td = round(th * 180 / π, digits=1)
            if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
                @printf("  θ=%+7.1f° | Legacy=%+8.2f | MLFMA=%+8.2f | Diff=%+6.3f dB\n",
                    td, RCS_legacy_dB[i], dBsm[i,1], diff[i])
            end
        end

        # Save
        output_dir = joinpath(@__DIR__, "..", "test_results", "emsuite_verification")
        mkpath(output_dir)
        df_out = DataFrame(
            Theta_Rad = θs,           
            RCS_MLFMA_Phi0_dB = dBsm[:, 1],
            RCS_Legacy_Phi0_dB = RCS_legacy_dB,
            Diff_dB = diff
        )
        CSV.write(joinpath(output_dir, "SCFIE_MLFMA_Sphere_vs_Legacy.csv"), df_out)
        println("\n  Saved to SCFIE_MLFMA_Sphere_vs_Legacy.csv")

        # Pass/Fail
        rmse = sqrt(mean(diff.^2))
        if rmse < 1.0
            println("\n✅ C3 PASS: SCFIE MLFMA vs Legacy RMSE = $(round(rmse, digits=4)) dB < 1.0 dB")
        else
            println("\n❌ C3 FAIL: SCFIE MLFMA vs Legacy RMSE = $(round(rmse, digits=4)) dB >= 1.0 dB")
        end
    else
        println("\n❌ Legacy baseline not found: $baseline_file")
    end

    println("\n--- Timing ---")
    println("  MLFMA total: $(round(t1+t2, digits=1))s (setup=$(round(t1,digits=1))s + GMRES=$(round(t2,digits=1))s)")
end

verify_SCFIE_mlfma_vs_legacy()
