"""
    verify_SCFIE_mlfma_sphere.jl

C3: S-CFIE MLFMA自洽性验证 (Sphere 600MHz)。
在 sphere_600MHz 网格上比较 CFIE Direct vs CFIE MLFMA。

用法: julia --project=. benchmark/verify_SCFIE_mlfma_sphere.jl
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

function verify_SCFIE_mlfma_sphere()
    println("=" ^ 72)
    println("  C3: S-CFIE MLFMA Self-Consistency — Sphere 600MHz")
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

    # 3. Basis and Equation
    set_frequency!(freq)
    basis_d = RWGBasis(mesh)
    basis_m = RWGBasis(mesh)
    N = num_basis(basis_d)
    println("  Unknowns: $N")

    cfie = CFIE(freq, alpha)
    println("  CFIE alpha=$alpha")

    # 4. Excitation
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

    # 5. Direct solve
    println("\n--- Direct ---")
    t1 = @elapsed Z = assemble_impedance_matrix(cfie, basis_d)
    println("  Assembly: $(round(t1, digits=2))s")

    V_d = excitation_vector(cfie, source, basis_d)
    t2 = @elapsed I_direct = solve!(LUSolver(), Z, V_d)
    println("  Solve: $(round(t2, digits=2))s")

    # 6. MLFMA solve
    println("\n--- MLFMA ---")
    leaf_size = 0.25 * lambda
    t3 = @elapsed mlfma_op = MLFMAOperator(cfie, basis_m, leaf_size)
    println("  MLFMA setup: $(round(t3, digits=2))s")

    V_m = excitation_vector(cfie, source, basis_m)
    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    t4 = @elapsed I_mlfma = solve!(solver, mlfma_op, V_m, Pl=P)
    println("  GMRES: $(round(t4, digits=2))s")

    # 7. Coefficient comparison
    coeff_err = norm(I_direct - I_mlfma) / norm(I_direct)
    println("\n--- Coefficient Comparison ---")
    println("  ||I_direct - I_mlfma|| / ||I_direct|| = $(round(coeff_err*100, digits=4))%")

    # 8. RCS comparison
    θs = collect(LinRange(-π, π, 721))
    ϕs = [0.0, π/2]

    RCS_d = radarCrossSection(θs, ϕs, I_direct, basis_d)
    RCS_m = radarCrossSection(θs, ϕs, I_mlfma, basis_m)

    dBsm_d = 10 * log10.(RCS_d[2])
    dBsm_m = 10 * log10.(RCS_m[2])

    diff = dBsm_d[:, 1] .- dBsm_m[:, 1]

    println("\n--- RCS Comparison (Phi=0) ---")
    println("  Mean |Diff|: $(round(mean(abs.(diff)), digits=4)) dB")
    println("  Max  |Diff|: $(round(maximum(abs.(diff)), digits=4)) dB")
    println("  RMSE:        $(round(sqrt(mean(diff.^2)), digits=4)) dB")

    println("\n--- 典型角度 (Phi=0) ---")
    for (i, th) in enumerate(θs)
        td = round(th * 180 / π, digits=1)
        if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
            @printf("  θ=%+7.1f° | Direct=%+8.2f | MLFMA=%+8.2f | Diff=%+6.3f dB\n",
                td, dBsm_d[i,1], dBsm_m[i,1], diff[i])
        end
    end

    # 9. Also compare MLFMA with Legacy baseline
    baseline_file = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SCFIE_Direct_Sphere.csv")
    if isfile(baseline_file)
        df_base = CSV.read(baseline_file, DataFrame)
        RCS_legacy_dB = df_base.RCS_Phi0_dBsm
        diff_legacy = dBsm_m[:, 1] .- RCS_legacy_dB
        println("\n--- MLFMA vs Legacy Direct ---")
        println("  Mean Diff:   $(round(mean(diff_legacy), digits=4)) dB")
        println("  RMSE:        $(round(sqrt(mean(diff_legacy.^2)), digits=4)) dB")
        println("  Max |Diff|:  $(round(maximum(abs.(diff_legacy)), digits=4)) dB")
    else
        println("\n  Legacy baseline not found, skipping comparison.")
    end

    println("\n--- Timing Summary ---")
    println("  Direct: Assembly=$(round(t1,digits=1))s + Solve=$(round(t2,digits=1))s = $(round(t1+t2,digits=1))s")
    println("  MLFMA:  Setup=$(round(t3,digits=1))s + GMRES=$(round(t4,digits=1))s = $(round(t3+t4,digits=1))s")

    # Pass/Fail
    if coeff_err < 0.01  # 1%
        println("\n✅ C3 PASS: SCFIE MLFMA coefficient error $(round(coeff_err*100, digits=4))% < 1%")
    else
        println("\n❌ C3 FAIL: SCFIE MLFMA coefficient error $(round(coeff_err*100, digits=4))% >= 1%")
    end
end

verify_SCFIE_mlfma_sphere()
