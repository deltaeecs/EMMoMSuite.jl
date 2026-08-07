"""
run_P1_P3_pmchw.jl — Phase 14 精度基准：PMCHW 介质球散射 P1-P3

测试说明：
  P1  PMCHW Direct (εᵣ=4, 无损, r=0.15m, 600 MHz)     vs Mie 介质级数 (RMSE < 2 dB)
  P3  PMCHW Direct (εᵣ=2.2-0.1j, 有损, r=0.15m, 300MHz) vs Mie 介质级数 (RMSE < 2 dB)

注: P2 (PMCHW MLFMA) 需先实现 PMCHWMLFMAOperator，当前跳过。

几何: dielectric sphere meshed with RWG surface triangles
网格: sphere_600MHz.nas  → r ≈ 0.15m（或 sphere, 重用球面 RWG 网格）

运行方式:
    julia -t auto --project=. benchmark/accuracy/run_P1_P3_pmchw.jl [P1] [P3]
  (不带参数时运行 P1)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using EMMoMSuite
using EMMoMSuite.Accuracy.ReferenceData: extract_sphere_radius, mie_dielectric_bistatic_rcs_dBsm
using LinearAlgebra, Printf, Statistics, Dates
using Base.Threads

function dense_complex_gib(dim)
    return dim * dim * sizeof(ComplexF64) / 2.0^30
end

# ─── 路径 ─────────────────────────────────────────────────────────────────────
const ROOT_DIR   = joinpath(@__DIR__, "..", "..")
const MESH_DIR   = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "meshfiles")
const RESULT_DIR = joinpath(ROOT_DIR, "test_results", "accuracy")
mkpath(RESULT_DIR)

enabled = isempty(ARGS) ? ["P1"] : ARGS

const θs_obs = collect(range(-π, π, length = 721))
const ϕs_obs = [0.0, π/2]
const θ_deg_obs = rad2deg.(θs_obs)

function report_runtime_configuration(case_label)
    julia_threads = nthreads()
    blas_threads = BLAS.get_num_threads()
    println("  Julia threads: $julia_threads")
    println("  BLAS threads: $blas_threads")
    if julia_threads == 1
        @warn "$case_label 当前以单线程运行，PMCHW 稠密组装会显著变慢。建议使用: julia -t auto --project=. benchmark/accuracy/run_P1_P3_pmchw.jl $case_label"
    end
    if blas_threads == 1
        @warn "$case_label 当前 BLAS 仅 1 线程，稠密 LU 求解时间会偏高。"
    end
end

function save_and_report_case(label, θdeg, rcs_ems, rcs_ref, threshold)
    result = compute_rcs_accuracy(rcs_ems, rcs_ref, θdeg, label; threshold)
    @printf("  %-48s  RMSE=%6.3f dB  %s\n", label, result.rmse_dB,
            result.pass ? "✓ PASS" : "✗ FAIL (门限 $(threshold) dB)")
    csv_path = joinpath(RESULT_DIR, label * ".csv")
    open(csv_path, "w") do io
        println(io, "theta_deg,rcs_ems_dBsm,rcs_mie_dBsm,diff_dB")
        for i in eachindex(θdeg)
            @printf(io, "%.2f,%.6f,%.6f,%.6f\n",
                    θdeg[i], rcs_ems[i], rcs_ref[i], rcs_ems[i] - rcs_ref[i])
        end
    end
    return result
end

# ─── 共用网格（球面三角网格用于 PMCHW） ──────────────────────────────────────
mesh_file = joinpath(MESH_DIR, "sphere_600MHz.nas")
isfile(mesh_file) || error("网格文件不存在: $mesh_file")

all_results = AccuracyResult[]

# ─── P1: PMCHW Direct, εᵣ=4, 无损, 600 MHz ───────────────────────────────────
if "P1" in enabled
    println("\n[P1] PMCHW Direct (εᵣ=4, 无损, 600 MHz)")
    report_runtime_configuration("P1")
    freq  = 6.0e8
    eps_r = 4.0 + 0im     # 无损介质
    mu_r  = 1.0 + 0im

    set_frequency!(freq)
    mesh  = read_nas_mesh(mesh_file, scale = 1.0)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)
    println("  RWG 未知量: $N (2N = $(2N))")
    @printf("  稠密矩阵理论占用: %.2f GiB\n", dense_complex_gib(2N))

    pmchw  = PMCHW(freq, eps_r, mu_r)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

    t1 = @elapsed Z = assemble_impedance_matrix(pmchw, basis)
    @printf("  PMCHW 组装: %.1fs  (%d×%d)\n", t1, size(Z)...)
    V = excitation_vector(pmchw, source, basis)
    t2 = @elapsed begin
        F = lu!(Z)
        I = similar(V)
        ldiv!(I, F, V)
    end
    @printf("  LU 求解: %.1fs\n", t2)
    F = nothing; Z = nothing; GC.gc()

    k0_val   = Float64(pmchw.k0)
    eta0_val = Float64(pmchw.eta0)
    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis, k0_val, eta0_val)
    rcs_phi0  = rcs_dB[:, 1]
    rcs_phi90 = rcs_dB[:, 2]

    # Mie 参考解
    radius_m  = extract_sphere_radius(mesh_file)
    @printf("  球半径: %.4f m\n", radius_m)
    mie_phi0  = mie_dielectric_bistatic_rcs_dBsm(radius_m, freq, eps_r, mu_r, θs_obs, 0.0, π/2, π, [0.0, 0.0, 1.0])
    mie_phi90 = mie_dielectric_bistatic_rcs_dBsm(radius_m, freq, eps_r, mu_r, θs_obs, π/2, π/2, π, [0.0, 0.0, 1.0])

    r1 = save_and_report_case("P1_PMCHW_Sphere_Direct_phi0_vs_Mie",
        θ_deg_obs, rcs_phi0, mie_phi0, 2.0)
    r2 = save_and_report_case("P1_PMCHW_Sphere_Direct_phi90_vs_Mie",
        θ_deg_obs, rcs_phi90, mie_phi90, 2.0)
    push!(all_results, r1, r2)
end

# ─── P2 跳过说明 ─────────────────────────────────────────────────────────────
if "P2" in enabled
    @warn "P2 (PMCHW MLFMA) 已跳过。原因: PMCHWMLFMAOperator 尚未实现。" *
          "待完成实现后在 Phase 14 后续子任务中运行。"
end

# ─── P3: PMCHW Direct, εᵣ=2.2-0.1j, 有损, 300 MHz ───────────────────────────
if "P3" in enabled
    println("\n[P3] PMCHW Direct (εᵣ=2.2-0.1j, 有损, 300 MHz)")
    report_runtime_configuration("P3")
    freq  = 3.0e8
    eps_r = 2.2 - 0.1im   # 有损介质
    mu_r  = 1.0 + 0im

    set_frequency!(freq)
    mesh  = read_nas_mesh(mesh_file, scale = 1.0)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)
    println("  RWG 未知量: $N (2N = $(2N))")
    @printf("  稠密矩阵理论占用: %.2f GiB\n", dense_complex_gib(2N))

    pmchw  = PMCHW(freq, eps_r, mu_r)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

    t1 = @elapsed Z = assemble_impedance_matrix(pmchw, basis)
    @printf("  PMCHW 组装: %.1fs  (%d×%d)\n", t1, size(Z)...)
    V = excitation_vector(pmchw, source, basis)
    t2 = @elapsed begin
        F = lu!(Z)
        I = similar(V)
        ldiv!(I, F, V)
    end
    @printf("  LU 求解: %.1fs\n", t2)
    F = nothing; Z = nothing; GC.gc()

    k0_val  = Float64(pmchw.k0)
    eta0_val = Float64(pmchw.eta0)
    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis, k0_val, eta0_val)
    rcs_phi0  = rcs_dB[:, 1]
    rcs_phi90 = rcs_dB[:, 2]

    # Mie 参考解（有损介质）
    radius_m = extract_sphere_radius(mesh_file)
    @printf("  球半径: %.4f m\n", radius_m)
    mie_phi0  = mie_dielectric_bistatic_rcs_dBsm(radius_m, freq, eps_r, mu_r, θs_obs, 0.0, π/2, π, [0.0, 0.0, 1.0])
    mie_phi90 = mie_dielectric_bistatic_rcs_dBsm(radius_m, freq, eps_r, mu_r, θs_obs, π/2, π/2, π, [0.0, 0.0, 1.0])

    r1 = save_and_report_case("P3_PMCHW_LossySphere_Direct_phi0_vs_Mie",
        θ_deg_obs, rcs_phi0, mie_phi0, 2.0)
    r2 = save_and_report_case("P3_PMCHW_LossySphere_Direct_phi90_vs_Mie",
        θ_deg_obs, rcs_phi90, mie_phi90, 2.0)
    push!(all_results, r1, r2)
end

# ─── 汇总 ────────────────────────────────────────────────────────────────────
if !isempty(all_results)
    println()
    print_accuracy_report(all_results; title = "P1/P3 PMCHW 介质球精度报告")

    open(joinpath(RESULT_DIR, "P1_P3_pmchw_report.md"), "w") do io
        println(io, "# P1/P3 PMCHW 介质球散射精度报告 (vs Mie 介质级数)")
        println(io, "\n生成时间: $(Dates.now())\n")
        println(io, "> P2 (PMCHW MLFMA) 待 PMCHWMLFMAOperator 实现后运行。\n")
        print_accuracy_report(io, all_results; title = "P1/P3 PMCHW")
    end
end
