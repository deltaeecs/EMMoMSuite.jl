"""
run_F5_F6_sphere.jl — Phase 14 精度基准：Sphere 散射 F5-F6

测试说明：
  F5  S-CFIE Direct  vs Feko + Mie PEC  (RMSE < 2 dB)
  F6  S-CFIE MLFMA   vs Feko + Mie PEC  (RMSE < 3 dB)
  X1  S-EFIE Direct  vs Mie PEC         (RMSE < 2 dB)  [可选]

几何: PEC 球, r ≈ 0.3m, 频率 = 600 MHz
基线: Feko (sphere_600MHzRCS.csv) + Mie 解析解

运行方式:
    julia -t auto --project=. benchmark/accuracy/run_F5_F6_sphere.jl [F5] [F6] [X1]
  (不带参数时运行 F5)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using EMMoMSuite
using LinearAlgebra, Printf, Statistics, Dates
using Base.Threads

struct LUPreconditioner; F; end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

# ─── 路径 ─────────────────────────────────────────────────────────────────────
const ROOT_DIR   = joinpath(@__DIR__, "..", "..")
const MESH_DIR   = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "meshfiles")
const FEKO_DIR   = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "deps", "compare_feko")
const RESULT_DIR = joinpath(ROOT_DIR, "test_results", "accuracy")
mkpath(RESULT_DIR)

enabled = isempty(ARGS) ? ["F5"] : ARGS

const θs_obs = collect(range(-π, π, length = 721))
const ϕs_obs = [0.0, π/2]

function report_runtime_configuration(case_label; require_threads = false)
    julia_threads = nthreads()
    blas_threads = BLAS.get_num_threads()
    println("  Julia threads: $julia_threads")
    println("  BLAS threads: $blas_threads")
    if require_threads && julia_threads == 1
        @warn "$case_label 当前以单线程运行，MLFMA 构建时间会明显偏高。建议使用: julia -t auto --project=. benchmark/accuracy/run_F5_F6_sphere.jl $case_label"
    end
end

function load_feko_phi_cuts(feko_file)
    isfile(feko_file) || error("Feko 文件不存在: $feko_file")
    theta_deg, phi_deg, _, rcs_dBsm = read_feko_rcs(feko_file)
    return split_phi_cuts(theta_deg, phi_deg, rcs_dBsm)
end

function save_and_report_case(label, θdeg, rcs_ems, rcs_ref, threshold)
    result = compute_rcs_accuracy(rcs_ems, rcs_ref, θdeg, label; threshold)
    @printf("  %-44s  RMSE=%6.3f dB  %s\n", label, result.rmse_dB,
            result.pass ? "✓ PASS" : "✗ FAIL (门限 $(threshold) dB)")
    csv_path = joinpath(RESULT_DIR, label * ".csv")
    open(csv_path, "w") do io
        println(io, "theta_deg,rcs_ems_dBsm,rcs_ref_dBsm,diff_dB")
        for i in eachindex(θdeg)
            @printf(io, "%.2f,%.6f,%.6f,%.6f\n",
                    θdeg[i], rcs_ems[i], rcs_ref[i], rcs_ems[i] - rcs_ref[i])
        end
    end
    return result
end

# ─── 共用设置 ─────────────────────────────────────────────────────────────────
freq = 6.0e8
λ    = 299792458.0 / freq
set_frequency!(freq)

mesh_file = joinpath(MESH_DIR, "sphere_600MHz.nas")
isfile(mesh_file) || error("网格不存在: $mesh_file")
mesh  = read_nas_mesh(mesh_file, scale = 1.0)
basis = RWGBasis(mesh)
N     = num_basis(basis)
println("Sphere 600MHz: $N 未知量")

source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

feko_file = joinpath(FEKO_DIR, "sphere_600MHzRCS.csv")
feko_cuts = load_feko_phi_cuts(feko_file)
θ_deg_feko = feko_cuts[0.0].theta    # 721 点 θ 向量（度）

# Mie 解析解（将全局观测角映射为相对入射方向的双站散射角）
sphere_nas = joinpath(MESH_DIR, "sphere_600MHz.nas")
mie_phi0 = mie_pec_bistatic_rcs_dBsm(
    sphere_nas, freq, θs_obs, 0.0, source.theta, source.phi, source.polarization)
mie_phi90 = mie_pec_bistatic_rcs_dBsm(
    sphere_nas, freq, θs_obs, π / 2, source.theta, source.phi, source.polarization)

all_results = AccuracyResult[]

# ─── X1: S-EFIE Direct vs Mie PEC ────────────────────────────────────────────
if "X1" in enabled
    println("\n[X1] S-EFIE Direct vs Mie PEC")
    efie = EFIE(freq)
    t1 = @elapsed Z = assemble_impedance_matrix(efie, basis)
    @printf("  组装: %.1fs\n", t1)
    V = excitation_vector(efie, source, basis)
    t2 = @elapsed I = Z \ V
    @printf("  LU: %.1fs\n", t2)
    Z = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    r1 = save_and_report_case("X1_SEFIE_Sphere_Direct_phi0_vs_Mie",
        rad2deg.(θs_obs), rcs_dB[:, 1], mie_phi0, 2.0)
    r2 = save_and_report_case("X1_SEFIE_Sphere_Direct_phi90_vs_Mie",
        rad2deg.(θs_obs), rcs_dB[:, 2], mie_phi90, 2.0)
    push!(all_results, r1, r2)
end

# ─── F5: S-CFIE Direct vs Feko + Mie ─────────────────────────────────────────
if "F5" in enabled
    println("\n[F5] S-CFIE Direct vs Feko + Mie")
    cfie = CFIE(freq, 0.5)
    t1 = @elapsed Z = assemble_impedance_matrix(cfie, basis)
    @printf("  组装: %.1fs\n", t1)
    V = excitation_vector(cfie, source, basis)
    t2 = @elapsed I = Z \ V
    @printf("  LU: %.1fs\n", t2)
    Z = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    rcs_phi0  = rcs_dB[:, 1]
    rcs_phi90 = rcs_dB[:, 2]
    θ_rad_deg = rad2deg.(θs_obs)

    # vs Feko
    r1f = save_and_report_case("F5_CFIE_Sphere_Direct_phi0_vs_Feko",
        θ_deg_feko, rcs_phi0, feko_cuts[0.0].rcs_dBsm, 2.0)
    r2f = save_and_report_case("F5_CFIE_Sphere_Direct_phi90_vs_Feko",
        θ_deg_feko, rcs_phi90, feko_cuts[90.0].rcs_dBsm, 2.0)
    # vs Mie
    r1m = save_and_report_case("F5_CFIE_Sphere_Direct_phi0_vs_Mie",
        θ_rad_deg, rcs_phi0, mie_phi0, 2.0)
    r2m = save_and_report_case("F5_CFIE_Sphere_Direct_phi90_vs_Mie",
        θ_rad_deg, rcs_phi90, mie_phi90, 2.0)
    push!(all_results, r1f, r2f, r1m, r2m)
end

# ─── F6: S-CFIE MLFMA vs Feko + Mie ─────────────────────────────────────────
if "F6" in enabled
    println("\n[F6] S-CFIE MLFMA+GMRES vs Feko + Mie")
    report_runtime_configuration("F6"; require_threads = true)
    cfie = CFIE(freq, 0.5)
    V = excitation_vector(cfie, source, basis)

    t1 = @elapsed mlfma_op = MLFMAOperator(cfie, basis, 0.35λ)
    @printf("  MLFMA构建: %.1fs\n", t1)
    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart = 50, maxiter = 200, tol = 1e-3, verbose = true)
    t2 = @elapsed I_mlfma = solve!(solver, mlfma_op, V; Pl = P)
    @printf("  GMRES: %.1fs\n", t2)
    mlfma_op = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I_mlfma, basis)
    rcs_phi0  = rcs_dB[:, 1]
    rcs_phi90 = rcs_dB[:, 2]
    θ_rad_deg = rad2deg.(θs_obs)

    r1f = save_and_report_case("F6_CFIE_Sphere_MLFMA_phi0_vs_Feko",
        θ_deg_feko, rcs_phi0, feko_cuts[0.0].rcs_dBsm, 3.0)
    r2f = save_and_report_case("F6_CFIE_Sphere_MLFMA_phi90_vs_Feko",
        θ_deg_feko, rcs_phi90, feko_cuts[90.0].rcs_dBsm, 3.0)
    r1m = save_and_report_case("F6_CFIE_Sphere_MLFMA_phi0_vs_Mie",
        θ_rad_deg, rcs_phi0, mie_phi0, 3.0)
    r2m = save_and_report_case("F6_CFIE_Sphere_MLFMA_phi90_vs_Mie",
        θ_rad_deg, rcs_phi90, mie_phi90, 3.0)
    push!(all_results, r1f, r2f, r1m, r2m)
end

# ─── 汇总 ────────────────────────────────────────────────────────────────────
if !isempty(all_results)
    println()
    print_accuracy_report(all_results; title = "F5-F6 Sphere 散射精度报告")

    open(joinpath(RESULT_DIR, "F5_F6_sphere_report.md"), "w") do io
        println(io, "# F5–F6 Sphere 散射精度报告 (vs Feko + Mie)")
        println(io, "\n生成时间: $(Dates.now())\n")
        print_accuracy_report(io, all_results; title = "F5-F6")
    end
end
