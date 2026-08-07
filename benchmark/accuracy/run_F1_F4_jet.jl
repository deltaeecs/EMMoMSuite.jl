"""
run_F1_F4_jet.jl — Phase 14 精度基准：Jet 散射 F1-F4

测试说明：
  F1  S-EFIE Direct  vs Feko  (RMSE < 2 dB)
  F2  S-CFIE Direct  vs Feko  (RMSE < 2 dB)
  F3  S-EFIE MLFMA   vs Feko  (RMSE < 3 dB)
  F4  S-CFIE MLFMA   vs Feko  (RMSE < 3 dB)

几何: Jet 战机 (PEC), 频率 = 100 MHz
基线: MoM_AllinOne/deps/compare_feko/jet_100MHzRCS.csv

运行方式:
  julia --project=. benchmark/accuracy/run_F1_F4_jet.jl [F1] [F2] [F3] [F4]
  (不带参数时运行 F1 和 F2)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using EMMoMSuite
using LinearAlgebra, Printf, Statistics, Dates

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

# ─── 路径 ─────────────────────────────────────────────────────────────────────
const ROOT_DIR  = joinpath(@__DIR__, "..", "..")
const MESH_DIR  = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "meshfiles")
const FEKO_DIR  = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "deps", "compare_feko")
const RESULT_DIR = joinpath(ROOT_DIR, "test_results", "accuracy")
mkpath(RESULT_DIR)

# ─── 命令行参数 ───────────────────────────────────────────────────────────────
enabled = isempty(ARGS) ? ["F1", "F2"] : ARGS

# ─── 观测角度 ─────────────────────────────────────────────────────────────────
const θs_obs = collect(range(-π, π, length = 721))
const ϕs_obs = [0.0, π/2]

# ─── 辅助：从 Feko 读取 φ=0° 和 φ=90° 切面 ───────────────────────────────────
function load_feko_phi_cuts(feko_file)
    isfile(feko_file) || error("Feko 文件不存在: $feko_file")
    theta_deg, phi_deg, _, rcs_dBsm = read_feko_rcs(feko_file)
    cuts = split_phi_cuts(theta_deg, phi_deg, rcs_dBsm)
    return cuts
end

# ─── 辅助：报告并保存单个用例结果 ─────────────────────────────────────────────
function save_and_report_case(label, θdeg, rcs_ems, rcs_ref, threshold; phi_tag = "")
    result = compute_rcs_accuracy(rcs_ems, rcs_ref, θdeg, label;
                                  threshold = threshold)
    @printf("  %-40s  RMSE=%6.3f dB  MaxErr=%6.3f dB  %s\n",
            label, result.rmse_dB, result.max_err_dB,
            result.pass ? "✓ PASS" : "✗ FAIL (门限 $(threshold) dB)")

    # 保存 CSV
    csv_path = joinpath(RESULT_DIR, label * ".csv")
    open(csv_path, "w") do io
        println(io, "theta_deg,rcs_emsuite_dBsm,rcs_feko_dBsm,diff_dB")
        for i in eachindex(θdeg)
            @printf(io, "%.2f,%.6f,%.6f,%.6f\n",
                    θdeg[i], rcs_ems[i], rcs_ref[i], rcs_ems[i] - rcs_ref[i])
        end
    end
    return result
end

# ─── 共用资源 ─────────────────────────────────────────────────────────────────
freq = 1.0e8
λ    = 299792458.0 / freq
set_frequency!(freq)

mesh_file = joinpath(MESH_DIR, "jet_100MHz.nas")
isfile(mesh_file) || error("网格文件不存在: $mesh_file")
mesh  = read_nas_mesh(mesh_file, scale = 1.0)
basis = RWGBasis(mesh)
N     = num_basis(basis)
println("Jet 100MHz: $N 未知量")

source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])

feko_file = joinpath(FEKO_DIR, "jet_100MHzRCS.csv")
feko_cuts = load_feko_phi_cuts(feko_file)
θ_deg_feko = feko_cuts[0.0].theta        # 标量 θ, 度

all_results = AccuracyResult[]

# ─── F1: S-EFIE Direct ───────────────────────────────────────────────────────
if "F1" in enabled
    println("\n[F1] S-EFIE Direct")
    efie = EFIE(freq)
    t1 = @elapsed Z = assemble_impedance_matrix(efie, basis)
    @printf("  组装: %.1fs\n", t1)
    V = excitation_vector(efie, source, basis)
    t2 = @elapsed I = Z \ V
    @printf("  LU求解: %.1fs\n", t2)
    Z = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    rcs_ems_phi0 = rcs_dB[:, 1]
    rcs_ems_phi90 = rcs_dB[:, 2]

    r1 = save_and_report_case("F1_SEFIE_Jet_Direct_phi0_vs_Feko",
        θ_deg_feko, rcs_ems_phi0, feko_cuts[0.0].rcs_dBsm, 2.0)
    r2 = save_and_report_case("F1_SEFIE_Jet_Direct_phi90_vs_Feko",
        θ_deg_feko, rcs_ems_phi90, feko_cuts[90.0].rcs_dBsm, 2.0)
    push!(all_results, r1, r2)
end

# ─── F2: S-CFIE Direct ───────────────────────────────────────────────────────
if "F2" in enabled
    println("\n[F2] S-CFIE Direct (α=0.5)")
    cfie = CFIE(freq, 0.5)
    t1 = @elapsed Z = assemble_impedance_matrix(cfie, basis)
    @printf("  组装: %.1fs\n", t1)
    V = excitation_vector(cfie, source, basis)
    t2 = @elapsed I = Z \ V
    @printf("  LU求解: %.1fs\n", t2)
    Z = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    rcs_ems_phi0  = rcs_dB[:, 1]
    rcs_ems_phi90 = rcs_dB[:, 2]

    r1 = save_and_report_case("F2_CFIE_Jet_Direct_phi0_vs_Feko",
        θ_deg_feko, rcs_ems_phi0, feko_cuts[0.0].rcs_dBsm, 2.0)
    r2 = save_and_report_case("F2_CFIE_Jet_Direct_phi90_vs_Feko",
        θ_deg_feko, rcs_ems_phi90, feko_cuts[90.0].rcs_dBsm, 2.0)
    push!(all_results, r1, r2)
end

# ─── F3: S-EFIE MLFMA ────────────────────────────────────────────────────────
if "F3" in enabled
    println("\n[F3] S-EFIE MLFMA+GMRES")
    efie = EFIE(freq)
    V = excitation_vector(efie, source, basis)

    t1 = @elapsed mlfma_op = MLFMAOperator(efie, basis, 0.35λ)
    @printf("  MLFMA构建: %.1fs\n", t1)

    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart = 50, maxiter = 200, tol = 1e-3, verbose = true)
    t2 = @elapsed I_mlfma = solve!(solver, mlfma_op, V; Pl = P)
    @printf("  GMRES: %.1fs\n", t2)
    mlfma_op = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I_mlfma, basis)
    r1 = save_and_report_case("F3_SEFIE_Jet_MLFMA_phi0_vs_Feko",
        θ_deg_feko, rcs_dB[:, 1], feko_cuts[0.0].rcs_dBsm, 3.0)
    r2 = save_and_report_case("F3_SEFIE_Jet_MLFMA_phi90_vs_Feko",
        θ_deg_feko, rcs_dB[:, 2], feko_cuts[90.0].rcs_dBsm, 3.0)
    push!(all_results, r1, r2)
end

# ─── F4: S-CFIE MLFMA ────────────────────────────────────────────────────────
if "F4" in enabled
    println("\n[F4] S-CFIE MLFMA+GMRES (α=0.5)")
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
    r1 = save_and_report_case("F4_CFIE_Jet_MLFMA_phi0_vs_Feko",
        θ_deg_feko, rcs_dB[:, 1], feko_cuts[0.0].rcs_dBsm, 3.0)
    r2 = save_and_report_case("F4_CFIE_Jet_MLFMA_phi90_vs_Feko",
        θ_deg_feko, rcs_dB[:, 2], feko_cuts[90.0].rcs_dBsm, 3.0)
    push!(all_results, r1, r2)
end

# ─── 汇总 ─────────────────────────────────────────────────────────────────────
if !isempty(all_results)
    println()
    print_accuracy_report(all_results; title = "F1-F4 Jet 散射精度报告")

    report_path = joinpath(RESULT_DIR, "F1_F4_jet_report.md")
    open(report_path, "w") do io
        println(io, "# F1–F4 Jet 散射精度报告 (vs Feko)")
        println(io, "\n生成时间: $(Dates.now())\n")
        print_accuracy_report(io, all_results; title = "F1-F4 Jet 散射精度")
    end
    println("报告已保存: $report_path")
end
