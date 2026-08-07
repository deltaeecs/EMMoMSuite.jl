"""
run_F7_F9_plate.jl — Phase 14 精度基准：介质板散射 F7, F9

测试说明：
  F7  V-EFIE Direct   (SWG)     vs Feko plate_1dot2GHz     (RMSE < 2 dB)
  F9  SCFIE Direct    (RWG+SWG) vs Feko plate_metal_1dot2GHz (RMSE < 2 dB)

注意: F8 (纯介质板 SCFIE) 需要表面网格提取，当前网格文件 plate_1dot2GHz.nas
      仅含四面体元素 (CTETRA)，无法直接构建 RWGBasis，故 F8 本次跳过。

几何:
  F7  — plate_1dot2GHz.nas (介质板, 纯 CTETRA)
  F9  — plate_and_metal_1dot2GHz.nas (CTRIA3 PEC 表面 + CTETRA 介质)

频率: 1.2 GHz
介质: ε_r = 2.0 × (1 - 0.0002j)（低损耗介质板）

运行方式:
  julia --project=. benchmark/accuracy/run_F7_F9_plate.jl [F7] [F9]
  (不带参数时运行 F7)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using EMMoMSuite
using LinearAlgebra, Printf, Statistics, Dates

# ─── 路径 ─────────────────────────────────────────────────────────────────────
const ROOT_DIR   = joinpath(@__DIR__, "..", "..")
const MESH_DIR   = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "meshfiles")
const FEKO_DIR   = joinpath(ROOT_DIR, "..", "deps", "fixtures", "AllinOne", "deps", "compare_feko")
const RESULT_DIR = joinpath(ROOT_DIR, "test_results", "accuracy")
mkpath(RESULT_DIR)

enabled = isempty(ARGS) ? ["F7"] : ARGS

const θs_obs = collect(range(-π, π, length = 721))
const ϕs_obs = [0.0, π/2]

function load_feko_phi_cuts(feko_file)
    isfile(feko_file) || error("Feko 文件不存在: $feko_file")
    theta_deg, phi_deg, _, rcs_dBsm = read_feko_rcs(feko_file)
    return split_phi_cuts(theta_deg, phi_deg, rcs_dBsm)
end

function save_and_report_case(label, θdeg, rcs_ems, rcs_ref, threshold)
    result = compute_rcs_accuracy(rcs_ems, rcs_ref, θdeg, label; threshold)
    @printf("  %-48s  RMSE=%6.3f dB  %s\n", label, result.rmse_dB,
            result.pass ? "✓ PASS" : "✗ FAIL (门限 $(threshold) dB)")
    csv_path = joinpath(RESULT_DIR, label * ".csv")
    open(csv_path, "w") do io
        println(io, "theta_deg,rcs_ems_dBsm,rcs_feko_dBsm,diff_dB")
        for i in eachindex(θdeg)
            @printf(io, "%.2f,%.6f,%.6f,%.6f\n",
                    θdeg[i], rcs_ems[i], rcs_ref[i], rcs_ems[i] - rcs_ref[i])
        end
    end
    return result
end

freq = 1.2e9
λ    = 299792458.0 / freq
set_frequency!(freq)

εr = ComplexF64(2.0 * (1 - 0.0002im))
# Legacy PlaneWave(θ=π/4, ϕ=0, α=0) maps to propagation along
# (-sinθ, 0, -cosθ) with polarization -θhat = (-cosθ, 0, sinθ).
source = PlaneWave(freq, 3π / 4, π, [-1.0, 0.0, 1.0])

all_results = AccuracyResult[]

# ─── F7: V-EFIE Direct ────────────────────────────────────────────────────────
if "F7" in enabled
    println("\n[F7] V-EFIE Direct (SWG, plate_1dot2GHz)")
    mesh_file = joinpath(MESH_DIR, "plate_1dot2GHz.nas")
    isfile(mesh_file) || error("网格不存在: $mesh_file")

    mesh      = read_nas_mesh(mesh_file, scale = 1.0)
    swg_basis = SWGBasis(mesh)
    n_vol     = num_basis(swg_basis)
    n_tet     = EMMoMSuite.CoreModule.num_elements(mesh)
    println("  SWG 基函数: $n_vol, 四面体: $n_tet")

    perms = fill(εr, n_tet)
    vefie = VEFIE(freq, perms)

    V = excitation_vector(vefie, source, swg_basis, perms)

    t1 = @elapsed Z = assemble_impedance_matrix(vefie, swg_basis)
    @printf("  VEFIE 组装: %.1fs  (%d×%d)\n", t1, size(Z)...)
    t2 = @elapsed I = Z \ V
    @printf("  LU求解: %.1fs\n", t2)
    Z = nothing; GC.gc()

    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, swg_basis, perms)
    θ_deg = rad2deg.(θs_obs)

    feko_cuts = load_feko_phi_cuts(joinpath(FEKO_DIR, "plate_1dot2GHzRCS.csv"))
    r1 = save_and_report_case("F7_VEFIE_Plate_Direct_phi0_vs_Feko",
        feko_cuts[0.0].theta, rcs_dB[:, 1], feko_cuts[0.0].rcs_dBsm, 2.0)
    r2 = save_and_report_case("F7_VEFIE_Plate_Direct_phi90_vs_Feko",
        feko_cuts[90.0].theta, rcs_dB[:, 2], feko_cuts[90.0].rcs_dBsm, 2.0)
    push!(all_results, r1, r2)
end

# ─── F9: SCFIE Direct (Plate+Metal) ──────────────────────────────────────────
if "F9" in enabled
    println("\n[F9] SCFIE Direct (RWG+SWG, plate_and_metal_1dot2GHz)")
    mesh_file = joinpath(MESH_DIR, "plate_and_metal_1dot2GHz.nas")
    isfile(mesh_file) || error("网格不存在: $mesh_file")

    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale = 1.0)
    rwg_basis = RWGBasis(surf_mesh)
    swg_basis = SWGBasis(vol_mesh)
    n_surf = num_basis(rwg_basis)
    n_vol  = num_basis(swg_basis)
    n_tet  = EMMoMSuite.CoreModule.num_elements(vol_mesh)
    println("  RWG: $n_surf, SWG: $n_vol, Total: $(n_surf+n_vol)")

    perms = fill(εr, n_tet)
    scfie = SCFIE(freq, perms; alpha = 0.5)

    V = excitation_vector(source, rwg_basis, swg_basis)

    t1 = @elapsed Z = assemble_impedance_matrix(scfie, rwg_basis, swg_basis)
    @printf("  SCFIE 组装: %.1fs  (%d×%d)\n", t1, size(Z)...)
    t2 = @elapsed I = Z \ V
    @printf("  LU求解: %.1fs\n", t2)
    Z = nothing; GC.gc()

    # SCFIE RCS (合并表面+体积远场)
    _, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, rwg_basis, swg_basis, perms)

    feko_cuts = load_feko_phi_cuts(joinpath(FEKO_DIR, "plate_metal_1dot2GHzRCS.csv"))
    r1 = save_and_report_case("F9_SCFIE_PlateMetal_Direct_phi0_vs_Feko",
        feko_cuts[0.0].theta, rcs_dB[:, 1], feko_cuts[0.0].rcs_dBsm, 2.0)
    r2 = save_and_report_case("F9_SCFIE_PlateMetal_Direct_phi90_vs_Feko",
        feko_cuts[90.0].theta, rcs_dB[:, 2], feko_cuts[90.0].rcs_dBsm, 2.0)
    push!(all_results, r1, r2)
end

# ─── F8 跳过说明 ─────────────────────────────────────────────────────────────
if "F8" in enabled
    @warn "F8 (SCFIE on 纯介质板) 已跳过。原因: plate_1dot2GHz.nas 仅含 CTETRA，" *
          "缺少 RWGBasis 所需的表面三角网格 (CTRIA3)。" *
          "如需测试，需先从 CTETRA 提取外表面三角面，或提供单独的表面网格文件。"
end

# ─── 汇总 ────────────────────────────────────────────────────────────────────
if !isempty(all_results)
    println()
    print_accuracy_report(all_results; title = "F7/F9 介质板散射精度报告")

    open(joinpath(RESULT_DIR, "F7_F9_plate_report.md"), "w") do io
        println(io, "# F7/F9 介质板散射精度报告 (vs Feko)")
        println(io, "\n生成时间: $(Dates.now())\n")
        println(io, "> F8 (纯介质 SCFIE) 需要表面网格提取，当前未实现，已跳过。\n")
        print_accuracy_report(io, all_results; title = "F7/F9")
    end
end
