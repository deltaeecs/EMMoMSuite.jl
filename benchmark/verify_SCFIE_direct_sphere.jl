"""
    verify_SCFIE_direct_sphere.jl

C1: S-CFIE Direct Sphere 验证脚本。
在 sphere_600MHz 网格上运行 CFIE Direct 求解 (alpha=0.6 匹配 Legacy)，
与 Legacy SCFIE_Direct_Sphere.csv 基线对比。

用法: julia --project=. benchmark/verify_SCFIE_direct_sphere.jl
"""

using EMSuite
using LinearAlgebra
using Printf
using Statistics

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")
const LEGACY_BASELINE_DIR = joinpath(@__DIR__, "..", "test_results", "legacy_baseline")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "test_results", "emsuite_verification")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

function verify_SCFIE_direct_sphere()
    println("=" ^ 72)
    println("  C1: S-CFIE Direct — Sphere 600MHz (闭体)")
    println("=" ^ 72)

    # 1. 参数
    freq = 6e8  # 600 MHz
    alpha = 0.6  # 匹配 Legacy 默认值

    # 2. 网格
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "sphere_600MHz.nas")
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file")
    end
    println("  Mesh: sphere_600MHz.nas")

    mesh = read_nas_mesh(mesh_file, scale=1.0)
    println("  Elements: $(num_elements(mesh))")

    # 3. 基函数
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("  Unknowns (RWG): $N")

    # 4. CFIE (alpha = 0.6)
    println("  CFIE alpha = $alpha")
    cfie = CFIE(freq, alpha)

    # 5. 组装
    println("  Assembling Z matrix...")
    t_asm = @elapsed Z = assemble_impedance_matrix(cfie, basis)
    @printf("  Z assembled: %.1f s, size %d x %d\n", t_asm, size(Z)...)

    # 6. 激励 (同 Legacy: θ=π/2, ϕ=π → -x 方向入射, z 极化)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    println("  |V|_max = $(maximum(abs.(V)))")

    # 7. 求解
    println("  Solving (LU)...")
    t_solve = @elapsed I = Z \ V
    @printf("  Solved: %.1f s\n", t_solve)

    # 8. RCS — 匹配 Legacy 采样: [-π, π], 721 点, 0.5° 间隔
    println("  Computing RCS...")
    θs_obs = collect(LinRange(-π, π, 721))
    ϕs_obs = [0.0, π/2]

    t_rcs = @elapsed RCS_res = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    @printf("  RCS computed: %.1f s\n", t_rcs)

    # 解析返回值
    RCS_total = RCS_res[2]  # (Nθ, Nϕ) 线性
    RCS_Phi0_linear = abs.(RCS_total[:, 1])
    RCS_Phi90_linear = abs.(RCS_total[:, 2])
    RCS_Phi0_dB = 10 .* log10.(max.(RCS_Phi0_linear, 1e-30))
    RCS_Phi90_dB = 10 .* log10.(max.(RCS_Phi90_linear, 1e-30))

    # 保存
    θ_deg = rad2deg.(θs_obs)
    output_file = joinpath(OUTPUT_DIR, "SCFIE_Direct_Sphere.csv")
    open(output_file, "w") do io
        println(io, "Theta_Deg,RCS_Phi0_dB,RCS_Phi90_dB")
        for i in 1:length(θs_obs)
            @printf(io, "%.6f,%.10f,%.10f\n", θ_deg[i], RCS_Phi0_dB[i], RCS_Phi90_dB[i])
        end
    end
    println("  Saved: $output_file")

    # ========================================================
    #  与 Legacy 基线对比
    # ========================================================
    legacy_file = joinpath(LEGACY_BASELINE_DIR, "SCFIE_Direct_Sphere.csv")
    if !isfile(legacy_file)
        println("  Legacy baseline not found: $legacy_file")
        return
    end

    println("\n  --- Legacy 对比 ---")
    legacy_lines = readlines(legacy_file)
    n_legacy = length(legacy_lines) - 1
    legacy_phi0 = Vector{Float64}(undef, n_legacy)
    legacy_phi90 = Vector{Float64}(undef, n_legacy)
    legacy_theta = Vector{Float64}(undef, n_legacy)

    for i in 1:n_legacy
        fields = split(legacy_lines[i+1], ',')
        legacy_theta[i] = parse(Float64, fields[2])
        legacy_phi0[i] = parse(Float64, fields[3])
        legacy_phi90[i] = parse(Float64, fields[4])
    end

    @assert length(legacy_theta) == length(θ_deg) "采样点数不匹配: Legacy=$(length(legacy_theta)), EMSuite=$(length(θ_deg))"

    # 逐切面统计
    for (label, leg, ems) in [("Phi=0", legacy_phi0, RCS_Phi0_dB),
                               ("Phi=90", legacy_phi90, RCS_Phi90_dB)]
        diff = ems .- leg
        abs_diff = abs.(diff)
        rmse = sqrt(mean(diff.^2))
        max_d = maximum(abs_diff)
        mean_d = mean(diff)
        pct_lt1 = count(d -> d < 1.0, abs_diff) / length(abs_diff) * 100

        @printf("  %-8s  RMSE=%.3f  MaxDiff=%.3f  MeanDiff=%+.3f  <1dB=%.1f%%\n",
                label, rmse, max_d, mean_d, pct_lt1)

        # 找最大误差点
        idx = argmax(abs_diff)
        @printf("           最大误差: theta=%.1f deg, Legacy=%.2f, EMSuite=%.2f dB\n",
                legacy_theta[idx], leg[idx], ems[idx])
    end

    # 综合
    all_diff = vcat(RCS_Phi0_dB .- legacy_phi0, RCS_Phi90_dB .- legacy_phi90)
    rmse_all = sqrt(mean(all_diff.^2))
    max_all = maximum(abs.(all_diff))
    @printf("\n  综合  RMSE=%.3f  MaxDiff=%.3f dB\n", rmse_all, max_all)

    # 判定
    passed = rmse_all < 1.0 && max_all < 5.0
    println("  判定: $(passed ? "PASS" : "FAIL")")

    println("\n  Total time: $(round(t_asm + t_solve + t_rcs, digits=1)) s")
end

verify_SCFIE_direct_sphere()
