"""
    generate_fullsphere_legacy_baselines.jl

全球面 RCS Legacy 基线生成器。
使用 MoM_AllinOne 在标准全球面采样网格上生成 8 个基线:
  1. SEFIE Direct  (Jet 100MHz)
  2. SEFIE MLFMA   (Jet 100MHz)
  3. SCFIE Direct  (Sphere 600MHz)
  4. SCFIE MLFMA   (Sphere 600MHz)
  5. VEFIE Direct  (plate 1.2GHz)
  6. VEFIE MLFMA   (plate 1.2GHz)
  7. VSEFIE Direct (plate+metal 1.2GHz)
  8. VSEFIE MLFMA  (plate+metal 1.2GHz)

全球面采样: θ ∈ [-π, π], 73 点 (5° 间隔) × φ ∈ [0, π), 18 条切面 (10° 间隔) = 1314 方向

依赖: MoM_AllinOne 需可加载
用法: julia generate_fullsphere_legacy_baselines.jl [case_name]
  不带参数: 运行全部 8 个用例
  带参数:   仅运行指定用例 (如 SEFIE_Direct_Jet)
"""

using MoM_AllinOne
using DelimitedFiles
using Printf

# ============================================================
#  路径配置
# ============================================================
const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "fullsphere")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

# ============================================================
#  标准全球面采样网格
# ============================================================

"""标准全球面 RCS 观测角 (Float32 版)"""
function standard_fullsphere_grid_f32()
    θs = LinRange{Float32}(-π, π, 73)   # 5° 间隔
    ϕs = LinRange{Float32}(0, Float32(170/180*π), 18)  # 0°:10°:170°
    return θs, ϕs
end

"""标准全球面 RCS 观测角 (Float64 版)"""
function standard_fullsphere_grid_f64()
    θs = LinRange{Float64}(-π, π, 73)
    ϕs = LinRange{Float64}(0, 170/180*π, 18)
    return θs, ϕs
end

# ============================================================
#  测试用例定义
# ============================================================

struct LegacyCase
    name::String
    mesh_file::String
    frequency::Float64
    ie_type::Symbol       # :EFIE or :CFIE
    sbfT::Symbol          # :RWG or :nothing
    vbfT::Symbol          # :SWG or :nothing
    solver_type::Symbol   # :direct or :gmres
    precision::Type       # Float32 or Float64
    source_theta::Float64 # PlaneWave θ 入射方向
    rtol::Float64
    restart::Int
end

const ALL_CASES = [
    # --- 表面 EFIE (Jet, 开体) ---
    LegacyCase("SEFIE_Direct_Jet",  "jet_100MHz.nas",  1e8,  :EFIE, :RWG, :nothing, :direct, Float32, π/2, 1e-3, 50),
    LegacyCase("SEFIE_MLFMA_Jet",   "jet_100MHz.nas",  1e8,  :EFIE, :RWG, :nothing, :gmres,  Float32, π/2, 1e-3, 50),

    # --- 表面 CFIE (Sphere, 闭体) ---
    LegacyCase("SCFIE_Direct_Sphere", "sphere_600MHz.nas", 6e8, :CFIE, :RWG, :nothing, :direct, Float32, π/2, 1e-3, 50),
    LegacyCase("SCFIE_MLFMA_Sphere", "sphere_600MHz.nas", 6e8, :CFIE, :RWG, :nothing, :gmres,  Float32, π/2, 1e-3, 50),

    # --- 体 EFIE (介质板) ---
    LegacyCase("VEFIE_Direct_Plate",  "plate_1dot2GHz.nas", 12e8, :EFIE, :nothing, :SWG, :direct, Float32, π/4, 1e-3, 50),
    LegacyCase("VEFIE_MLFMA_Plate",   "plate_1dot2GHz.nas", 12e8, :EFIE, :nothing, :SWG, :gmres,  Float64, π/4, 1e-3, 50),

    # --- 面体 EFIE (PEC+介质) ---
    LegacyCase("VSEFIE_Direct_PlateMetal", "plate_and_metal_1dot2GHz.nas", 12e8, :EFIE, :RWG, :SWG, :direct, Float32, π/4, 1e-3, 50),
    LegacyCase("VSEFIE_MLFMA_PlateMetal",  "plate_and_metal_1dot2GHz.nas", 12e8, :EFIE, :RWG, :SWG, :gmres,  Float32, π/4, 1e-3, 50),
]

# ============================================================
#  运行单个用例
# ============================================================

function run_legacy_case(case::LegacyCase)
    println("\n" * "="^72)
    println("  Legacy Case: $(case.name)")
    println("  Mesh: $(case.mesh_file), Freq: $(case.frequency/1e6) MHz")
    println("  IE: $(case.ie_type), Solver: $(case.solver_type), Precision: $(case.precision)")
    println("="^72)

    # 选择全球面采样网格
    if case.precision == Float64
        θs_obs, ϕs_obs = standard_fullsphere_grid_f64()
    else
        θs_obs, ϕs_obs = standard_fullsphere_grid_f32()
    end

    filename = joinpath(MOM_ALLINONE_DIR, "meshfiles", case.mesh_file)
    if !isfile(filename)
        @warn "Mesh file not found: $filename — skipping"
        return nothing
    end

    solver_script = case.solver_type == :direct ? "direct_solver.jl" : "fast_solver.jl"
    solver_path = joinpath(MOM_ALLINONE_DIR, "src", solver_script)

    # 用 FT 版本的 PlaneWave
    if case.precision == Float64
        source = PlaneWave(case.source_theta, 0, 0.0, 1.0)
    else
        source = PlaneWave(Float32(case.source_theta), 0f0, 0f0, 1f0)
    end

    t_start = time()

    # 通过 @eval Main 在全局作用域中设置变量并 include 求解脚本
    result = @eval Main begin
        # 设置精度
        setPrecision!($(case.precision))
        SimulationParams.SHOWIMAGE = false

        # 参数
        filename  = $filename
        meshUnit  = :m
        frequency = $(case.frequency)
        ieT       = $(QuoteNode(case.ie_type))
        sbfT      = $(QuoteNode(case.sbfT))
        vbfT      = $(QuoteNode(case.vbfT))
        solverT   = $(QuoteNode(case.solver_type))
        rtol      = $(case.rtol)
        restart   = $(case.restart)
        source    = $source
        θs_obs    = $θs_obs
        ϕs_obs    = $ϕs_obs

        include($solver_path)
    end

    t_total = time() - t_start

    # result = (RCSθsϕs, RCSθsϕsdB, RCS, RCSdB)
    # RCS: (Nθ × Nφ), 线性值
    RCS_total = result[3]     # 线性
    RCS_dB    = result[4]     # dB
    RCSθsϕs   = result[1]    # (2, Nθ, Nφ)

    # 保存全球面 CSV
    save_fullsphere_csv(case.name, θs_obs, ϕs_obs, RCSθsϕs, RCS_total, RCS_dB)

    @printf("  ✓ 完成: %.1f 秒, %d×%d = %d 观测方向\n",
            t_total, length(θs_obs), length(ϕs_obs), length(θs_obs)*length(ϕs_obs))

    return RCS_dB
end

# ============================================================
#  保存结果
# ============================================================

function save_fullsphere_csv(case_name, θs_obs, ϕs_obs, RCSθsϕs, RCS_total, RCS_dB)
    output_file = joinpath(OUTPUT_DIR, "$(case_name).csv")

    # 构建平铺表格: theta_deg, phi_deg, RCS_theta_dB, RCS_phi_dB, RCS_total_dB
    Nθ = length(θs_obs)
    Nφ = length(ϕs_obs)
    n_rows = Nθ * Nφ

    data = Matrix{Float64}(undef, n_rows, 5)
    idx = 0
    for jφ in 1:Nφ, iθ in 1:Nθ
        idx += 1
        rcs_θ = max(RCSθsϕs[1, iθ, jφ], 1e-30)  # 防止 log(0)
        rcs_φ = max(RCSθsϕs[2, iθ, jφ], 1e-30)
        rcs_t = max(RCS_total[iθ, jφ], 1e-30)

        data[idx, 1] = rad2deg(Float64(θs_obs[iθ]))
        data[idx, 2] = rad2deg(Float64(ϕs_obs[jφ]))
        data[idx, 3] = 10log10(Float64(rcs_θ))
        data[idx, 4] = 10log10(Float64(rcs_φ))
        data[idx, 5] = 10log10(Float64(rcs_t))
    end

    # 写入 CSV (带表头)
    open(output_file, "w") do io
        println(io, "theta_deg,phi_deg,RCS_theta_dB,RCS_phi_dB,RCS_total_dB")
        for i in 1:n_rows
            @printf(io, "%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    data[i, 1], data[i, 2], data[i, 3], data[i, 4], data[i, 5])
        end
    end

    println("  → 保存: $output_file ($n_rows 行)")
end

# ============================================================
#  主入口
# ============================================================

function main()
    # 可选命令行参数指定运行哪些 case
    if length(ARGS) > 0
        target_names = Set(ARGS)
        cases = filter(c -> c.name in target_names, ALL_CASES)
        if isempty(cases)
            println("可用 case 名称:")
            for c in ALL_CASES
                println("  $(c.name)")
            end
            error("未找到匹配的 case: $(ARGS)")
        end
    else
        cases = ALL_CASES
    end

    println("\n全球面 Legacy 基线生成")
    println("输出目录: $OUTPUT_DIR")
    println("用例数量: $(length(cases))")
    println("采样网格: 73θ × 18φ = 1314 方向\n")

    results = Dict{String, Any}()
    for case in cases
        try
            rcs = run_legacy_case(case)
            results[case.name] = rcs
        catch e
            @warn "Case $(case.name) failed" exception=(e, catch_backtrace())
            results[case.name] = nothing
        end
        GC.gc()  # 释放内存
    end

    # 汇总
    println("\n" * "="^72)
    println("  基线生成完成")
    println("="^72)
    for case in cases
        status = results[case.name] === nothing ? "❌ FAILED" : "✅ OK"
        println("  $status  $(case.name)")
    end
end

main()
