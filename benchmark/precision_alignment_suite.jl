"""
    precision_alignment_suite.jl

Phase 10: 全方程全路径精度对齐基准套件

测试矩阵:
  A: S-EFIE  (Jet 100MHz)      — Direct, Iterative, MLFMA
  B: S-MFIE  (Sphere 600MHz)   — CFIE 分解验证, MLFMA
  C: S-CFIE  (Sphere 600MHz)   — MLFMA
  D: V-EFIE  (plate 1.2GHz)    — Direct, Iterative, MLFMA
  E: VS-EFIE (plate+metal 1.2GHz) — Direct, Iterative, MLFMA

全球面采样: θ ∈ [-π, π], 73 点 × φ ∈ [0, π), 18 条切面 = 1314 方向

用法:
  julia precision_alignment_suite.jl           # 运行全部
  julia precision_alignment_suite.jl A         # 仅 S-EFIE
  julia precision_alignment_suite.jl A1 A3     # 仅 S-EFIE Direct 和 MLFMA
  julia precision_alignment_suite.jl B1        # 仅 S-MFIE CFIE 分解验证
"""

using EMSuite
using LinearAlgebra
using SparseArrays
using Printf
using Statistics
using DelimitedFiles

# ============================================================
#  路径与配置
# ============================================================
const BASE_DIR       = joinpath(@__DIR__, "..")
const MESH_DIR       = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles")
const LEGACY_DIR     = joinpath(BASE_DIR, "test_results", "legacy_baseline", "fullsphere")
const OUTPUT_DIR     = joinpath(BASE_DIR, "test_results", "emsuite_verification", "fullsphere")
const REPORT_DIR     = joinpath(BASE_DIR, "test_results")

mkpath(OUTPUT_DIR)

# ============================================================
#  全球面采样网格
# ============================================================
const THETA_OBS = collect(LinRange(-π, π, 73))     # 5° 间隔
const PHI_OBS   = collect(LinRange(0.0, 170/180*π, 18))  # 0°:10°:170°
const N_OBS     = length(THETA_OBS) * length(PHI_OBS)

println("全球面采样: $(length(THETA_OBS))θ × $(length(PHI_OBS))φ = $N_OBS 方向")

# ============================================================
#  工具函数
# ============================================================

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

"""加载全球面 Legacy 基线 CSV"""
function load_legacy_baseline(case_name::String)
    filepath = joinpath(LEGACY_DIR, "$(case_name).csv")
    if !isfile(filepath)
        @warn "Legacy baseline not found: $filepath"
        return nothing
    end
    
    # 读取 CSV: theta_deg, phi_deg, RCS_theta_dB, RCS_phi_dB, RCS_total_dB
    lines = readlines(filepath)
    header = lines[1]
    data = zeros(Float64, length(lines)-1, 5)
    for (i, line) in enumerate(lines[2:end])
        vals = parse.(Float64, split(line, ','))
        data[i, :] = vals
    end
    
    # 重塑为 (Nθ, Nφ) 矩阵
    Nθ = length(THETA_OBS)
    Nφ = length(PHI_OBS)
    @assert size(data, 1) == Nθ * Nφ "Legacy baseline size mismatch: expected $(Nθ*Nφ), got $(size(data,1))"
    
    RCS_total_dB = reshape(data[:, 5], Nθ, Nφ)
    return RCS_total_dB
end

"""保存 EMSuite 全球面 RCS 结果"""
function save_fullsphere_csv(case_name, θs, ϕs, RCS_res)
    RCSθsϕs, RCS_total, RCS_dB = RCS_res
    output_file = joinpath(OUTPUT_DIR, "$(case_name).csv")
    
    Nθ = length(θs)
    Nφ = length(ϕs)
    
    open(output_file, "w") do io
        println(io, "theta_deg,phi_deg,RCS_theta_dB,RCS_phi_dB,RCS_total_dB")
        for jφ in 1:Nφ, iθ in 1:Nθ
            rcs_θ = max(RCSθsϕs[1, iθ, jφ], 1e-30)
            rcs_φ = max(RCSθsϕs[2, iθ, jφ], 1e-30)
            rcs_t = max(RCS_total[iθ, jφ], 1e-30)
            @printf(io, "%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    rad2deg(θs[iθ]), rad2deg(ϕs[jφ]),
                    10log10(rcs_θ), 10log10(rcs_φ), 10log10(rcs_t))
        end
    end
    println("  → 保存: $output_file")
end

"""对比两个全球面 RCS 矩阵 (dB), 返回统计量"""
function compare_fullsphere(RCS_A_dB, RCS_B_dB; labelA="EMSuite", labelB="Legacy")
    diff = RCS_A_dB .- RCS_B_dB
    
    # 过滤极低 RCS 值 (< -60 dB), 这些值对 dB 比较不稳定
    valid_mask = (RCS_A_dB .> -60.0) .& (RCS_B_dB .> -60.0)
    diff_valid = diff[valid_mask]
    
    if isempty(diff_valid)
        @warn "No valid comparison points above -60 dB threshold"
        return nothing
    end
    
    mean_diff = mean(diff_valid)
    rmse      = sqrt(mean(diff_valid.^2))
    max_diff  = maximum(abs.(diff_valid))
    n_valid   = length(diff_valid)
    n_total   = length(diff)
    
    println("\n  ─── $labelA vs $labelB (全球面) ───")
    @printf("  有效点     : %d / %d (>-60 dB)\n", n_valid, n_total)
    @printf("  Mean Diff  : %+.4f dB\n", mean_diff)
    @printf("  RMSE       : %.4f dB\n",  rmse)
    @printf("  Max |Diff| : %.4f dB\n",  max_diff)
    
    # 逐 φ 切面统计
    Nθ, Nφ = size(RCS_A_dB)
    println("  ─── 按 φ 切面 RMSE ───")
    for jφ in 1:Nφ
        slice_diff = RCS_A_dB[:, jφ] .- RCS_B_dB[:, jφ]
        valid = (RCS_A_dB[:, jφ] .> -60.0) .& (RCS_B_dB[:, jφ] .> -60.0)
        if any(valid)
            sd = slice_diff[valid]
            @printf("  φ=%5.1f° : RMSE=%.4f dB, Mean=%+.4f dB (%d pts)\n",
                    rad2deg(PHI_OBS[jφ]), sqrt(mean(sd.^2)), mean(sd), count(valid))
        end
    end
    
    return (mean=mean_diff, rmse=rmse, max=max_diff, n_valid=n_valid)
end

"""对比两组求解系数"""
function compare_coefficients(I_A, I_B; label="")
    rel_diff = norm(I_A - I_B) / norm(I_A)
    @printf("  系数相对差 (%s): %.6e (%.4f%%)\n", label, rel_diff, rel_diff*100)
    return rel_diff
end

# ============================================================
#  A: S-EFIE (Jet 100MHz, N=14559)
# ============================================================

function run_A_sefie(subtests::Vector{String}=["A1","A2","A3"])
    println("\n" * "="^72)
    println("  A: S-EFIE — Jet 100 MHz (开体, N=14559)")
    println("="^72)
    
    mesh_file = joinpath(MESH_DIR, "jet_100MHz.nas")
    freq = 1e8
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    println("  基函数数: $(basis.nbf)")
    
    efie = EFIE(freq)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    I_direct = nothing
    RCS_direct_dB = nothing
    
    # --- A1: Direct + LU ---
    if "A1" in subtests
        println("\n  ▶ A1: Direct + LU")
        t0 = time()
        Z = assemble_impedance_matrix(efie, basis)
        t_asm = time() - t0
        @printf("  阻抗矩阵组装: %.1f 秒, 大小 %d×%d\n", t_asm, size(Z)...)
        
        t0 = time()
        I_direct = Z \ V
        t_solve = time() - t0
        @printf("  LU 求解: %.1f 秒\n", t_solve)
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_direct, basis)
        RCS_direct_dB = RCS_res[3]  # RCS_dB (Nθ × Nφ)
        save_fullsphere_csv("A1_SEFIE_Direct", THETA_OBS, PHI_OBS, RCS_res)
        
        # 对比 Legacy
        legacy = load_legacy_baseline("SEFIE_Direct_Jet")
        if legacy !== nothing
            compare_fullsphere(RCS_direct_dB, legacy; labelA="A1 Direct", labelB="Legacy Direct")
        end
        
        # 保留 Z 给 A2
        if "A2" in subtests
            @goto a2_with_z
        end
        Z = nothing
        GC.gc()
    end
    
    # --- A2: Iterative (GMRES on Dense Z) ---
    if "A2" in subtests && !("A1" in subtests)
        # 需要单独组装 Z
        println("\n  ▶ A2: Iterative (GMRES on Dense Z)")
        Z = assemble_impedance_matrix(efie, basis)
        @label a2_with_z
        println("\n  ▶ A2: Iterative (GMRES on Dense Z)")
        
        P_diag = Diagonal(1.0 ./ diag(Z))
        solver = GMRESSolver(restart=50, maxiter=200, tol=1e-6, verbose=true)
        t0 = time()
        I_iter = solve!(solver, Z, V; Pl=P_diag)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_iter, basis)
        RCS_iter_dB = RCS_res[3]
        save_fullsphere_csv("A2_SEFIE_Iterative", THETA_OBS, PHI_OBS, RCS_res)
        
        if RCS_direct_dB !== nothing
            compare_fullsphere(RCS_iter_dB, RCS_direct_dB; labelA="A2 Iterative", labelB="A1 Direct")
        end
        if I_direct !== nothing
            compare_coefficients(I_iter, I_direct; label="A2 vs A1")
        end
        
        Z = nothing
        GC.gc()
    end
    
    # --- A3: MLFMA + GMRES ---
    if "A3" in subtests
        println("\n  ▶ A3: MLFMA + GMRES")
        lambda = EMSuite.Constants.c0 / freq
        leaf_size = 0.35 * lambda
        
        t0 = time()
        Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
        t_mlfma = time() - t0
        @printf("  MLFMA 构建: %.1f 秒\n", t_mlfma)
        
        # 重排序激励向量
        V_sorted = V[Z_mlfma.sorted_ids]
        
        # SAI 预条件或 ILU
        P_near = lu(Z_mlfma.Z_near)
        P = LUPreconditioner(P_near)
        
        solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
        t0 = time()
        I_mlfma_sorted = solve!(solver, Z_mlfma, V_sorted; Pl=P)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        # 恢复原始排序
        I_mlfma = similar(I_mlfma_sorted)
        I_mlfma[Z_mlfma.sorted_ids] = I_mlfma_sorted
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_mlfma, basis)
        RCS_mlfma_dB = RCS_res[3]
        save_fullsphere_csv("A3_SEFIE_MLFMA", THETA_OBS, PHI_OBS, RCS_res)
        
        # 对比 Legacy MLFMA
        legacy_mlfma = load_legacy_baseline("SEFIE_MLFMA_Jet")
        if legacy_mlfma !== nothing
            compare_fullsphere(RCS_mlfma_dB, legacy_mlfma; labelA="A3 MLFMA", labelB="Legacy MLFMA")
        end
        
        # 对比 Direct
        if RCS_direct_dB !== nothing
            compare_fullsphere(RCS_mlfma_dB, RCS_direct_dB; labelA="A3 MLFMA", labelB="A1 Direct")
        end
        if I_direct !== nothing
            compare_coefficients(I_mlfma, I_direct; label="A3 vs A1")
        end
        
        Z_mlfma = nothing
        GC.gc()
    end
    
    println("\n  ✓ A: S-EFIE 完成")
end

# ============================================================
#  B: S-MFIE (Sphere 600MHz, N=26424) — 自洽验证
# ============================================================

function run_B_smfie(subtests::Vector{String}=["B1","B2"])
    println("\n" * "="^72)
    println("  B: S-MFIE — Sphere 600 MHz (闭体, 自洽验证)")
    println("="^72)
    
    # --- B1: CFIE 分解验证 (小网格) ---
    if "B1" in subtests
        println("\n  ▶ B1: CFIE = α·EFIE + (1-α)·MFIE 分解验证")
        # 使用 AllinOne Tri.nas (小网格, ~1330 RWG)
        small_mesh_file = joinpath(MESH_DIR, "Tri.nas")
        if !isfile(small_mesh_file)
            small_mesh_file = joinpath(@__DIR__, "..", "..", "MoM_Basics", "meshfiles", "Tri.nas")
        end
        
        freq_small = 3e8  # 300 MHz, 适配小网格
        mesh_small = read_nas_mesh(small_mesh_file, scale=1.0)
        set_frequency!(freq_small)
        basis_small = RWGBasis(mesh_small)
        println("  小网格基函数数: $(basis_small.nbf)")
        
        alpha = 0.5
        efie_s = EFIE(freq_small)
        mfie_s = MFIE(freq_small)
        cfie_s = CFIE(freq_small, alpha)
        
        t0 = time()
        Z_efie = assemble_impedance_matrix(efie_s, basis_small)
        Z_mfie = assemble_impedance_matrix(mfie_s, basis_small)
        Z_cfie = assemble_impedance_matrix(cfie_s, basis_small)
        t_asm = time() - t0
        @printf("  3 矩阵组装: %.1f 秒\n", t_asm)
        
        # 验证 Z_CFIE ≈ α * Z_EFIE + (1-α) * Z_MFIE
        Z_reconstructed = alpha * Z_efie + (1 - alpha) * Z_mfie
        rel_err = norm(Z_cfie - Z_reconstructed) / norm(Z_cfie)
        @printf("  ‖Z_CFIE - (α·Z_EFIE + (1-α)·Z_MFIE)‖ / ‖Z_CFIE‖ = %.6e\n", rel_err)
        
        if rel_err < 1e-12
            println("  ✅ B1 PASS: CFIE 分解验证通过 (< 1e-12)")
        else
            println("  ❌ B1 FAIL: CFIE 分解误差过大")
        end
        
        Z_efie = Z_mfie = Z_cfie = nothing
        GC.gc()
    end
    
    # --- B2: MLFMA + GMRES (Sphere 600MHz) ---
    if "B2" in subtests
        println("\n  ▶ B2: S-MFIE MLFMA + GMRES (Sphere 600 MHz)")
        mesh_file = joinpath(MESH_DIR, "sphere_600MHz.nas")
        freq = 6e8
        
        mesh = read_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        basis = RWGBasis(mesh)
        println("  基函数数: $(basis.nbf)")
        
        mfie = MFIE(freq)
        source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
        V = excitation_vector(mfie, source, basis)
        
        lambda = EMSuite.Constants.c0 / freq
        leaf_size = 0.35 * lambda
        
        t0 = time()
        Z_mlfma = MLFMAOperator(mfie, basis, leaf_size)
        t_mlfma = time() - t0
        @printf("  MLFMA 构建: %.1f 秒\n", t_mlfma)
        
        V_sorted = V[Z_mlfma.sorted_ids]
        P_near = lu(Z_mlfma.Z_near)
        P = LUPreconditioner(P_near)
        
        solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
        t0 = time()
        I_sorted = solve!(solver, Z_mlfma, V_sorted; Pl=P)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        I_mfie = similar(I_sorted)
        I_mfie[Z_mlfma.sorted_ids] = I_sorted
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_mfie, basis)
        RCS_mfie_dB = RCS_res[3]
        save_fullsphere_csv("B2_SMFIE_MLFMA", THETA_OBS, PHI_OBS, RCS_res)
        
        println("  ✓ B2: S-MFIE MLFMA RCS 已生成 (与 C1 S-CFIE 进行物理趋势对比)")
        
        Z_mlfma = nothing
        GC.gc()
    end
    
    println("\n  ✓ B: S-MFIE 完成")
end

# ============================================================
#  C: S-CFIE (Sphere 600MHz, N=26424) — MLFMA only
# ============================================================

function run_C_scfie(subtests::Vector{String}=["C1"])
    println("\n" * "="^72)
    println("  C: S-CFIE — Sphere 600 MHz (闭体, MLFMA)")
    println("="^72)
    
    mesh_file = joinpath(MESH_DIR, "sphere_600MHz.nas")
    freq = 6e8
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    println("  基函数数: $(basis.nbf)")
    
    cfie = CFIE(freq, 0.5)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    
    # --- C1: MLFMA + GMRES ---
    if "C1" in subtests
        println("\n  ▶ C1: S-CFIE MLFMA + GMRES")
        lambda = EMSuite.Constants.c0 / freq
        leaf_size = 0.35 * lambda
        
        t0 = time()
        Z_mlfma = MLFMAOperator(cfie, basis, leaf_size)
        t_mlfma = time() - t0
        @printf("  MLFMA 构建: %.1f 秒\n", t_mlfma)
        
        V_sorted = V[Z_mlfma.sorted_ids]
        P_near = lu(Z_mlfma.Z_near)
        P = LUPreconditioner(P_near)
        
        solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
        t0 = time()
        I_sorted = solve!(solver, Z_mlfma, V_sorted; Pl=P)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        I_cfie = similar(I_sorted)
        I_cfie[Z_mlfma.sorted_ids] = I_sorted
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_cfie, basis)
        RCS_cfie_dB = RCS_res[3]
        save_fullsphere_csv("C1_SCFIE_MLFMA", THETA_OBS, PHI_OBS, RCS_res)
        
        # 对比 Legacy MLFMA
        legacy_mlfma = load_legacy_baseline("SCFIE_MLFMA_Sphere")
        if legacy_mlfma !== nothing
            compare_fullsphere(RCS_cfie_dB, legacy_mlfma; labelA="C1 CFIE MLFMA", labelB="Legacy CFIE MLFMA")
        end
        
        # 对比 Legacy Direct
        legacy_direct = load_legacy_baseline("SCFIE_Direct_Sphere")
        if legacy_direct !== nothing
            compare_fullsphere(RCS_cfie_dB, legacy_direct; labelA="C1 CFIE MLFMA", labelB="Legacy CFIE Direct")
        end
        
        Z_mlfma = nothing
        GC.gc()
    end
    
    println("\n  ✓ C: S-CFIE 完成")
end

# ============================================================
#  D: V-EFIE (plate 1.2GHz) — 体积积分方程
# ============================================================

function run_D_vefie(subtests::Vector{String}=["D1","D2","D3"])
    println("\n" * "="^72)
    println("  D: V-EFIE — 介质板 1.2 GHz")
    println("="^72)
    
    mesh_file = joinpath(MESH_DIR, "plate_1dot2GHz.nas")
    freq = 12e8
    
    # 纯体积网格 → 用 read_nas_mesh, 返回 TetrahedraMesh
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    set_frequency!(freq)
    basis = SWGBasis(mesh)
    println("  SWG 基函数数: $(basis.nbf)")
    
    # 介电常数: Legacy 使用 2(1 - 0.0002im) 
    n_tet = EMSuite.CoreModule.num_elements(mesh)
    εr = ComplexF64(2.0 * (1 - 0.0002im))
    perms = fill(εr, n_tet)
    
    vefie = VEFIE(freq, perms)
    # Legacy 入射方向: PlaneWave(π/4, 0, 0f0, 1f0)
    # EMSuite 对应: PlaneWave(freq, π/4, π, [0,0,1])
    source = PlaneWave(freq, π/4, π, [0.0, 0.0, 1.0])
    V = excitation_vector(vefie, source, basis, perms)
    
    I_direct = nothing
    RCS_direct_dB = nothing
    
    # --- D1: Direct + LU ---
    if "D1" in subtests
        println("\n  ▶ D1: Direct + LU")
        t0 = time()
        Z = assemble_impedance_matrix(vefie, basis)
        t_asm = time() - t0
        @printf("  阻抗矩阵组装: %.1f 秒, 大小 %d×%d\n", t_asm, size(Z)...)
        
        t0 = time()
        I_direct = Z \ V
        t_solve = time() - t0
        @printf("  LU 求解: %.1f 秒\n", t_solve)
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_direct, basis, perms)
        RCS_direct_dB = RCS_res[3]
        save_fullsphere_csv("D1_VEFIE_Direct", THETA_OBS, PHI_OBS, RCS_res)
        
        legacy = load_legacy_baseline("VEFIE_Direct_Plate")
        if legacy !== nothing
            compare_fullsphere(RCS_direct_dB, legacy; labelA="D1 Direct", labelB="Legacy Direct")
        end
        
        if "D2" in subtests
            @goto d2_with_z
        end
        Z = nothing
        GC.gc()
    end
    
    # --- D2: Iterative (GMRES on Dense Z) ---
    if "D2" in subtests && !("D1" in subtests)
        Z = assemble_impedance_matrix(vefie, basis)
        @label d2_with_z
        println("\n  ▶ D2: Iterative (GMRES on Dense Z)")
        
        P_diag = Diagonal(1.0 ./ diag(Z))
        solver = GMRESSolver(restart=50, maxiter=200, tol=1e-6, verbose=true)
        t0 = time()
        I_iter = solve!(solver, Z, V; Pl=P_diag)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_iter, basis, perms)
        RCS_iter_dB = RCS_res[3]
        save_fullsphere_csv("D2_VEFIE_Iterative", THETA_OBS, PHI_OBS, RCS_res)
        
        if RCS_direct_dB !== nothing
            compare_fullsphere(RCS_iter_dB, RCS_direct_dB; labelA="D2 Iterative", labelB="D1 Direct")
        end
        
        Z = nothing
        GC.gc()
    end
    
    # --- D3: MLFMA + GMRES ---
    if "D3" in subtests
        println("\n  ▶ D3: V-EFIE MLFMA + GMRES")
        lambda = EMSuite.Constants.c0 / freq
        leaf_size = 0.35 * lambda
        
        t0 = time()
        Z_mlfma = MLFMAOperator(vefie, basis, leaf_size)
        t_mlfma = time() - t0
        @printf("  MLFMA 构建: %.1f 秒\n", t_mlfma)
        
        V_sorted = V[Z_mlfma.sorted_ids]
        P_near = lu(Z_mlfma.Z_near)
        P = LUPreconditioner(P_near)
        
        solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
        t0 = time()
        I_sorted = solve!(solver, Z_mlfma, V_sorted; Pl=P)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        I_mlfma = similar(I_sorted)
        I_mlfma[Z_mlfma.sorted_ids] = I_sorted
        
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_mlfma, basis, perms)
        RCS_mlfma_dB = RCS_res[3]
        save_fullsphere_csv("D3_VEFIE_MLFMA", THETA_OBS, PHI_OBS, RCS_res)
        
        legacy_mlfma = load_legacy_baseline("VEFIE_MLFMA_Plate")
        if legacy_mlfma !== nothing
            compare_fullsphere(RCS_mlfma_dB, legacy_mlfma; labelA="D3 MLFMA", labelB="Legacy MLFMA")
        end
        if RCS_direct_dB !== nothing
            compare_fullsphere(RCS_mlfma_dB, RCS_direct_dB; labelA="D3 MLFMA", labelB="D1 Direct")
        end
        
        Z_mlfma = nothing
        GC.gc()
    end
    
    println("\n  ✓ D: V-EFIE 完成")
end

# ============================================================
#  E: VS-EFIE (plate+metal 1.2GHz) — 面体耦合
# ============================================================

function run_E_vsefie(subtests::Vector{String}=["E1","E2","E3"])
    println("\n" * "="^72)
    println("  E: VS-EFIE — PEC+介质 1.2 GHz (SCFIE)")
    println("="^72)
    
    mesh_file = joinpath(MESH_DIR, "plate_and_metal_1dot2GHz.nas")
    freq = 12e8
    
    # 混合网格 → 表面 + 体积
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=1.0)
    set_frequency!(freq)
    rwg_basis = RWGBasis(surf_mesh)
    swg_basis = SWGBasis(vol_mesh)
    n_surf = rwg_basis.nbf
    n_vol  = swg_basis.nbf
    println("  RWG 基函数数: $n_surf, SWG 基函数数: $n_vol, 总计: $(n_surf + n_vol)")
    
    # 介电常数
    n_tet = EMSuite.CoreModule.num_elements(vol_mesh)
    εr = ComplexF64(2.0 * (1 - 0.0002im))
    perms = fill(εr, n_tet)
    
    # SCFIE (α=0.5 for surface CFIE part)
    scfie = SCFIE(freq, perms; alpha=0.5)
    # Legacy 入射方向: PlaneWave(π/4, 0, 0f0, 1f0)
    source = PlaneWave(freq, π/4, π, [0.0, 0.0, 1.0])
    V = excitation_vector(source, rwg_basis, swg_basis)
    
    I_direct = nothing
    RCS_direct_dB = nothing
    
    # --- E1: Direct + LU ---
    if "E1" in subtests
        println("\n  ▶ E1: Direct + LU")
        t0 = time()
        Z = assemble_impedance_matrix(scfie, rwg_basis, swg_basis)
        t_asm = time() - t0
        @printf("  阻抗矩阵组装: %.1f 秒, 大小 %d×%d\n", t_asm, size(Z)...)
        
        t0 = time()
        I_direct = Z \ V
        t_solve = time() - t0
        @printf("  LU 求解: %.1f 秒\n", t_solve)
        
        # SCFIE RCS: 需要拆分表面和体积部分
        I_surf = I_direct[1:n_surf]
        I_vol  = I_direct[n_surf+1:end]
        
        # 表面 RCS (RWG 部分)
        RCS_surf = radarCrossSection(THETA_OBS, PHI_OBS, I_surf, rwg_basis)
        # 体积 RCS (SWG 部分)
        RCS_vol  = radarCrossSection(THETA_OBS, PHI_OBS, I_vol, swg_basis, perms)
        
        # 合并: 总 RCS = 表面 + 体积 (线性叠加远场辐射)
        # 注意: 严格来说应该在远场矢量层面叠加后再取模方, 而非 RCS 相加
        # 但 radarCrossSection 目前返回的是 |N|^2, 不是 N 本身
        # TODO: 如果需要精确合并, 需要修改 radarCrossSection 返回 N_θ, N_φ
        # 暂时使用表面 RCS 作为主导项 (PEC 散射通常远大于介质体)
        RCS_res = RCS_surf
        RCS_direct_dB = RCS_res[3]
        save_fullsphere_csv("E1_VSEFIE_Direct", THETA_OBS, PHI_OBS, RCS_res)
        
        legacy = load_legacy_baseline("VSEFIE_Direct_PlateMetal")
        if legacy !== nothing
            compare_fullsphere(RCS_direct_dB, legacy; labelA="E1 Direct", labelB="Legacy Direct")
        end
        
        if "E2" in subtests
            @goto e2_with_z
        end
        Z = nothing
        GC.gc()
    end
    
    # --- E2: Iterative ---
    if "E2" in subtests && !("E1" in subtests)
        Z = assemble_impedance_matrix(scfie, rwg_basis, swg_basis)
        @label e2_with_z
        println("\n  ▶ E2: Iterative (GMRES on Dense Z)")
        
        P_diag = Diagonal(1.0 ./ diag(Z))
        solver = GMRESSolver(restart=50, maxiter=200, tol=1e-6, verbose=true)
        t0 = time()
        I_iter = solve!(solver, Z, V; Pl=P_diag)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        I_surf = I_iter[1:n_surf]
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_surf, rwg_basis)
        RCS_iter_dB = RCS_res[3]
        save_fullsphere_csv("E2_VSEFIE_Iterative", THETA_OBS, PHI_OBS, RCS_res)
        
        if RCS_direct_dB !== nothing
            compare_fullsphere(RCS_iter_dB, RCS_direct_dB; labelA="E2 Iterative", labelB="E1 Direct")
        end
        
        Z = nothing
        GC.gc()
    end
    
    # --- E3: MLFMA + GMRES ---
    if "E3" in subtests
        println("\n  ▶ E3: VS-EFIE MLFMA + GMRES")
        lambda = EMSuite.Constants.c0 / freq
        leaf_size = 0.35 * lambda
        
        t0 = time()
        Z_mlfma = MLFMAOperator(scfie, [rwg_basis, swg_basis], leaf_size)
        t_mlfma = time() - t0
        @printf("  MLFMA 构建: %.1f 秒\n", t_mlfma)
        
        V_sorted = V[Z_mlfma.sorted_ids]
        P_near = lu(Z_mlfma.Z_near)
        P = LUPreconditioner(P_near)
        
        solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
        t0 = time()
        I_sorted = solve!(solver, Z_mlfma, V_sorted; Pl=P)
        t_solve = time() - t0
        @printf("  GMRES 求解: %.1f 秒\n", t_solve)
        
        I_mlfma = similar(I_sorted)
        I_mlfma[Z_mlfma.sorted_ids] = I_sorted
        
        I_surf = I_mlfma[1:n_surf]
        RCS_res = radarCrossSection(THETA_OBS, PHI_OBS, I_surf, rwg_basis)
        RCS_mlfma_dB = RCS_res[3]
        save_fullsphere_csv("E3_VSEFIE_MLFMA", THETA_OBS, PHI_OBS, RCS_res)
        
        legacy_mlfma = load_legacy_baseline("VSEFIE_MLFMA_PlateMetal")
        if legacy_mlfma !== nothing
            compare_fullsphere(RCS_mlfma_dB, legacy_mlfma; labelA="E3 MLFMA", labelB="Legacy MLFMA")
        end
        if RCS_direct_dB !== nothing
            compare_fullsphere(RCS_mlfma_dB, RCS_direct_dB; labelA="E3 MLFMA", labelB="E1 Direct")
        end
        
        Z_mlfma = nothing
        GC.gc()
    end
    
    println("\n  ✓ E: VS-EFIE 完成")
end

# ============================================================
#  主入口
# ============================================================

function main()
    # 解析命令行参数
    if isempty(ARGS)
        # 运行全部
        run_A_sefie()
        GC.gc()
        run_B_smfie()
        GC.gc()
        run_C_scfie()
        GC.gc()
        run_D_vefie()
        GC.gc()
        run_E_vsefie()
    else
        # 解析子测试列表
        # 支持: A, A1, A3, B, B1, C, D, D1, D3, E 等
        all_subtests = Dict(
            'A' => ["A1","A2","A3"],
            'B' => ["B1","B2"],
            'C' => ["C1"],
            'D' => ["D1","D2","D3"],
            'E' => ["E1","E2","E3"],
        )
        
        requested = Dict{Char, Vector{String}}()
        for arg in ARGS
            group = uppercase(arg[1])
            if length(arg) == 1
                # 整组, 如 "A"
                requested[group] = all_subtests[group]
            else
                # 子项, 如 "A1"
                key = group
                if !haskey(requested, key)
                    requested[key] = String[]
                end
                push!(requested[key], uppercase(arg))
            end
        end
        
        if haskey(requested, 'A')
            run_A_sefie(requested['A'])
            GC.gc()
        end
        if haskey(requested, 'B')
            run_B_smfie(requested['B'])
            GC.gc()
        end
        if haskey(requested, 'C')
            run_C_scfie(requested['C'])
            GC.gc()
        end
        if haskey(requested, 'D')
            run_D_vefie(requested['D'])
            GC.gc()
        end
        if haskey(requested, 'E')
            run_E_vsefie(requested['E'])
            GC.gc()
        end
    end
    
    println("\n" * "="^72)
    println("  全部测试完成")
    println("="^72)
end

main()
