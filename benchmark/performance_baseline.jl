"""
    performance_baseline.jl

Phase 8.0: 性能基线测量 (EMMoMSuite 侧)

对 6 个用例分阶段计时:
1. PEC Plate 300MHz — EFIE Direct (N≈2640)
2. Jet 100MHz — EFIE Direct (N=14559)
3. Jet 100MHz — EFIE MLFMA+GMRES (N=14559)
4. Sphere 600MHz — CFIE MLFMA+GMRES (N=26424)
5. Plate 1.2GHz — VEFIE Direct (N≈986)
6. Plate+Metal 1.2GHz — SCFIE Direct (N≈1071)

每个用例计时:
  - 网格读取 + 基函数构建
  - 阻抗矩阵组装 (或 MLFMA setup)
  - 求解 (LU 或 GMRES)
  - RCS 计算

用法:
  julia --project=. --threads=4 benchmark/performance_baseline.jl        # 全部
  julia --project=. --threads=4 benchmark/performance_baseline.jl 1 2    # 仅用例 1,2
"""

using EMMoMSuite
using LinearAlgebra
using SparseArrays
using Printf
using Statistics
using Dates

const MESH_DIR = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne", "meshfiles")
const REPORT_FILE = joinpath(@__DIR__, "..", "test_results", "PERFORMANCE_BASELINE.md")
const REPORT_DIR = joinpath(@__DIR__, "..", "test_results", "reports")
const REPORT_CSV_FILE = joinpath(REPORT_DIR, "PERFORMANCE_BASELINE.csv")

# --- 简单 RCS 观测点 (8 方向, 不做全球面) ---
const THETA_PERF = collect(LinRange(0, π, 37))  # 5° 间隔
const PHI_PERF   = [0.0, π/2]  # E-plane + H-plane

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

# =================================================================
#  计时结果结构
# =================================================================
mutable struct TimingResult
    case_name::String
    equation::String
    solver::String
    N::Int
    t_mesh::Float64      # 网格 + 基函数
    t_assembly::Float64  # Z 组装 / MLFMA setup
    t_precond::Float64   # 预条件器 (MLFMA only)
    t_solve::Float64     # LU / GMRES
    t_rcs::Float64       # RCS 后处理
    t_total::Float64     # 总耗时
    notes::String
end

function TimingResult(name, eq, solver, N)
    TimingResult(name, eq, solver, N, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "")
end

const RESULTS = TimingResult[]

# =================================================================
#  用例 1: PEC Plate 300MHz — EFIE Direct
# =================================================================
function bench_1_plate_efie_direct()
    println("\n" * "="^60)
    println("  Case 1: PEC Plate 300MHz — EFIE Direct")
    println("="^60)
    
    mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
    freq = 3e8
    
    r = TimingResult("Plate EFIE Direct", "EFIE", "LU", 0)
    
    # 网格 + 基函数
    r.t_mesh = @elapsed begin
        mesh = read_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        basis = RWGBasis(mesh)
    end
    r.N = num_basis(basis)
    println("  N = $(r.N), 网格+基函数: $(round(r.t_mesh, digits=3))s")
    
    efie = EFIE(freq)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    # 阻抗矩阵组装
    r.t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(efie, basis)
    end
    println("  Z 组装: $(round(r.t_assembly, digits=3))s, 大小 $(size(Z))")
    
    # LU 求解
    r.t_solve = @elapsed begin
        I = Z \ V
    end
    println("  LU 求解: $(round(r.t_solve, digits=3))s")
    
    # RCS
    r.t_rcs = @elapsed begin
        radarCrossSection(THETA_PERF, PHI_PERF, I, basis)
    end
    println("  RCS: $(round(r.t_rcs, digits=3))s")
    
    r.t_total = r.t_mesh + r.t_assembly + r.t_solve + r.t_rcs
    println("  总计: $(round(r.t_total, digits=3))s")
    
    push!(RESULTS, r)
    GC.gc()
end

# =================================================================
#  用例 2: Jet 100MHz — EFIE Direct
# =================================================================
function bench_2_jet_efie_direct()
    println("\n" * "="^60)
    println("  Case 2: Jet 100MHz — EFIE Direct")
    println("="^60)
    
    mesh_file = joinpath(MESH_DIR, "jet_100MHz.nas")
    freq = 1e8
    
    r = TimingResult("Jet EFIE Direct", "EFIE", "LU", 0)
    
    r.t_mesh = @elapsed begin
        mesh = read_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        basis = RWGBasis(mesh)
    end
    r.N = num_basis(basis)
    println("  N = $(r.N), 网格+基函数: $(round(r.t_mesh, digits=3))s")
    
    efie = EFIE(freq)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    r.t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(efie, basis)
    end
    println("  Z 组装: $(round(r.t_assembly, digits=3))s, 大小 $(size(Z))")
    
    r.t_solve = @elapsed begin
        I = Z \ V
    end
    println("  LU 求解: $(round(r.t_solve, digits=3))s")
    
    r.t_rcs = @elapsed begin
        radarCrossSection(THETA_PERF, PHI_PERF, I, basis)
    end
    println("  RCS: $(round(r.t_rcs, digits=3))s")
    
    r.t_total = r.t_mesh + r.t_assembly + r.t_solve + r.t_rcs
    println("  总计: $(round(r.t_total, digits=3))s")
    
    push!(RESULTS, r)
    
    # 计算 CFIE 用于用例 2b 对比 (CFIE=EFIE+MFIE)
    println("\n  --- Case 2b: Jet 100MHz — CFIE Direct ---")
    r2 = TimingResult("Jet CFIE Direct", "CFIE", "LU", r.N)
    r2.t_mesh = r.t_mesh
    
    cfie = CFIE(freq, 0.5)
    
    r2.t_assembly = @elapsed begin
        Z_cfie = assemble_impedance_matrix(cfie, basis)
    end
    println("  Z_CFIE 组装: $(round(r2.t_assembly, digits=3))s")
    
    r2.t_solve = @elapsed begin
        I_cfie = Z_cfie \ V
    end
    println("  LU 求解: $(round(r2.t_solve, digits=3))s")
    
    r2.t_rcs = r.t_rcs  # same complexity
    r2.t_total = r2.t_mesh + r2.t_assembly + r2.t_solve + r2.t_rcs
    r2.notes = "CFIE/EFIE assembly ratio = $(round(r2.t_assembly / r.t_assembly, digits=1))×"
    println("  CFIE/EFIE 组装比: $(round(r2.t_assembly / r.t_assembly, digits=1))×")
    println("  总计: $(round(r2.t_total, digits=3))s")
    
    push!(RESULTS, r2)
    GC.gc()
end

# =================================================================
#  用例 3: Jet 100MHz — EFIE MLFMA+GMRES
# =================================================================
function bench_3_jet_efie_mlfma()
    println("\n" * "="^60)
    println("  Case 3: Jet 100MHz — EFIE MLFMA+GMRES")
    println("="^60)
    
    mesh_file = joinpath(MESH_DIR, "jet_100MHz.nas")
    freq = 1e8
    
    r = TimingResult("Jet EFIE MLFMA", "EFIE", "MLFMA+GMRES", 0)
    
    r.t_mesh = @elapsed begin
        mesh = read_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        basis = RWGBasis(mesh)
    end
    r.N = num_basis(basis)
    println("  N = $(r.N)")
    
    efie = EFIE(freq)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    lambda = EMMoMSuite.Constants.c0 / freq
    leaf_size = 0.35 * lambda
    
    # MLFMA 构建 (octree + near-field Z + translation)
    r.t_assembly = @elapsed begin
        Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
    end
    println("  MLFMA setup: $(round(r.t_assembly, digits=3))s")
    
    # 预条件器 (Sparse LU of Z_near — UMFPACK, optimal for ill-conditioned EFIE)
    r.t_precond = @elapsed begin
        P_near = lu(Z_mlfma.Z_near)
    end
    P = LUPreconditioner(P_near)
    println("  预条件器 (LU of Z_near): $(round(r.t_precond, digits=3))s")
    
    # GMRES
    solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
    r.t_solve = @elapsed begin
        I_mlfma = solve!(solver, Z_mlfma, V; Pl=P)
    end
    println("  GMRES: $(round(r.t_solve, digits=3))s")
    
    r.t_rcs = @elapsed begin
        radarCrossSection(THETA_PERF, PHI_PERF, I_mlfma, basis)
    end
    println("  RCS: $(round(r.t_rcs, digits=3))s")
    
    r.t_total = r.t_mesh + r.t_assembly + r.t_precond + r.t_solve + r.t_rcs
    println("  总计: $(round(r.t_total, digits=3))s")
    
    push!(RESULTS, r)
    GC.gc()
end

# =================================================================
#  用例 4: Sphere 600MHz — CFIE MLFMA+GMRES
# =================================================================
function bench_4_sphere_cfie_mlfma()
    println("\n" * "="^60)
    println("  Case 4: Sphere 600MHz — CFIE MLFMA+GMRES")
    println("="^60)
    
    mesh_file = joinpath(MESH_DIR, "sphere_600MHz.nas")
    freq = 6e8
    
    r = TimingResult("Sphere CFIE MLFMA", "CFIE", "MLFMA+GMRES", 0)
    
    r.t_mesh = @elapsed begin
        mesh = read_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        basis = RWGBasis(mesh)
    end
    r.N = num_basis(basis)
    println("  N = $(r.N)")
    
    cfie = CFIE(freq, 0.5)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    
    lambda = EMMoMSuite.Constants.c0 / freq
    leaf_size = 0.35 * lambda
    
    r.t_assembly = @elapsed begin
        Z_mlfma = MLFMAOperator(cfie, basis, leaf_size)
    end
    println("  MLFMA setup: $(round(r.t_assembly, digits=3))s")
    
    r.t_precond = @elapsed begin
        P_near = lu(Z_mlfma.Z_near)
    end
    P = LUPreconditioner(P_near)
    println("  预条件器 (LU of Z_near): $(round(r.t_precond, digits=3))s")
    
    solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
    r.t_solve = @elapsed begin
        I_mlfma = solve!(solver, Z_mlfma, V; Pl=P)
    end
    println("  GMRES: $(round(r.t_solve, digits=3))s")
    
    r.t_rcs = @elapsed begin
        radarCrossSection(THETA_PERF, PHI_PERF, I_mlfma, basis)
    end
    println("  RCS: $(round(r.t_rcs, digits=3))s")
    
    r.t_total = r.t_mesh + r.t_assembly + r.t_precond + r.t_solve + r.t_rcs
    println("  总计: $(round(r.t_total, digits=3))s")
    
    push!(RESULTS, r)
    GC.gc()
end

# =================================================================
#  用例 5: Plate 1.2GHz — VEFIE Direct 
# =================================================================
function bench_5_plate_vefie_direct()
    println("\n" * "="^60)
    println("  Case 5: Plate 1.2GHz — VEFIE Direct")
    println("="^60)
    
    mesh_file = joinpath(MESH_DIR, "plate_1dot2GHz.nas")
    freq = 12e8
    
    r = TimingResult("Plate VEFIE Direct", "VEFIE", "LU", 0)
    
    r.t_mesh = @elapsed begin
        mesh = read_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        basis = SWGBasis(mesh)
    end
    r.N = num_basis(basis)
    println("  N = $(r.N)")
    
    n_tet = EMMoMSuite.CoreModule.num_elements(mesh)
    εr = ComplexF64(2.0 * (1 - 0.0002im))
    perms = fill(εr, n_tet)
    
    vefie = VEFIE(freq, perms)
    source = PlaneWave(freq, π/4, π, [0.0, 0.0, 1.0])
    V = excitation_vector(vefie, source, basis, perms)
    
    r.t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(vefie, basis)
    end
    println("  Z 组装: $(round(r.t_assembly, digits=3))s, 大小 $(size(Z))")
    
    r.t_solve = @elapsed begin
        I = Z \ V
    end
    println("  LU 求解: $(round(r.t_solve, digits=3))s")
    
    r.t_rcs = @elapsed begin
        radarCrossSection(THETA_PERF, PHI_PERF, I, basis, perms)
    end
    println("  RCS: $(round(r.t_rcs, digits=3))s")
    
    r.t_total = r.t_mesh + r.t_assembly + r.t_solve + r.t_rcs
    println("  总计: $(round(r.t_total, digits=3))s")
    
    push!(RESULTS, r)
    GC.gc()
end

# =================================================================
#  用例 6: Plate+Metal 1.2GHz — SCFIE Direct
# =================================================================
function bench_6_plate_scfie_direct()
    println("\n" * "="^60)
    println("  Case 6: Plate+Metal 1.2GHz — SCFIE Direct")
    println("="^60)
    
    mesh_file = joinpath(MESH_DIR, "plate_and_metal_1dot2GHz.nas")
    freq = 12e8
    
    r = TimingResult("PlateMetal SCFIE Direct", "SCFIE", "LU", 0)
    
    r.t_mesh = @elapsed begin
        surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=1.0)
        set_frequency!(freq)
        rwg_basis = RWGBasis(surf_mesh)
        swg_basis = SWGBasis(vol_mesh)
    end
    n_surf = num_basis(rwg_basis)
    n_vol  = num_basis(swg_basis)
    r.N = n_surf + n_vol
    println("  N = $(r.N) (RWG=$n_surf, SWG=$n_vol)")
    
    n_tet = EMMoMSuite.CoreModule.num_elements(vol_mesh)
    εr = ComplexF64(2.0 * (1 - 0.0002im))
    perms = fill(εr, n_tet)
    
    scfie = SCFIE(freq, perms; alpha=0.5)
    source = PlaneWave(freq, π/4, π, [0.0, 0.0, 1.0])
    V = excitation_vector(source, rwg_basis, swg_basis)
    
    r.t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(scfie, rwg_basis, swg_basis)
    end
    println("  Z 组装: $(round(r.t_assembly, digits=3))s, 大小 $(size(Z))")
    
    r.t_solve = @elapsed begin
        I = Z \ V
    end
    println("  LU 求解: $(round(r.t_solve, digits=3))s")
    
    I_surf = I[1:n_surf]
    r.t_rcs = @elapsed begin
        radarCrossSection(THETA_PERF, PHI_PERF, I_surf, rwg_basis)
    end
    println("  RCS: $(round(r.t_rcs, digits=3))s")
    
    r.t_total = r.t_mesh + r.t_assembly + r.t_solve + r.t_rcs
    println("  总计: $(round(r.t_total, digits=3))s")
    
    push!(RESULTS, r)
    GC.gc()
end

# =================================================================
#  报告生成
# =================================================================
function generate_report()
    mkpath(dirname(REPORT_FILE))
    mkpath(REPORT_DIR)
    
    open(REPORT_FILE, "w") do io
        println(io, "# EMMoMSuite 性能基线报告")
        println(io, "")
        println(io, "> 生成时间: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"))")
        println(io, "> Julia 版本: $(VERSION)")
        println(io, "> 线程数: $(Threads.nthreads())")
        println(io, "> 操作系统: $(Sys.iswindows() ? "Windows" : Sys.islinux() ? "Linux" : "macOS")")
        println(io, "")
        println(io, "## 汇总表")
        println(io, "")
        println(io, "| 用例 | N | 方程 | 求解器 | 网格+BF (s) | 组装 (s) | 预条件 (s) | 求解 (s) | RCS (s) | 总计 (s) | 备注 |")
        println(io, "|------|---|------|--------|------------|----------|-----------|---------|---------|---------|------|")
        
        for r in RESULTS
            @printf(io, "| %s | %d | %s | %s | %.2f | %.2f | %.2f | %.2f | %.2f | **%.2f** | %s |\n",
                    r.case_name, r.N, r.equation, r.solver,
                    r.t_mesh, r.t_assembly, r.t_precond, r.t_solve, r.t_rcs, r.t_total, r.notes)
        end
        
        println(io, "")
        println(io, "## 关键比值")
        println(io, "")
        
        # 找 CFIE / EFIE 组装比
        efie_jet = findfirst(r -> r.case_name == "Jet EFIE Direct", RESULTS)
        cfie_jet = findfirst(r -> r.case_name == "Jet CFIE Direct", RESULTS)
        if efie_jet !== nothing && cfie_jet !== nothing
            ratio = RESULTS[cfie_jet].t_assembly / RESULTS[efie_jet].t_assembly
            @printf(io, "- CFIE/EFIE 组装比 (Jet N=14559): **%.1f×** (目标 ≤ 2.5×)\n", ratio)
        end
        
        # MLFMA setup 占比
        for r in RESULTS
            if r.solver == "MLFMA+GMRES"
                pct = r.t_assembly / r.t_total * 100
                @printf(io, "- %s: MLFMA setup 占比 **%.0f%%** (目标 ≤ 70%%)\n", r.case_name, pct)
            end
        end
        
        println(io, "")
        println(io, "## 与 Legacy 对比 (待补充)")
        println(io, "")
        println(io, "| 用例 | Legacy (s) | EMMoMSuite (s) | Ratio | 状态 |")
        println(io, "|------|-----------|-------------|-------|------|")
        for r in RESULTS
            println(io, "| $(r.case_name) | TBD | $(round(r.t_total, digits=2)) | TBD | — |")
        end
    end

    open(REPORT_CSV_FILE, "w") do io
        println(io, "case_name,equation,solver,N,t_mesh,t_assembly,t_precond,t_solve,t_rcs,t_total,notes")
        for r in RESULTS
            notes = replace(r.notes, "," => ";")
            @printf(io, "%s,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%s\n",
                r.case_name, r.equation, r.solver, r.N,
                r.t_mesh, r.t_assembly, r.t_precond, r.t_solve, r.t_rcs, r.t_total, notes)
        end
    end
    
    println("\n报告已生成: $REPORT_FILE")
    println("性能 CSV 已生成: $REPORT_CSV_FILE")
end

# =================================================================
#  主入口
# =================================================================
function main()
    println("="^60)
    println("  EMMoMSuite 性能基线测量 (Phase 8.0)")
    println("  线程: $(Threads.nthreads()), Julia: $(VERSION)")
    println("="^60)
    
    # 解析命令行参数
    cases = if isempty(ARGS)
        [1, 2, 3, 4, 5, 6]
    else
        [parse(Int, a) for a in ARGS]
    end
    
    bench_funcs = [
        bench_1_plate_efie_direct,
        bench_2_jet_efie_direct,
        bench_3_jet_efie_mlfma,
        bench_4_sphere_cfie_mlfma,
        bench_5_plate_vefie_direct,
        bench_6_plate_scfie_direct,
    ]
    
    for c in cases
        if 1 ≤ c ≤ length(bench_funcs)
            bench_funcs[c]()
        else
            println("未知用例: $c")
        end
    end
    
    # 生成报告
    if !isempty(RESULTS)
        generate_report()
    end
end

main()
