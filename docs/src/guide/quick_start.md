# Quick Start

本指南演示如何用 EMSuite 计算 PEC 金属板的双站 RCS，
完整流程从加载网格到输出 RCS（单位：dBsm）。

## 1. 准备网格文件

EMSuite 支持 Nastran（`.nas`）和 Gmsh（`.msh`）格式。
本示例使用 `plate.nas`——一块位于 XY 平面的矩形金属板。

## 2. 仿真脚本

新建文件 `run_plate.jl`：

```julia
using EMSuite

# ─── 1. 频率与全局参数 ────────────────────────────────────────────────────
freq = 300e6        # 300 MHz
set_frequency!(freq)   # 更新全局 k0、η0（供 RCS/激励向量内部使用）

# ─── 2. 加载网格 ──────────────────────────────────────────────────────────
mesh  = read_nas_mesh("plate.nas")
println("网格: $(num_elements(mesh)) 个三角形，$(num_vertices(mesh)) 个节点")

# ─── 3. 基函数 ────────────────────────────────────────────────────────────
basis = RWGBasis(mesh)
println("未知量数: $(num_basis(basis))")

# ─── 4. 积分方程算子 + 阻抗矩阵 ──────────────────────────────────────────
efie = EFIE(freq)
Z    = assemble_impedance_matrix(efie, basis)

# ─── 5. 激励源（+z 方向入射，x 极化）──────────────────────────────────────
inc_wave = PlaneWave(freq, π/2, π, [1.0, 0.0, 0.0])
V        = excitation_vector(inc_wave, basis)

# ─── 6. 求解 ─────────────────────────────────────────────────────────────
I_coeff = solve!(LUSolver(), Z, V)

# ─── 7. RCS 计算 ─────────────────────────────────────────────────────────
θ_obs  = collect(range(0.0, π, length = 181))   # 0° ~ 180°
ϕ_obs  = [0.0]                                  # phi = 0 平面

_, rcs_total, rcs_db = radarCrossSection(θ_obs, ϕ_obs, I_coeff, basis)
println("单站 RCS (θ=90°): ", rcs_db[91, 1], " dBsm")

# ─── 8. 保存结果 ─────────────────────────────────────────────────────────
save_RCS_txt("plate_rcs.txt", θ_obs, ϕ_obs, rcs_db)
save_results_hdf5("plate_results.h5", I_coeff, rcs_db)
println("完成！结果已保存至 plate_rcs.txt 和 plate_results.h5")
```

## 3. 运行

```bash
julia --project=. run_plate.jl
```

## 4. 使用 MLFMA 加速（大规模问题）

对于 N > 10,000 个未知量，推荐使用 MLFMA 替代直接装配：

```julia
# 替换步骤 4–6
op       = EFIE(freq)
mlfma_op = MLFMAOperator(op, basis)

V        = excitation_vector(inc_wave, basis)
precon   = BlockJacobiPreconditioner(mlfma_op, basis)
I_coeff  = solve!(GMRESSolver(tol = 1e-4, maxiter = 500, restart = 50),
                  mlfma_op, V, Pl = precon)
```

## 5. 可视化电流分布（ParaView）

```julia
save_vtk("plate_currents", mesh, abs.(I_coeff))
```

生成 `plate_currents.vtu`，可在 ParaView 中打开查看表面电流密度。

