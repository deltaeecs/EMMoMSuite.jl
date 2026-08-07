# Quick Start

本指南演示如何使用 EMMoMSuite 计算 PEC 金属板的双站 RCS，流程覆盖网格读取、积分方程装配、求解和结果导出。

## 1. 准备网格

EMMoMSuite 支持 Nastran `.nas` 和 Gmsh `.msh` 等常见网格格式。下面示例假设当前目录中已有一块位于 XY 平面的矩形金属板网格 `plate.nas`。

## 2. 编写求解脚本

```julia
using EMMoMSuite

freq = 300e6
set_frequency!(freq)

mesh = read_nas_mesh("plate.nas")
println("网格: $(num_elements(mesh)) 个三角形, $(num_vertices(mesh)) 个节点")

basis = RWGBasis(mesh)
println("未知量数: $(num_basis(basis))")

efie = EFIE(freq)
Z = assemble_impedance_matrix(efie, basis)

inc_wave = PlaneWave(freq, π / 2, π, [1.0, 0.0, 0.0])
V = excitation_vector(inc_wave, basis)

I_coeff = solve!(LUSolver(), Z, V)

theta_obs = collect(range(0.0, π, length = 181))
phi_obs = [0.0]

_, rcs_total, rcs_db = radarCrossSection(theta_obs, phi_obs, I_coeff, basis)
println("单站 RCS (theta = 90 deg): ", rcs_db[91, 1], " dBsm")

save_RCS_txt("plate_rcs.txt", theta_obs, phi_obs, rcs_db)
save_results_hdf5("plate_results.h5", I_coeff, rcs_db)
println("完成，结果已保存到 plate_rcs.txt 和 plate_results.h5")
```

## 3. 运行脚本

```bash
julia --project=. run_plate.jl
```

## 4. 使用 MLFMA 加速大规模问题

当未知量规模较大时，可以把直接矩阵装配替换为 MLFMA 算子与迭代求解。

```julia
op = EFIE(freq)
mlfma_op = MLFMAOperator(op, basis)

V = excitation_vector(inc_wave, basis)
precon = BlockJacobiPreconditioner(mlfma_op, basis)
I_coeff = solve!(
    GMRESSolver(tol = 1e-4, maxiter = 500, restart = 50),
    mlfma_op,
    V,
    Pl = precon,
)
```

## 5. 导出可视化结果

```julia
save_vtk("plate_currents", mesh, abs.(I_coeff))
```

生成的 `plate_currents.vtu` 可以在 ParaView 中查看表面电流分布。

