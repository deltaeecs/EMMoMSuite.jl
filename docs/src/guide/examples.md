# 示例 (Examples)

## 1. PEC 球面 Mie 级数验证（CFIE）

PEC 球面散射是最基础的验证算例，可与 Mie 级数精确解对比。

```julia
using EMSuite

freq  = 1e9                          # 1 GHz
r_m   = 0.1                          # 球半径 0.1 m ≈ 0.33λ
set_frequency!(freq)

mesh  = generate_sphere_mesh(r_m, 8) # 生成球面三角网格
basis = RWGBasis(mesh)
println("N = $(num_basis(basis)) DOF")

cfie = CFIE(freq, 0.5)               # α=0.5 避免内部谐振
Z    = assemble_impedance_matrix(cfie, basis)

inc  = PlaneWave(freq, π/2, π, [1.0, 0.0, 0.0])   # +z 方向入射，x 极化
V    = excitation_vector(inc, basis)
I    = solve!(LUSolver(), Z, V)

θ = collect(range(0.0, π, length = 181))
_, _, rcs_db = radarCrossSection(θ, [0.0], I, basis)
println("单站 RCS (θ=90°): $(rcs_db[91,1]) dBsm")
```

## 2. 金属目标散射（EFIE + MLFMA）

大规模问题使用 MLFMA 加速矩阵向量乘。

```julia
using EMSuite

freq = 3e9
set_frequency!(freq)

mesh  = read_nas_mesh("jet.nas")
basis = RWGBasis(mesh)
println("N = $(num_basis(basis)) RWG 基函数")

efie     = EFIE(freq)
mlfma_op = MLFMAOperator(efie, basis, leafCubeEdgel = 0.3)

inc     = PlaneWave(freq, π/2, π, [1.0, 0.0, 0.0])
V       = excitation_vector(inc, basis)
precon  = BlockJacobiPreconditioner(mlfma_op, basis)
I_coeff = solve!(GMRESSolver(tol = 1e-4, maxiter = 500, restart = 50),
                 mlfma_op, V, Pl = precon)

θ = collect(range(0.0, π, length = 181))
_, _, rcs_db = radarCrossSection(θ, [0.0], I_coeff, basis)
save_RCS_txt("jet_rcs.txt", θ, [0.0], rcs_db)
```

## 3. 介质体散射（VEFIE + SWG）

```julia
using EMSuite

freq = 1e9
set_frequency!(freq)

mesh  = read_nas_mesh("dielectric_sphere.nas")   # 四面体网格
εr    = [ComplexF64(4.0, -0.1)]
basis = SWGBasis(mesh)

vefie = VEFIE(freq)
Z     = assemble_impedance_matrix(vefie, basis, εr)
inc   = PlaneWave(freq, π/2, π, [1.0, 0.0, 0.0])
V     = excitation_vector(inc, basis, εr)
I     = solve!(LUSolver(), Z, V)

θ = collect(range(0.0, π, length = 181))
_, _, rcs_db = radarCrossSection(θ, [0.0], I, basis, εr)
save_RCS_txt("dielectric_rcs.txt", θ, [0.0], rcs_db)
```

## 4. MPI 分布式并行

```bash
# 启动 4 进程
mpiexecjl -n 4 julia --project=. run_parallel.jl
```

```julia
# run_parallel.jl
using EMSuite, MPI

MPI.Init()
init_parallel!()

freq = 1e9;  set_frequency!(freq)
mesh  = read_nas_mesh("large_plate.nas")
basis = RWGBasis(mesh)
efie  = EFIE(freq)

Z_mpi = assemble_impedance_matrix_parallel(efie, basis)  # MPIMatrix
inc   = PlaneWave(freq, π/2, π, [1.0, 0.0, 0.0])
V     = excitation_vector(inc, basis)
I     = mpi_gmres(Z_mpi, V; tol = 1e-4, maxiter = 300)

if mpi_rank() == 0
    θ = collect(range(0.0, π, length = 181))
    _, _, rcs_db = radarCrossSection(θ, [0.0], I, basis)
    save_RCS_txt("parallel_rcs.txt", θ, [0.0], rcs_db)
end
MPI.Finalize()
```

