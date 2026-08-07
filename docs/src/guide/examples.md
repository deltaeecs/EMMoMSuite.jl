# 示例

## 1. PEC 球面散射与 Mie 级数对照

PEC 球面散射是最基础的验证算例，可与 Mie 级数精确解进行比较。

```julia
using EMMoMSuite

freq = 1e9
radius = 0.1
set_frequency!(freq)

mesh = generate_sphere_mesh(radius, 8, 16)
basis = RWGBasis(mesh)
println("N = $(num_basis(basis)) DOF")

cfie = CFIE(freq, 0.5)
Z = assemble_impedance_matrix(cfie, basis)

inc = PlaneWave(freq, π / 2, π, [1.0, 0.0, 0.0])
V = excitation_vector(inc, basis)
I = solve!(LUSolver(), Z, V)

theta = collect(range(0.0, π, length = 181))
_, _, rcs_db = radarCrossSection(theta, [0.0], I, basis)
println("Monostatic RCS at 90 deg = $(rcs_db[91, 1]) dBsm")
```

## 2. 金属目标散射与 MLFMA

大规模金属散射问题通常需要 MLFMA 加速矩阵向量乘。

```julia
using EMMoMSuite

freq = 3e9
set_frequency!(freq)

mesh = read_nas_mesh("jet.nas")
basis = RWGBasis(mesh)
println("N = $(num_basis(basis)) RWG basis functions")

efie = EFIE(freq)
mlfma_op = MLFMAOperator(efie, basis, leafCubeEdgel = 0.3)

inc = PlaneWave(freq, π / 2, π, [1.0, 0.0, 0.0])
V = excitation_vector(inc, basis)
precon = BlockJacobiPreconditioner(mlfma_op, basis)
I_coeff = solve!(
    GMRESSolver(tol = 1e-4, maxiter = 500, restart = 50),
    mlfma_op,
    V,
    Pl = precon,
)

theta = collect(range(0.0, π, length = 181))
_, _, rcs_db = radarCrossSection(theta, [0.0], I_coeff, basis)
save_RCS_txt("jet_rcs.txt", theta, [0.0], rcs_db)
```

## 3. 介质体散射与 VEFIE

```julia
using EMMoMSuite

freq = 1e9
set_frequency!(freq)

mesh = read_nas_mesh("dielectric_sphere.nas")
eps_r = [ComplexF64(4.0, -0.1)]
basis = SWGBasis(mesh)

vefie = VEFIE(freq)
Z = assemble_impedance_matrix(vefie, basis, eps_r)
inc = PlaneWave(freq, π / 2, π, [1.0, 0.0, 0.0])
V = excitation_vector(inc, basis, eps_r)
I = solve!(LUSolver(), Z, V)

theta = collect(range(0.0, π, length = 181))
_, _, rcs_db = radarCrossSection(theta, [0.0], I, basis, eps_r)
save_RCS_txt("dielectric_rcs.txt", theta, [0.0], rcs_db)
```

## 4. MPI 分布式并行

```bash
mpiexecjl -n 4 julia --project=. run_parallel.jl
```

```julia
# run_parallel.jl
using EMMoMSuite, MPI

MPI.Init()
init_parallel!()

freq = 1e9
set_frequency!(freq)

mesh = read_nas_mesh("large_plate.nas")
basis = RWGBasis(mesh)
efie = EFIE(freq)

Z_mpi = assemble_impedance_matrix_parallel(efie, basis)
inc = PlaneWave(freq, π / 2, π, [1.0, 0.0, 0.0])
V = excitation_vector(inc, basis)
I = mpi_gmres(Z_mpi, V; tol = 1e-4, maxiter = 300)

if mpi_rank() == 0
    theta = collect(range(0.0, π, length = 181))
    _, _, rcs_db = radarCrossSection(theta, [0.0], I, basis)
    save_RCS_txt("parallel_rcs.txt", theta, [0.0], rcs_db)
end

MPI.Finalize()
```

