# 安装指南

## 系统要求

- Julia 1.10 或更高版本，推荐 1.12+。
- MPI 运行时是可选依赖，仅在分布式并行场景中需要，例如 MPICH、OpenMPI 或 Windows 上的 Microsoft MPI。

## 方式 1：本地开发环境

1. 克隆仓库。

```bash
git clone <your-repo-url>/EMMoMSuite.jl.git
cd EMMoMSuite.jl
```

2. 进入项目环境并安装依赖。

```bash
julia --project=. --startup-file=no
```

```julia
using Pkg
Pkg.instantiate()
```

## 方式 2：以开发包形式接入现有环境

```julia
using Pkg
Pkg.develop(path = "/path/to/EMMoMSuite.jl")
```

## 验证安装

在 Julia REPL 中执行以下代码，确认主包可导入且基础网格与基函数构造正常。

```julia
using EMMoMSuite

mesh = generate_sphere_mesh(1.0, 5)
basis = RWGBasis(mesh)
println("RWG DOF: ", num_basis(basis))
```

若输出的自由度数量大于 0，则说明核心环境已经可用。

也可以直接运行项目测试入口。

```bash
julia --project=. --startup-file=no run_tests.jl
```

## MPI 并行支持

若需要 MPI 并行，先安装 MPI.jl 并配置启动器。

```julia
using Pkg
Pkg.add("MPI")

using MPI
MPI.install_mpiexecjl()
```

之后可按如下方式启动 4 进程任务。

```bash
mpiexecjl -n 4 julia --project=. my_simulation.jl
```

