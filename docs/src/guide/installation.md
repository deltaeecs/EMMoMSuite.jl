# 安装指南 (Installation Guide)

## 系统要求

*   **Julia**: 1.10 或更高（推荐 1.12+）。从 [julialang.org](https://julialang.org/downloads/) 下载。
*   **MPI**（可选）：并行计算功能需要 MPI 库（如 MPICH、OpenMPI 或 Windows 上的 Microsoft MPI）。

## 安装方式

### 方式 1：本地开发（推荐）

1. 克隆仓库：
    ```bash
    git clone <your-repo-url>/EMSuite.jl.git
    cd EMSuite.jl
    ```

2. 在项目目录启动 Julia：
    ```bash
    julia --project=. --startup-file=no
    ```

3. 安装依赖：
    ```julia
    using Pkg
    Pkg.instantiate()
    ```

### 方式 2：以开发包形式添加

```julia
using Pkg
Pkg.develop(path = "/path/to/EMSuite.jl")
```

## 验证安装

在 Julia REPL 中运行：

```julia
using EMSuite
mesh = generate_sphere_mesh(1.0, 5)      # 生成球面三角网格
basis = RWGBasis(mesh)
println("RWG DOF: ", num_basis(basis))   # 应输出 > 0
```

或运行完整测试套件：

```bash
julia --project=EMSuite --startup-file=no EMSuite/run_tests.jl
```

## MPI 并行支持

EMSuite 支持通过 MPI.jl 进行分布式并行。安装 MPI 后初始化：

```julia
using Pkg
Pkg.add("MPI")

using MPI
MPI.install_mpiexecjl()  # 安装 mpiexecjl 启动器到 ~/.julia/bin/
```

启动 4 进程并行仿真：

```bash
mpiexecjl -n 4 julia --project=. my_simulation.jl
```

