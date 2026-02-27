# EMSuite.jl

**EMSuite.jl** 是一个基于 Julia 语言开发的综合性计算电磁学（Computational Electromagnetics, CEM）矩量法（Method of Moments, MoM）框架。旨在提供高效、灵活且易于扩展的电磁仿真工具。

## 主要特性 (Features)

- **积分方程 (Integral Equations)**
  - 电场积分方程 (EFIE)
  - 磁场积分方程 (MFIE)
  - 混合场积分方程 (CFIE)
  
- **基函数 (Basis Functions)**
  - RWG (Rao-Wilton-Glisson) 基函数
  - SWG (Schaubert-Wilton-Glisson) 基函数
  - RBF (Rooftop Basis Functions) 屋顶基函数
  - PWC (Piecewise Constant) 脉冲基函数

- **快速算法 (Fast Algorithms)**
  - 多层快速多极子算法 (MLFMA)
  - Lebedev 正交积分

- **求解器 (Solvers)**
  - 直接求解器 (Direct Solvers): LU 分解
  - 迭代求解器 (Iterative Solvers): GMRES, BiCGSTAB (集成 IterativeSolvers.jl)

- **并行计算 (Parallel Computing)**
  - 支持 MPI 分布式内存并行
  - 支持 Julia 原生多线程 (Multi-threading) 并行

- **后处理与 I/O (Post-Processing & I/O)**
  - 雷达散射截面 (RCS) 计算
  - 远场/近场计算
  - 结果导出: 
    - VTK (.vtu) 格式 (支持 ParaView/Tecplot 可视化)
    - HDF5 格式
    - TXT 文本格式

## 安装 (Installation)

EMSuite 是一个 Julia 包。在 Julia REPL 中按 `]` 进入包管理模式，然后运行：

```julia
pkg> add EMSuite
```

或者使用 `Pkg` 模块：

```julia
using Pkg
Pkg.add("EMSuite")
```

## 快速开始 (Quick Start)

以下是一个计算简单几何体散射问题的完整示例：

```julia
using EMSuite
using LinearAlgebra

# 1. 设置仿真参数
freq = 300e6 # 频率 300 MHz
set_frequency!(freq)

# 2. 创建或加载网格
# mesh = read_nas_mesh("plate.nas") 
# 或者手动创建简单的几何 (例如一个三角形):
nodes = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]' # 3x3 节点矩阵
elements = [1 2 3]' # 3x1 单元矩阵
mesh = TriangleMesh(1, nodes, elements)

# 3. 定义基函数 (RWG)
basis = RWGBasis(mesh)

# 4. 定义积分方程 (EFIE)
ie = EFIE(freq)

# 5. 组装阻抗矩阵
println("正在组装阻抗矩阵...")
Z = assemble_impedance_matrix(ie, basis)

# 6. 定义激励 (例如平面波)
# 这里仅作为示例使用全1向量，实际应用中需计算入射场向量 V_m = <E_inc, f_m>
n_unknowns = num_basis(basis)
V = ones(ComplexF64, n_unknowns)

# 7. 求解线性方程组
println("正在求解...")
I = solve!(LUSolver(), Z, V)

# 8. 后处理与结果导出
# 导出电流分布到 VTK 文件 (可使用 ParaView 查看)
save_vtk("simulation_result", mesh, vec(abs.(I)); data_name="CurrentMagnitude")

println("计算完成！结果已保存为 simulation_result.vtu")
```

## 模块结构 (Module Structure)

EMSuite 采用模块化设计，主要包含以下组件：

*   **Core**: 定义核心抽象类型 (`AbstractMesh`, `AbstractBasisFunction` 等) 和接口。
*   **Geometry**: 网格处理 (`TriangleMesh`, `TetrahedraMesh`)、几何变换和网格文件读取 (`.nas`, `.msh`)。
*   **BasisFunctions**: 各种基函数的定义与实现 (`RWG`, `SWG` 等)。
*   **IntegralEquations**: 积分方程算子 (`EFIE`, `MFIE`, `CFIE`) 与矩阵填充逻辑。
*   **Solvers**: 线性方程组求解器接口与封装。
*   **FastAlgorithms**: MLFMA 及其相关组件 (Octree, Aggregation, Translation 等)。
*   **PostProcessing**: RCS、近远场计算等物理量提取。
*   **IO**: 结果数据的导入导出 (`VTKExport`, `ResultIO`)。
*   **Utilities**: 日志记录、物理常数管理等辅助工具。

## 文档 (Documentation)

更多详细信息、API 文档和高级用法，请参阅 [在线文档](https://deltaeecs.github.io/EMSuite.jl/dev)。

## 许可证 (License)

本项目采用 MIT 许可证。
