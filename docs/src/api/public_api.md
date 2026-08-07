# API 参考

本页面从源代码 docstring 自动生成，覆盖 EMMoMSuite 的公开导出类型与函数。

模块分类速查：

| 类别 | 关键类型/函数 |
|------|---------------|
| **几何** | `TriangleMesh`, `TetrahedraMesh`, `HexahedraMesh`, `read_nas_mesh`, `generate_sphere_mesh` |
| **基函数** | `RWGBasis`, `SWGBasis`, `PWCBasis`, `PWCHexBasis`, `RBFBasis` |
| **积分方程** | `EFIE`, `MFIE`, `CFIE`, `VEFIE`, `SCFIE`, `assemble_impedance_matrix` |
| **MLFMA** | `MLFMAOperator`, `MLFMAOperatorMPI`, `get_leaf_intervals` |
| **求解器** | `LUSolver`, `GMRESSolver`, `BiCGSTABSolver`, `BlockJacobiPreconditioner` |
| **激励源** | `PlaneWave`, `DeltaGapSource`, `excitation_vector` |
| **后处理** | `radarCrossSection`, `farField`, `calculate_near_field`, `save_vtk` |
| **并行** | `init_parallel!`, `MPIMatrix`, `mpi_gmres!`, `mpi_rank` |

---

## 完整 API 文档

```@autodocs
Modules = [EMMoMSuite]
```

---

## 符号索引

```@index
```
