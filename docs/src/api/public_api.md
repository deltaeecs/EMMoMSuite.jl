# API 参考

本页面从源代码 docstring 自动生成，覆盖 EMMoMSuite 及其全部子模块的公开导出类型与函数。

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
Modules = [EMMoMSuite, EMMoMSuite.Accuracy, EMMoMSuite.Accuracy.AccuracyMetrics, EMMoMSuite.Accuracy.FekoReader, EMMoMSuite.Accuracy.ReferenceData, EMMoMSuite.BasisFunctions, EMMoMSuite.CoreModule, EMMoMSuite.CoreModule.Configuration, EMMoMSuite.CoreModule.Constants, EMMoMSuite.CoreModule.Materials, EMMoMSuite.CoreModule.Results, EMMoMSuite.CoreModule.Sources, EMMoMSuite.Driver, EMMoMSuite.FastAlgorithms, EMMoMSuite.FastAlgorithms.Lebedev, EMMoMSuite.FastAlgorithms.Lebedev.LVI, EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints, EMMoMSuite.FastAlgorithms.Lebedev.SHInterp, EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator, EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW, EMMoMSuite.FastAlgorithms.MLFMA, EMMoMSuite.FastAlgorithms.MLFMA.Aggregation, EMMoMSuite.FastAlgorithms.MLFMA.Disaggregation, EMMoMSuite.FastAlgorithms.MLFMA.Interpolation, EMMoMSuite.FastAlgorithms.MLFMA.Level, EMMoMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule, EMMoMSuite.FastAlgorithms.MLFMA.Octree, EMMoMSuite.FastAlgorithms.MLFMA.OctreeBuilder, EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule, EMMoMSuite.FastAlgorithms.MLFMA.Precomputations, EMMoMSuite.FastAlgorithms.MLFMA.Translation, EMMoMSuite.Geometry, EMMoMSuite.IO, EMMoMSuite.IO.ResultIO, EMMoMSuite.IO.VTKExport, EMMoMSuite.IntegralEquations, EMMoMSuite.IntegralEquations.CFIEModule, EMMoMSuite.IntegralEquations.EFIEModule, EMMoMSuite.IntegralEquations.EFIEModule.Singularities, EMMoMSuite.IntegralEquations.Excitation, EMMoMSuite.IntegralEquations.Impedance, EMMoMSuite.IntegralEquations.Kernels, EMMoMSuite.IntegralEquations.MFIEModule, EMMoMSuite.IntegralEquations.NMullerModule, EMMoMSuite.IntegralEquations.PMCHWBlockOperatorsModule, EMMoMSuite.IntegralEquations.PMCHWModule, EMMoMSuite.IntegralEquations.SCFIEModule, EMMoMSuite.IntegralEquations.VEFIEModule, EMMoMSuite.IntegralEquations.VEFIEModule.FastExpModule, EMMoMSuite.MaterialsModule, EMMoMSuite.Parallel, EMMoMSuite.Parallel.Assembly, EMMoMSuite.Parallel.MPIArrays, EMMoMSuite.PortsModule, EMMoMSuite.PostProcessing, EMMoMSuite.PostProcessing.Absorption, EMMoMSuite.PostProcessing.AntennaMetrics, EMMoMSuite.PostProcessing.CurrentOnGeos, EMMoMSuite.PostProcessing.FarField, EMMoMSuite.PostProcessing.FarFieldPatternModule, EMMoMSuite.PostProcessing.FieldCut, EMMoMSuite.PostProcessing.MLFMAFastPostModule, EMMoMSuite.PostProcessing.NearField, EMMoMSuite.PostProcessing.NearFieldAdvancedModule, EMMoMSuite.PostProcessing.RCS, EMMoMSuite.PostProcessing.RCSBatchModule, EMMoMSuite.PostProcessing.RadiationIntegral, EMMoMSuite.PostProcessing.SolverResultModule, EMMoMSuite.Solvers, EMMoMSuite.Solvers.DirectSolvers, EMMoMSuite.Solvers.IterativeWrappers, EMMoMSuite.Solvers.Preconditioners, EMMoMSuite.Utilities, EMMoMSuite.Utilities.HDF5Utils, EMMoMSuite.Utilities.LightweightSupport, EMMoMSuite.Utilities.LoggingUtils, EMMoMSuite.Utilities.MieSeries, EMMoMSuite.Utilities.Parameters]
Public = true
Private = false
```

---

## 符号索引

```@index
```
