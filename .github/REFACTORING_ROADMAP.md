# EMSuite 重构路线图

> 最后更新: 2026-03-01

## 项目概述

将分散的 Legacy 代码 (`MoM_Basics`, `MoM_Kernels`, `MoM_MPI`, `MoM_Lebedev`, `MoM_AllinOne`, `MPIArray4MoMs`) 重构为统一的 `EMSuite.jl` 包。

## 架构

```
EMSuite.jl/src/
├── Core/            # 抽象接口、类型、常数、材料、激励源、配置
├── Geometry/        # 网格类型、I/O (Nastran/Gmsh)、坐标变换、高斯求积
├── BasisFunctions/  # RWG, SWG, RBF, PWC
├── IntegralEquations/ # EFIE, MFIE, CFIE, VEFIE, SCFIE, 奇异性处理
├── FastAlgorithms/  # MLFMA (Octree, Agg, Trans, Disagg), Lebedev
├── Solvers/         # Direct (LU), Iterative (GMRES, BiCGSTAB), 预条件器
├── PostProcessing/  # RCS, FarField, NearField, 电流分布
├── Parallel/        # MPI, Threading
├── IO/              # VTK, HDF5, CSV, 结果文件
├── Utilities/       # 日志、参数、Mie 级数
└── Driver.jl        # 统一仿真入口
```

---

## 阶段划分

### Phase 1: 基础架构 ✅
- 项目结构、Core 模块、CI/CD、文档框架

### Phase 2: 几何与基函数 ✅
- 网格类型 (Triangle, Tetrahedra, Hexahedra)
- Nastran/Gmsh I/O、坐标变换、高斯求积
- RWG, SWG, RBF, PWC 基函数

### Phase 3: 积分方程与矩阵填充 ✅
- EFIE, MFIE, CFIE (表面)
- VEFIE (体积)
- SCFIE (面体耦合)
- 奇异性处理 (解析提取)

### Phase 4: MLFMA 快速算法 ✅
- 八叉树、聚合、转移、解聚
- Lebedev 球面插值
- 近场/远场精度验证

### Phase 5: 求解器与并行 ✅
- Direct (LU), GMRES, BiCGSTAB
- ILU, SPAI 预条件器
- MPI 分布式 + 多线程并行

### Phase 6: 后处理与 I/O ✅
- RCS, FarField, NearField
- VTK 导出、结果文件 I/O

### Phase 7: 验证与对齐 ✅
- [x] **SCFIE MLFMA**: `MLFMAOperator` 支持混合基 (RWG+SWG) — 近场 1.58e-15, 远场 0.85%
- [ ] **MoM_AllinOne 对标**: 完成全部 8 个算例对比
  - [x] SEFIE Direct / MLFMA
  - [x] VEFIE Direct / MLFMA
  - [x] SCFIE Direct (examples 对标)
  - [x] SCFIE MLFMA (examples 对标) — 7 bugs fixed, near-field 1.58e-15, far-field 0.85%
- [x] **回归测试锁定**: 138 tests 全部通过, 覆盖核心路径

### Phase 8: 性能优化 (延期)
- [ ] Signed Geometric Properties 优化 (消除内循环分支) — 收益较小，风险较高，延期
- [ ] 矩阵组装 SIMD/线程优化
- [ ] MLFMA 大规模扩展性

### Phase 9: 代码质量与发布
- [x] 测试套件清理: 修复所有预存测试失败, 138/138 全部通过
- [ ] JuliaFormatter 统一代码风格
- [ ] 测试覆盖率 > 80%
- [ ] API 文档完善 (Documenter.jl)
- [ ] 用户教程和理论文档
- [ ] 发布到 Julia General Registry

### Phase 10: 全方程全路径精度对齐 (**当前**)

> 目标: 对 5 类积分方程 × 4 种求解路径进行全球面 RCS 精密对比，定量证明 EMSuite 与 Legacy 一致。

#### 10.0 历史验证结果 (Phase 7)

| 算子 | 求解器 | 对标基准 | 容差 | 状态 |
|------|--------|----------|------|------|
| S-EFIE | Direct | Legacy | < 1 dB | ✅ RMSE 0.29 dB (φ=0), 0.09 dB (φ=90) |
| S-EFIE | MLFMA | Direct | Rel Err < 2% | ✅ MatVec 0.84% |
| S-CFIE | Direct | EFIE Direct | 主瓣 < 1 dB | ✅ 主瓣差 0.34 dB |
| S-CFIE | MLFMA | Legacy Direct | < 2 dB | ⚠️ φ=0 RMSE 2.30 dB, φ=90 1.26 dB |
| VEFIE | MLFMA | Direct | Rel Err < 2% | ✅ 1.80% |
| VSIE (SCFIE) | MLFMA | Direct | Near<1e-12, Far<2% | ✅ 0%/0.75% |

> **P0 已修复** (2026-02-28): `edgev̂` 方向反转 + `calc_near_interaction!` 面积归一化。

#### 10.1 全球面采样方案

替换旧的 "φ=0 + φ=π/2 两条切面" 方案，采用:

| 参数 | 值 | 说明 |
|------|----|----- |
| θ 范围 | [-π, π] | 等同双站 RCS 全角度扫描 |
| θ 采样 | 73 点 (5° 间隔) | 兼容 Legacy 721 点子集 |
| φ 范围 | [0, π) | 半球对称, 避免冗余 |
| φ 采样 | 18 条切面 (10° 间隔) | 包含 φ=0°/90° |
| 总观测方向 | 73 × 18 = **1314** | 覆盖全球面 |

输出格式: CSV 文件, 列 = `theta_deg, phi_deg, RCS_theta_dB, RCS_phi_dB, RCS_total_dB`

#### 10.2 测试矩阵

| 编号 | 方程类型 | EMSuite 算子 | 几何体 | 网格文件 | 频率 | 基函数 | N (approx) | Legacy ieT |
|------|----------|-------------|--------|----------|------|--------|------------|-----------|
| **A** | S-EFIE | `EFIE(freq)` | Jet (开体) | `jet_100MHz.nas` | 100 MHz | RWG | 14559 | `:EFIE` |
| **B** | S-MFIE | `MFIE(freq)` | Sphere (闭体) | `sphere_600MHz.nas` | 600 MHz | RWG | 26424 | N/A¹ |
| **C** | S-CFIE | `CFIE(freq, 0.5)` | Sphere (闭体) | `sphere_600MHz.nas` | 600 MHz | RWG | 26424 | `:CFIE` |
| **D** | V-EFIE | `VEFIE(freq, perms)` | 介质板 | `plate_1dot2GHz.nas` | 1.2 GHz | SWG | ~6000 | `:EFIE, vbf=SWG` |
| **E** | VS-EFIE | `SCFIE(freq, perms)` | PEC+介质 | `plate_and_metal_1dot2GHz.nas` | 1.2 GHz | RWG+SWG | ~6500 | `:EFIE, s+vbf` |

¹ Legacy 无独立 MFIE 示例, 仅做自洽验证与 CFIE 分解验证。

#### 10.3 求解器路径

| 路径 | 说明 | 可用方程 | 约束 |
|------|------|---------|------|
| **Direct** | Dense Z → LU | A, D, E | B/C 的 N=26424 Dense Z 需 11 GB, 不可行 |
| **Iterative** | Dense Z → GMRES (预条件=对角) | A, D, E | 同上, 排除 B/C |
| **MLFMA** | MLFMAOperator → GMRES + SAI | A, B, C, D, E | 全部可用 |
| **MPI** | `assemble_impedance_matrix_parallel` | A, C | 仅支持 RWG 表面方程; D/E (SWG) 无 MPI |

#### 10.4 各测试用例详细方案

##### A: S-EFIE — Jet 100 MHz (开体, N=14559)

| 子项 | 求解路径 | 对比基准 | 通过准则 | 状态 |
|------|---------|---------|---------|------|
| A1 | Direct + LU | Legacy Direct (全球面) | RMSE < 1 dB | [ ] |
| A2 | Iterative (GMRES on Dense Z) | A1 | RMSE < 0.1 dB | [ ] |
| A3 | MLFMA + GMRES | Legacy MLFMA, A1 | vs Legacy < 3 dB, vs A1 < 2 dB | [ ] |
| A4 | MPI (2 进程) | A1 | 机器精度 (< 1e-10 dB) | [ ] |

##### B: S-MFIE — Sphere 600 MHz (闭体, N=26424)

| 子项 | 验证方式 | 通过准则 | 状态 |
|------|---------|---------|------|
| B1 | CFIE 分解验证² | ‖Z_CFIE - α·Z_EFIE - (1-α)·Z_MFIE‖/‖Z_CFIE‖ < 1e-12 | [ ] |
| B2 | MLFMA + GMRES → RCS | RCS 趋势与 S-CFIE MLFMA 物理一致 | [ ] |
| B3 | MPI (2 进程) | vs B2 机器精度 | [ ] |

² 在小网格 (`Tri.nas`, N~1330) 上组装 Z_EFIE、Z_MFIE、Z_CFIE，验证线性组合关系。

##### C: S-CFIE — Sphere 600 MHz (闭体, N=26424)

| 子项 | 求解路径 | 对比基准 | 通过准则 | 状态 |
|------|---------|---------|---------|------|
| C1 | MLFMA + GMRES | Legacy MLFMA (全球面) | RMSE < 2 dB | [ ] |
| C2 | MLFMA + GMRES | Legacy Direct (全球面) | RMSE < 3 dB (含 Float32/64 差异) | [ ] |
| C3 | MPI (2 进程) | C1 | 机器精度 | [ ] |

##### D: V-EFIE — 介质板 1.2 GHz (N~6000)

| 子项 | 求解路径 | 对比基准 | 通过准则 | 状态 |
|------|---------|---------|---------|------|
| D1 | Direct + LU | Legacy Direct (全球面) | RMSE < 1 dB | [ ] |
| D2 | Iterative (GMRES on Dense Z) | D1 | RMSE < 0.1 dB | [ ] |
| D3 | MLFMA + GMRES | Legacy MLFMA, D1 | vs Legacy < 3 dB, vs D1 < 2 dB | [ ] |

##### E: VS-EFIE — PEC+介质 1.2 GHz (N~6500)

| 子项 | 求解路径 | 对比基准 | 通过准则 | 状态 |
|------|---------|---------|---------|------|
| E1 | Direct + LU | Legacy Direct (全球面) | RMSE < 1 dB | [ ] |
| E2 | Iterative (GMRES on Dense Z) | E1 | RMSE < 0.1 dB | [ ] |
| E3 | MLFMA + GMRES | Legacy MLFMA, E1 | vs Legacy < 3 dB, vs E1 < 2 dB | [ ] |

#### 10.5 传递准则汇总

| 对比类型 | 准则 | 说明 |
|---------|------|------|
| EMSuite Direct vs Legacy Direct | RMSE < 1 dB (全球面) | 核心指标, Float64 vs Float32 |
| EMSuite Iterative vs EMSuite Direct | RMSE < 0.1 dB | 收敛精度决定 |
| EMSuite MLFMA vs EMSuite Direct | RMSE < 2 dB | MLFMA 截断误差 |
| EMSuite MLFMA vs Legacy MLFMA | RMSE < 3 dB | 允许精度差+截断差 |
| EMSuite MPI vs EMSuite Serial | 机器精度 | 并行分块正确性 |
| S-CFIE Z 分解 | ‖Z_CFIE - α·Z_E - (1-α)·Z_M‖_F / ‖Z_CFIE‖_F < 1e-12 | MFIE 算子正确性 |

#### 10.6 实施步骤

| 步骤 | 内容 | 产出文件 | 状态 |
|------|------|---------|------|
| Step 1 | 全球面 Legacy 基线生成 (8 用例) | `test_results/legacy_baseline/fullsphere/*.csv` | [ ] |
| Step 2 | S-MFIE CFIE 分解验证 (B1) | `test/test_mfie_decomposition.jl` | [ ] |
| Step 3 | EMSuite 全方程 Direct 基准 (A1, D1, E1) | `benchmark/precision_alignment_direct.jl` | [ ] |
| Step 4 | EMSuite Iterative 基准 (A2, D2, E2) | 同上脚本 Iterative 模式 | [ ] |
| Step 5 | EMSuite MLFMA 基准 (A3, B2, C1, D3, E3) | `benchmark/precision_alignment_mlfma.jl` | [ ] |
| Step 6 | MPI 基准 (A4, B3, C3) | `benchmark/precision_alignment_mpi.jl` | [ ] |
| Step 7 | 生成全球面误差热力图 + 更新报告 | `test_results/ACCURACY_EFFICIENCY_REPORT.md` v3 | [ ] |

#### 10.7 EMSuite API 覆盖矩阵

| 功能 | EFIE | MFIE | CFIE | VEFIE | SCFIE |
|------|------|------|------|-------|-------|
| 直接装配 | ✅ RWG | ✅ RWG | ✅ RWG | ✅ SWG | ✅ RWG+SWG |
| 激励向量 | ✅ | ✅ | ✅ | ✅ | ✅ (拼接) |
| MLFMAOperator | ✅ | ✅ | ✅ | ✅ | ✅ |
| MPI 并行装配 | ✅ | ✅ | ✅ | ❌ | ❌ |
| RCS 计算 | ✅ RWG | ✅ RWG | ✅ RWG | ✅ SWG | ⚠️ 需手动拆分 |
| LU 求解 | ✅ | ✅ | ✅ | ✅ | ✅ |
| GMRES 求解 | ✅ | ✅ | ✅ | ✅ | ✅ |

---

## 关键参考

- **Legacy 代码**: `MoM_Basics/`, `MoM_Kernels/`, `MoM_AllinOne/`
- **验证脚本**: `EMSuite/benchmark/verify_*.jl`, `EMSuite/scripts/verification/`
- **理论**: Harrington "Field Computation by Moment Methods"; Chew et al. "Fast and Efficient Algorithms in CEM"
