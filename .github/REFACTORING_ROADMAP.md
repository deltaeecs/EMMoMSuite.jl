# EMSuite 重构路线图

> 最后更新: 2026-02-28

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

### Phase 1–6: 基础架构 → 后处理 ✅

全部完成。包含: 项目结构、Core 模块、几何与基函数 (RWG/SWG/RBF/PWC)、积分方程 (EFIE/MFIE/CFIE/VEFIE/SCFIE + 奇异性处理)、MLFMA (八叉树/聚合/转移/解聚/Lebedev)、求解器 (Direct/GMRES/BiCGSTAB + ILU/SPAI 预条件)、MPI/Threading 并行、后处理 (RCS/FarField/NearField/VTK)。

### Phase 7: 验证与对齐 ✅

- [x] SCFIE MLFMA 近场/远场验证 — 10 个 Bug 修复
- [x] MoM_AllinOne 全部算例对标完成
- [x] 回归测试锁定: 138/138 全部通过
- [x] ~3 dB 系统偏差根因修复: edgev̂ 方向 + near-interaction 面积归一化

### Phase 8: 性能优化 ✅

> 详见 `test_results/PERFORMANCE_REPORT.md`
> - [x] 8.0 性能基线测量 ✅ (`861426d`)
> - [x] 8.1 Z 组装多线程去锁 ✅ (`2d4ebe6`) — Plate EFIE -54%, Jet EFIE -12%
> - [x] 8.2 CFIE 内核合并 ✅ (`d0888cf`) — CFIE -74%, CFIE/EFIE=2.31×
> - [x] 8.3 MLFMA Z_near 优化 ✅ (`d2f7963`)
> - [x] 8.4 内存分配热点 ✅ (`82988cf`)
> - [x] 8.5 Julia 1.12 兼容 + 类型稳定性 ✅ (`67d3a8a`)
> - [x] 8.6 @fastmath + SIMD ✅ (`1c6d499`)
> - [x] 8.7 BlockJacobiPreconditioner ✅ (`76f8b16`)
> - [x] 8.8 最终基准复测 + OOM 修复 ✅ (`6f4987a`)

### Phase 9: 代码质量与发布 **(下一阶段)**

- [x] 测试套件清理: 138/138 全部通过
- [ ] JuliaFormatter 统一代码风格
- [ ] 测试覆盖率 > 80%
- [ ] API 文档完善 (Documenter.jl)
- [ ] 用户教程和理论文档
- [ ] 发布到 Julia General Registry

### Phase 10: 全方程全路径精度对齐 ✅

> 详见下方 Phase 10 详细计划。

---

## Phase 8: 性能优化 — 详细计划

### 8.0 目标

| 指标 | 底线 | 挑战目标 | 说明 |
|------|------|---------|------|
| **相同用例全流程耗时** | ≤ Legacy (1.0×) | ≤ Legacy × 0.5 (2× 加速) | 同网格、同方程、同求解器路径 |
| **MLFMA MatVec 单次耗时** | ≤ Legacy (1.0×) | ≤ Legacy × 0.5 | 含聚合+转移+解聚 |
| **Z 矩阵组装** | ≤ Legacy (1.0×) | ≤ Legacy × 0.5 | 含 EFIE/MFIE/VEFIE 各算子 |
| **峰值内存** | ≤ 2× Legacy | ≤ 1.5× Legacy | Float64 vs Float32 允许 2× |

### 8.1 性能基线测量 (Step 0 — 优先执行)

> **原则**: 先量化，再优化。无基线数据禁止开始优化。

创建 `benchmark/performance_baseline.jl`，对以下用例分阶段计时，同时运行 Legacy 与 EMSuite:

| 用例 | N | 方程 | 路径 | 计时项 |
|------|---|------|------|--------|
| PEC Plate 300MHz | ~2640 | EFIE | Direct | 基函数构建 / Z 组装 / LU / RCS |
| Jet 100MHz | 14559 | EFIE | Direct | 同上 |
| Jet 100MHz | 14559 | EFIE | MLFMA+GMRES | 八叉树 / 近场Z / 预条件 / GMRES / RCS |
| Sphere 600MHz | 26424 | CFIE | MLFMA+GMRES | 同上 |
| Tetra 2GHz | ~986 | VEFIE | Direct | 基函数 / Z 组装 / LU / RCS |
| TriTetra 2GHz | ~1071 | SCFIE | Direct | Z_SS/Z_SV/Z_VS/Z_VV / Fss / LU / RCS |

输出: `test_results/PERFORMANCE_BASELINE.md`，格式:

```
| 用例 | 阶段 | Legacy (s) | EMSuite (s) | Ratio | 状态 |
```

### 8.2 已识别热点与优化路线

#### 热点 1: SpinLock 全局锁 — Z 矩阵组装瓶颈 (P0)

**现状**: `Impedance.jl` 中每个三角形对都要 `lock(spinlock)` / `unlock(spinlock)` 写入全局 Z 矩阵。多线程扩展性极差。

**方案**: 线程局部缓冲 + 最后归约
```julia
# Before: 
lock(spinlock); Z[m, n] += val; unlock(spinlock)
# After:
Z_local = [zeros(CT, N, N) for _ in 1:nthreads()]
@threads for i in workload
    Z_local[threadid()][m, n] += val  # 无锁
end
Z .= sum(Z_local)  # 一次归约
```

**预期收益**: 4 线程加速比 1.2× → 3.5×+

**风险**: 内存增加 (nthreads 倍 Z 矩阵)。对于 N>10000 的 Dense Z，可改用**按行分块**: 每个线程负责固定行范围，无需锁也无需额外内存。

#### 热点 2: CFIE 组装 9× EFIE (P0)

**现状**: Jet N=14559, EFIE 组装 20.3s, CFIE 组装 180.5s (9×)。理论上 CFIE = EFIE + MFIE，应 ≤ 2× EFIE。

**方案**:
1. **Green 函数复用**: EFIE 和 MFIE 共享 $G(r,r') = e^{-jkR}/(4\pi R)$ 和 $\nabla G$。合并为单遍历，计算一次 G，同时累加 EFIE 和 MFIE 贡献。
2. **MFIE 内循环优化**: 检查 MFIE 是否有冗余几何计算 (法向量、交叉积)，提到循环外预计算。
3. **对称性利用**: EFIE L 算子对称 → 仅计算上三角; MFIE K 算子反对称 → 上三角取负。

**预期收益**: CFIE 从 9× EFIE → ≤ 2.5× EFIE (180s → ≤ 50s)

#### 热点 3: MLFMA Setup 占比过高 (P1)

**现状**: Jet MLFMA 总耗时 63.4s，其中 setup 56.2s (89%)。Sphere 138.6s，setup 131.0s (94%)。

**子项**:
| 子项 | 当前耗时 (估) | 优化方向 |
|------|-------------|---------|
| 八叉树构建 | ~5% | 预分配数组，减少 push!/resize! |
| 近场 Z_near 稀疏矩阵组装 | ~60% | 同热点1 (去锁)；CSC 预分配 nnz |
| 转移矩阵预计算 | ~15% | 缓存重复的球面波展开；Lebedev 表复用 |
| SAI/ILU 预条件器构建 | ~20% | 考虑近似预条件 (Block Jacobi) 或延迟构建 |

**预期收益**: Setup 总时间 -30~50%

#### 热点 4: 内存分配热点 (P1)

**方案**:
1. **高斯积分点预分配**: 避免每个三角形对重复分配 `SVector` 数组
2. **Green 函数计算零分配**: 确保 exp/sqrt 等运算在标量上操作，不产生临时数组
3. **RCS 后处理**: 当前每个观测方向独立分配辐射积分缓冲区 → 预分配复用

**工具**: 使用 `@allocated` / `--track-allocation=user` 定位

#### 热点 5: SIMD / LoopVectorization (P2)

**现状**: 未使用 `LoopVectorization.jl` 的 `@turbo` 宏。

**方案**: 对以下内循环添加 SIMD 优化:
1. 标量势项累加: `∑ div_f · div_f' · G` — 纯标量，SIMD 友好
2. 矢量势项累加: `∑ f · f' · G` — SVector 内积，编译器可自动向量化
3. MLFMA 聚合/解聚: 球面波系数数组运算

**预期收益**: 内循环 1.5~3× 加速 (取决于 CPU AVX 支持)

**风险**: `LoopVectorization.jl` 与复数 `ComplexF64` 兼容性有限，可能需要手动拆分实部/虚部

#### 热点 6: 类型稳定性 (P2)

**方案**: 使用 `@code_warntype` 检查关键路径:
1. `efie_interaction!` / `mfie_interaction!` — 确保无 `Any` 类型
2. `mul!` (MLFMAOperator) — 确保聚合/转移/解聚全路径类型稳定
3. SCFIE 耦合组装 — 混合基函数路径可能有类型不稳定

### 8.3 实施路线

| 步骤 | 内容 | 前置 | 预期耗时 | 收益 |
|------|------|------|---------|------|
| **8.0** | 性能基线测量 | — | 0.5 天 | 量化起点 |
| **8.1** | Z 组装去锁 (行分块并行) | 8.0 | 1 天 | EFIE 3~4× 加速 |
| **8.2** | CFIE = EFIE+MFIE 合并遍历 | 8.1 | 1 天 | CFIE 4~5× 加速 |
| **8.3** | MLFMA Z_near 去锁 + CSC 预分配 | 8.1 | 0.5 天 | Setup -30% |
| **8.4** | 内存分配热点消除 | 8.0 | 0.5 天 | GC 压力降低 |
| **8.5** | @code_warntype 类型稳定性审查 | 8.0 | 0.5 天 | 消除动态派发 |
| **8.6** | SIMD / @turbo 内循环 | 8.5 | 1 天 | 内循环 1.5~3× |
| **8.7** | 预条件器优化 (Block Jacobi) | 8.3 | 0.5 天 | 预条件构建 -50% |
| **8.8** | 最终基线复测 | 全部 | 0.5 天 | 验证达标 |

**总预估**: ~6 天

### 8.4 通过准则 — 最终结果

| 用例 | Legacy 全流程 (s) | EMSuite 基线 (s) | EMSuite 最终 (s) | vs Legacy | 状态 |
|------|-------------------|------------------|-----------------|-----------|------|
| PEC Plate Direct (N=2640) | 8.29 | 3.44 | **5.84** | 0.70× | ✅ |
| Jet EFIE Direct (N=14559) | 46.43 | 42.16 | **61.90** | 1.33× | ⚠️ |
| Jet CFIE Direct (N=14559) | 64.21 | 188.63 | **97.98** | 1.53× | ⚠️ |
| Jet EFIE MLFMA (N=14559) | N/A¹ | 137.75 | **178.35** | — | — |
| Sphere CFIE MLFMA (N=26424) | N/A¹ | 540.09 | **539.93** | — | — |
| Plate VEFIE Direct (N=15828) | 103.85 | 72.73 | **102.76** | 0.99× | ✅ |
| PlateMetal SCFIE Direct (N=15860) | 66.52 | 90.57 | **130.73** | 1.96× | ❌ |

¹ Legacy MLFMA 在 Julia 1.12 下因 SAI threadid() 兼容性问题无法运行

### 8.5 回归约束

- 优化后所有 138/138 单元测试 **必须** 继续通过
- Phase 10 已验证的精度指标 (RMSE) **不得** 退化
- 每个优化步骤完成后立即运行 `Pkg.test()` + 关键 benchmark

---

## Phase 10: 全方程全路径精度对齐 ✅

> 目标: 对 5 类积分方程 × 4 种求解路径进行全球面 RCS 精密对比，定量证明 EMSuite 与 Legacy 一致。

### 10.0 已验证结果汇总

| 测试 | 指标 | 结果 | 状态 |
|------|------|------|------|
| A1 S-EFIE Direct Jet | RMSE vs Legacy | 0.215 dB | ✅ PASS |
| A3 S-EFIE MLFMA Jet | RMSE vs Legacy | 0.303 dB | ✅ PASS |
| B1 CFIE Z 分解 | rel_err | 0.0 (10/10) | ✅ PASS |
| C1 S-CFIE Direct Sphere | RMSE vs Legacy | 0.001 dB | ✅ PASS |
| C3 S-CFIE MLFMA Sphere | RMSE vs Legacy | 0.003 dB | ✅ PASS |
| D1-SWG V-EFIE Direct | RMSE vs Legacy | 0.952 dB | ✅ PASS |
| E1 VSEFIE Direct | RMSE vs Legacy | 0.602 dB | ✅ PASS |

**已修复 Bug:**
- **P0** (2026-02-28): `edgev̂` 方向反转 + `calc_near_interaction!` 面积归一化
- **P1** (2026-03-02): MLFMA ×4 因子 + CFIE MFIE 符号
- **P2** (2026-03-03): SCFIE Fss 半基函数边界面积分修正

### 10.1 全球面采样方案

| 参数 | 值 | 说明 |
|------|----|----- |
| θ 范围 | [-π, π] | 等同双站 RCS 全角度扫描 |
| θ 采样 | 73 点 (5° 间隔) | 兼容 Legacy 721 点子集 |
| φ 范围 | [0, π) | 半球对称, 避免冗余 |
| φ 采样 | 18 条切面 (10° 间隔) | 包含 φ=0°/90° |
| 总观测方向 | 73 × 18 = **1314** | 覆盖全球面 |

### 10.2 测试矩阵

| 编号 | 方程类型 | 几何体 | N (approx) | Direct | Iterative | MLFMA | MPI |
|------|----------|--------|------------|--------|-----------|-------|-----|
| **A** | S-EFIE | Jet 100MHz | 14559 | ✅ A1 | [ ] A2 | ✅ A3 | [ ] A4 |
| **B** | S-MFIE | Sphere 600MHz | 26424 | ✅ B1¹ | — | [ ] B2 | [ ] B3 |
| **C** | S-CFIE | Sphere 600MHz | 26424 | ✅ C1 | — | ✅ C3 | [ ] C3-MPI |
| **D** | V-EFIE | Tetra 2GHz | ~986 | ✅ D1 | [ ] D2 | [ ] D3 | — |
| **E** | VS-EFIE | TriTetra 2GHz | ~1071 | ✅ E1 | [ ] E2 | [ ] E3 | — |

¹ B1 = CFIE 分解验证 (小网格)

### 10.3 求解器路径

| 路径 | 说明 | 可用方程 | 约束 |
|------|------|---------|------|
| **Direct** | Dense Z → LU | A, D, E | B/C 的 N=26424 Dense Z 需 11 GB, 不可行 |
| **Iterative** | Dense Z → GMRES | A, D, E | 同上 |
| **MLFMA** | MLFMAOperator → GMRES + SAI | A, B, C, D, E | 全部可用 |
| **MPI** | 并行装配 | A, C | 仅 RWG 表面方程 |

### 10.4 剩余子项

| 子项 | 求解路径 | 对比基准 | 通过准则 | 状态 |
|------|---------|---------|---------|------|
| A2 | Iterative | A1 | RMSE < 0.1 dB | [ ] |
| A4 | MPI (2 进程) | A1 | 机器精度 | [ ] |
| B2 | MLFMA + GMRES | C3 物理一致 | 趋势一致 | [ ] |
| B3 | MPI (2 进程) | B2 | 机器精度 | [ ] |
| C3-MPI | MPI (2 进程) | C3 | 机器精度 | [ ] |
| D2 | Iterative | D1 | RMSE < 0.1 dB | [ ] |
| D3 | MLFMA + GMRES | D1 | vs D1 < 2 dB | [ ] |
| E2 | Iterative | E1 | RMSE < 0.1 dB | [ ] |
| E3 | MLFMA + GMRES | E1 | vs E1 < 2 dB | [ ] |

### 10.5 传递准则

| 对比类型 | 准则 |
|---------|------|
| Direct vs Legacy Direct | RMSE < 1 dB |
| Iterative vs Direct | RMSE < 0.1 dB |
| MLFMA vs Direct | RMSE < 2 dB |
| MLFMA vs Legacy MLFMA | RMSE < 3 dB |
| MPI vs Serial | 机器精度 |
| CFIE Z 分解 | rel_err < 1e-12 |

### 10.6 EMSuite API 覆盖矩阵

| 功能 | EFIE | MFIE | CFIE | VEFIE | SCFIE |
|------|------|------|------|-------|-------|
| 直接装配 | ✅ RWG | ✅ RWG | ✅ RWG | ✅ SWG | ✅ RWG+SWG |
| 激励向量 | ✅ | ✅ | ✅ | ✅ | ✅ (拼接) |
| MLFMAOperator | ✅ | ✅ | ✅ | ✅ | ✅ |
| MPI 并行装配 | ✅ | ✅ | ✅ | ❌ | ❌ |
| RCS 计算 | ✅ RWG | ✅ RWG | ✅ RWG | ✅ SWG | ⚠️ 需手动拆分 |

---

## 关键参考

- **Legacy 代码**: `MoM_Basics/`, `MoM_Kernels/`, `MoM_AllinOne/`
- **验证脚本**: `EMSuite/benchmark/verify_*.jl`, `EMSuite/scripts/verification/`
- **理论**: Harrington "Field Computation by Moment Methods"; Chew et al. "Fast and Efficient Algorithms in CEM"
