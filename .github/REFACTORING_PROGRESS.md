# EMSuite 重构进度

> 最后更新: 2026-03-04

## 当前阶段: Phase 8 (性能优化)

---

## 已完成 ✅

### Phase 1: 基础架构 (2025-12)
- [x] 项目结构 `EMSuite.jl` 创建
- [x] `Project.toml` 依赖管理
- [x] Core 模块: `Interfaces.jl`, `Types.jl`, `Constants.jl`, `Materials.jl`, `Sources.jl`
- [x] Utilities: `Logging.jl`, `Parameters.jl`
- [x] CI/CD: `.github/workflows/CI.yml`
- [x] 文档框架: Documenter.jl 配置

### Phase 2: 几何与基函数 (2025-12)
- [x] 网格类型: `TriangleMesh`, `TetrahedraMesh`, `HexahedraMesh`
- [x] 网格 I/O: Nastran (`.nas`), Gmsh (`.msh`)
- [x] 坐标变换, 高斯求积
- [x] RWG 基函数 (与 Legacy 100% 匹配, 包括边排序逻辑)
- [x] SWG, RBF, PWC 基函数

### Phase 3: 积分方程 (2025-12 ~ 2026-01)
- [x] EFIE: PEC Plate/Sphere 验证, 奇异项 ($F_1$, $F_2$) 修正
- [x] MFIE: `mfie_interaction!` in-place 组装
- [x] CFIE: PEC Sphere 验证 (RMSE 2.09 dB vs Mie)
- [x] VEFIE: 体积积分方程, SWG 基函数支持
- [x] SCFIE: 面体耦合积分方程

### Phase 4: MLFMA (2026-01)
- [x] 八叉树构建 (`Octree.jl`, `OctreeBuilder.jl`)
- [x] 聚合 (`Aggregation.jl`), 含标量势项
- [x] 转移 (`Translation.jl`), Legacy 因子 $-jk/16\pi^2$ 对齐
- [x] 解聚 (`Disaggregation.jl`)
- [x] Lebedev 球面插值集成
- [x] 近场一致性: Max Diff < 1e-12
- [x] 远场精度: 修正 1/4 因子, 近邻缓冲区 = 4

### Phase 5: 求解器与并行 (2026-01)
- [x] Direct Solver (LU)
- [x] GMRES (误差 2.7e-7)
- [x] BiCGSTAB
- [x] ILU 预条件器
- [x] SPAI 预条件器
- [x] MPI 分布式并行 (n=2 vs n=1 机器精度匹配)
- [x] 多线程并行 (4 线程加速验证)

### Phase 6: 后处理与 I/O (2026-01)
- [x] RCS 计算
- [x] FarField / NearField 计算
- [x] VTK 导出 (ParaView)
- [x] 结果文件 I/O (HDF5, CSV, TXT)
- [x] 电流分布后处理

### Phase 7: 验证与对齐 (2026-02-28) ✅
- [x] SCFIE MLFMA 近场验证: Rel Err = 1.58e-15 (机器精度)
- [x] SCFIE MLFMA 远场验证: Overall 0.85%, Surface 0.85%, Volume 6.1%
- [x] Standalone EFIE MLFMA: 0.66% (4GHz, TriTetra.nas)
- [x] Standalone VEFIE MLFMA: 1.39% (4GHz, TriTetra.nas)
- [x] MoM_AllinOne 全部算例对标完成

**SCFIE MLFMA 修复的 7 个 Bug:**
1. `Disaggregation.jl`: SCFIE `efie_factor` 从 `1.0+0im` → `jkη/(16π)`
2. `MLFMAOperator.jl`: 近场 SS 块移除多余 `eta` (避免 MFIE 项双重乘 η)
3. `MLFMAOperator.jl`: VV 块从 `vefie_element_interaction` (c1=-jωμ₀κ) → `vefie_element_interaction_kernel` (c1=+jωμ₀κ)
4. `MLFMAOperator.jl`: VV 块添加缺失的 mass matrix (自交互项)
5. `MLFMAOperator.jl`: VV/SV/VS 块创建 `distribute_term_nosign!` 避免 bfsSign 双重计数
6. `Disaggregation.jl`: SWG `const_factor` 从 `-jkη` → `jkη/(4π)` (VEFIE 用 G=e^{-jkR}/(4πR))
7. `MLFMAOperator.jl`: 添加 VEFIE 缓存预计算 (TetBasisCache + precompute_vefie_basis)
8. `SCFIE.jl`: Z_SV 符号修正 `(term1 - term2)` → `(term1 + term2)` (恢复 $L$ 算子互易性)
9. `SCFIE.jl`: Z_VS 系数修正 `c1_vs = -jωμ₀` → `+jωμ₀` (统一与 Legacy 一致的正号约定)
10. 添加 `test/test_scfie.jl` 回归测试 (互易性、Direct 组装、MLFMA 近场/远场)

### Phase 7.5: ~3 dB 系统偏差根因修复 (2026-02-28) ✅

**根因 1 — `edgev̂`/`edgen̂` 方向反转** (`BasisUtilities.jl`):
- EMSuite 计算 `e1 = v2 - v3`（v3→v2 方向），Legacy 计算 `v3 - v2`（v2→v3 方向）
- 导致 `edgen̂ = cross(edgev̂, facen̂)` 指向三角形内部而非外部
- 影响: `faceSingularityIgIvecg` 中 `p02jl`（投影距离）符号错误 → 近奇异积分结果错误
- 修复: 将边向量改为 `e1 = v3 - v2`, `e2 = v1 - v3`, `e3 = v2 - v1`

**根因 2 — `calc_near_interaction!` 经验因子 + 面积归一化** (`EFIE.jl`):
- 存在经验因子 `* 1.25`（应为 1.0）
- `inv_areas = 1.0 / tri_test.area`（应为 `1.0 / (tri_test.area * tri_source.area)`）
- 在 edgev 方向错误时，这两个错误部分互相补偿（Z_near ≈ 0.0485× 正确值 → 近乎可忽略）
- 修复: 移除 `* 1.25`，改为正确的双面积归一化

**验证结果:**
- Z 矩阵对角线 ratio: 1.000000 (2640×2640 板网格 vs Legacy)
- Frobenius 范数比: 1.000010
- RCS (Jet 100MHz): Mean Diff 0.05 dB / RMSE 0.29 dB (Phi=0), Mean Diff 0.008 dB / RMSE 0.09 dB (Phi=90)
- 全部 138/138 单元测试通过（无回归）

### Phase 10.A: MLFMA 因子修复 (2026-03-02) ✅

**Bug 1 — MLFMA far-field ×4 因子** (`MLFMAOperator.jl`, `Disaggregation.jl`):
- 根因: `efie.factor = jkη/(16π)` 包含 `1/4` (来自 RWG `l²/4` 归一化), 但 MLFMA 聚合/解聚各用 `l/2`，乘积 `l²/4` 已自然包含该因子 → **双重计数 1/4**
- Legacy 避免此问题: translation 用 `-jk/(16π²)` + disagg 用 `jkη`（不同的因子分解方式）
- 修复: `y_far *= 4 * operator.factor` (EFIE mul!), CFIE/SCFIE disagg `efie_factor *= 4`
- 验证: 自洽性系数误差 65.7% → 0.30%, RCS RMSE 3.1 dB → 0.028 dB
- A3 S-EFIE MLFMA vs Legacy: Mean Diff 0.048 dB, RMSE 0.303 dB

**Bug 2 — CFIE MLFMA MFIE 符号错误** (`Disaggregation.jl`):
- 根因: MFIE K 算子使用 $\nabla_{r'}G$（源梯度），远场近似为 $+jk\hat{k}G$;
  代码错误地使用了 $\nabla_r G$（场梯度）的 $-jk\hat{k}G$，导致 MFIE 项符号反转
- Legacy 处理: 聚合/解聚分离, 统一乘以 `jkη`, EFIE 和 MFIE 远场系数同号 (+jkη)
- 修复: `(-efie_factor)` → `(+efie_factor)` (MFIE 项)
- 验证: C3 CFIE MLFMA vs Legacy: RMSE 3.45 dB → **0.003 dB** (1000× 改善)
- GMRES 迭代次数: 50 → 7（算子准确度提升后收敛加速）

**已验证测试结果汇总:**

| 测试 | 指标 | 结果 |
|------|------|------|
| 单元测试 | 138/138 | ✅ PASS |
| A1 S-EFIE Direct Jet | RMSE vs Legacy | 0.215 dB |
| A2 S-EFIE Iterative | — | ⏭ 跳过 (N=14559, 全 GMRES 需 3.4GB) |
| A3 S-EFIE MLFMA Jet | RMSE vs Legacy | 0.303 dB |
| A3 self-consistency | 系数误差 | 0.30% |
| B1 CFIE 分解 | rel_err | 0.0 (10/10) |
| B2 S-MFIE MLFMA Sphere | RCS 趋势 vs C3 | 物理一致 (GMRES 32 iter) |
| C1 S-CFIE Direct Sphere | RMSE vs Legacy | 0.001 dB |
| C3 S-CFIE MLFMA Sphere | RMSE vs Legacy | **0.003 dB** |
| D1-SWG V-EFIE Direct | RMSE vs Legacy | 0.952 dB |
| D2 V-EFIE Iterative | RMSE vs D1 | **0.000089 dB** |
| D3 V-EFIE MLFMA | RMSE vs D1 | **0.0000 dB** |
| E1 VSEFIE Direct | RMSE vs Legacy | **0.602 dB** |
| E2 VS-EFIE Iterative | RMSE vs E1 | **0.000327 dB** |
| E3 VS-EFIE MLFMA | RMSE vs E1 | **0.0000 dB** |
| EFIE MLFMA Sphere | RMSE vs Legacy SCFIE | 0.041 dB |

### Phase 10.C: MFIE MLFMA 支持 + 迭代求解器验证 (2026-03-04) ✅

**Bug 1 — MLFMA 近场缺少 MFIE 分派** (`MLFMAOperator.jl`):
- 根因: 近场 surface-surface 组装只有 `efie_interaction!` 分支, MFIE 算子会错误调用 EFIE
- 修复: 添加 `elseif operator isa MFIE` → `mfie_interaction!()` 分支

**Bug 2 — MLFMA 远场缺少 MFIE 测试路径** (`Disaggregation.jl`):
- 根因: 解聚只有 EFIE 的 ρ·E 测试, 缺少 MFIE 的 (ρ×n̂)·H 磁场测试
- 修复: 添加 `is_mfie` 标志, MFIE 因子 `4·jkη/(16π)`, 磁场测试 `H = (E_θ·φ̂ - E_φ·θ̂)·phase/η`

**迭代求解器验证:**
- A2 跳过: N=14559, 全 GMRES 需 3.4GB Krylov 基, 不可行; GMRES(50)+Diagonal 预条件器收敛差
- D2 VEFIE (全 GMRES, restart=N, tol=1e-6): RMSE = 0.000089 dB ✅ (~680 次迭代收敛)
- E2 SCFIE (全 GMRES, restart=N, tol=1e-6): RMSE = 0.000327 dB ✅ (~676 次迭代收敛)

**MLFMA 体积方程自洽性:**
- D3 VEFIE MLFMA: RMSE = 0.0000 dB ✅ (GMRES 1 次迭代收敛)
- E3 SCFIE MLFMA: RMSE = 0.0000 dB ✅ (GMRES 1 次迭代收敛)

**B2 MFIE MLFMA (Sphere 600MHz, N=26424):**
- MLFMA setup 375.3s, GMRES 32 次迭代 / 130.6s
- 残差: 6.04e-2 → 6.36e-4
- RCS 物理趋势与 C3 CFIE 一致

### Phase 10.B: SCFIE Fss 边界修正 (2026-03-03) ✅

**Bug — 缺失半基函数边界面积分修正 (Fss)** (`SCFIE.jl`, `MLFMAOperator.jl`):
- 根因: 边界 SWG 基函数（半基函数，仅有一个四面体支撑）的标量势缺失表面积分修正项
- Legacy 在 `EFIEVSIERWGSWG.jl` 中通过 `Fss` 项处理：
  - `isbdn=true` 时: Z_VS[n,m] += jωμ₀/(4πk²) × l_m × |A_n| × ∫∫ G dS_tri dS_face
  - `δκ≠0` 时: Z_SV[m,n] += δκ × (同上)
- EMSuite 完全缺失此修正 → 耦合矩阵 Z_SV/Z_VS 偏差 22%, Z_VV 偏差 48%
- 修复: 在 `SCFIE.jl` 中添加 `assemble_fss_boundary_correction!` 和 `assemble_fss_boundary_correction_sparse`
  - 直接求解路径: 在 `assemble_coupling_blocks!` 后调用
  - MLFMA 路径: 以稀疏矩阵形式加到 Z_near
- 验证: E1-VSEFIE RMSE 从 **5.3 dB → 0.60 dB** (PASS), 138/138 测试全通过

### 验证里程碑
- [x] **SEFIE Direct**: `verify_SEFIE_direct.jl` (18.8s 4线程 / 30.3s 1线程, Legacy 31.2s)
- [x] **SEFIE MLFMA**: `verify_SEFIE_mlfma.jl` (Ratio 1.0000, Rel Err 1.5%)
- [x] **VEFIE Direct**: `verify_VEFIE_direct.jl` (Legacy Parity)
- [x] **VEFIE MLFMA**: `verify_VEFIE_mlfma.jl` (Ratio 1.0000, Rel Err 0.04%)
- [x] **SCFIE Direct**: `verify_SCFIE_direct.jl` (VSIE Plate+Metal, RCS -15.35 dBsm)
- [x] **SCFIE MLFMA**: `quick_scfie_mlfma_test.jl` (Near-field 1.58e-15, Far-field 0.85%)
- [x] **Standalone EFIE MLFMA**: `test_efie_vefie_farfield.jl` (Rel Err 0.66%)
- [x] **Standalone VEFIE MLFMA**: `test_efie_vefie_farfield.jl` (Rel Err 1.39%)
- [x] **MPI**: `benchmark_parallel_sphere.jl` (Consistency check passed)
- [x] **Threading**: 4 threads vs 1 thread speedup 验证

---

## 进行中 🔧

### Phase 8: 性能优化 (当前)

> **目标**: 相同用例全流程耗时 ≤ Legacy (保底一致)，争取 ≤ 0.5× Legacy (2× 加速)
> 详细计划见 `REFACTORING_ROADMAP.md` Phase 8

---

## 已完成 (最近) ✅

### Phase 10: 全方程全路径精度对齐 (2026-03-04) ✅

> 详细计划见 `REFACTORING_ROADMAP.md` Phase 10

#### 10.0 全球面采样方案
- **标准网格**: θ ∈ [-π, π] 73 点 × φ ∈ [0, π) 18 条切面 = 1314 观测方向
- **输出格式**: CSV (theta_deg, phi_deg, RCS_theta_dB, RCS_phi_dB, RCS_total_dB)
- **状态**: 设计完成，待实施

#### 10.1 测试矩阵 (5 方程 × 4 路径)

| 编号 | 方程 | 几何 | Direct | Iterative | MLFMA | MPI |
|------|------|------|--------|-----------|-------|-----|
| A | S-EFIE | Jet 100MHz | ✅ A1 | ⏭ A2³ | ✅ A3 | [ ] A4 |
| B | S-MFIE | Sphere 600MHz | ✅ B1³ | — | ✅ B2 | [ ] B3 |
| C | S-CFIE | Sphere 600MHz | ✅ C1 | — | ✅ C3 | [ ] C3-MPI |
| D | V-EFIE | Tetra 2GHz | ✅ D1 | ✅ D2 | ✅ D3 | — |
| E | VS-EFIE | TriTetra 2GHz | ✅ E1 | ✅ E2 | ✅ E3 | — |

³ B1 = CFIE 分解验证 (小网格), 非 Direct 求解; A2 跳过 (N=14559 全 GMRES 需 3.4GB Krylov 基, 不可行)

#### 10.2 实施进度

| 步骤 | 内容 | 状态 |
|------|------|------|
| Step 1 | 全球面 Legacy 基线生成 (8 用例) | ✅ 完成 |
| Step 2 | S-MFIE CFIE 分解验证 | ✅ B1 PASS |
| Step 3 | EMSuite Direct 基准 (A1, D1, E1) | ✅ 全部 PASS |
| Step 4 | EMSuite Iterative 基准 (D2, E2) | ✅ 全部 PASS |
| Step 5 | EMSuite MLFMA 基准 (A3, B2, C3, D3, E3) | ✅ 全部 PASS |
| Step 6 | MPI 基准 (A4, B3, C3-MPI) | ⏭ 延后 |
| Step 7 | 全球面误差热力图 + 报告 v3 | ✅ 完成 |

### Phase 9: 代码质量与发布 (暂缓)

#### 9.1 测试套件清理 ✅
- **状态**: 完成
- **结果**: 138/138 测试全部通过 (之前: 119 pass, 4 fail, 3 error)
- **修复内容**:
  1. `test_basis_functions.jl` (SWG): 修正基函数数量期望值 (1→7), 修正内部 BF 查找逻辑
  2. `test_postprocessing.jl`: 修正 API 调用 (旧: `trianglesInfo, RWG` → 新: `basis`)
  3. `test_mlfma.jl`: 修正 `TriangleMesh` 构造 (`trinum=1` → `size(elements,2)`)
  4. `test_integration.jl`: 修正 `radarCrossSection`/`geoElectricJCal` API 调用
  5. `Aggregation.jl`/`Disaggregation.jl`: 修正 `Vector{AbstractBasisFunction}` → `Vector{<:AbstractBasisFunction}` 类型派发
  6. `Disaggregation.jl`: 添加 `disaggregate_leaf!` 单基函数便捷包装器

#### 9.2 JuliaFormatter
- **状态**: 待开发

#### 9.3 API 文档
- **状态**: 待开发

---

## 待开始 📋

### Phase 8: 性能优化 (当前)

> **目标**: 相同用例全流程耗时 ≤ Legacy (保底一致)，争取 ≤ 0.5× Legacy (2× 加速)

#### 8.0 性能基线测量 ✅
- [x] `benchmark/performance_baseline.jl` 已创建
- [x] EMSuite 7 用例 + Legacy 5 用例分阶段计时
- [x] `test_results/PERFORMANCE_BASELINE.md` 已生成

#### 8.1 Z 组装去锁 ✅
- [x] `Impedance.jl` 全局 SpinLock → Per-row SpinLock 数组 (N 把锁)
- [x] 争用概率: ~nthreads/N ≈ 0.03%，缓存友好的 3×3 立即写入
- [x] 138/138 测试通过
- [x] 实测: Plate EFIE **-54%**, Jet EFIE **-12%**, Jet CFIE -2%

#### 8.2 CFIE/MFIE 内核优化 ✅
- [x] MFIE 预计算高斯点 (消除 ~94M 堆分配)
- [x] MFIE 循环重排 (i,j)外(m,n)内, 消除 9× 冗余 rvec/R 计算
- [x] MFIE 积分阶数 7→4 点 (对齐 Legacy GQPNTri=4)
- [x] 138/138 测试通过
- [x] CFIE 168.29s → **43.48s** (-74%), CFIE/EFIE = **2.31×** ≤ 2.5× ✅

#### 8.3 MLFMA Z_near 去锁 + CSC 预分配
- [ ] 近场稀疏矩阵组装并行化
- [ ] CSC 格式 nnz 预估 + 预分配

#### 8.4-8.7 后续优化
- [ ] 内存分配热点消除
- [ ] 类型稳定性审查
- [ ] SIMD / LoopVectorization
- [ ] 预条件器优化

#### 8.8 最终复测
- [ ] 全部用例重新计时，对比 Phase 8.0 基线
- [ ] 生成 `test_results/PERFORMANCE_REPORT.md`

**已知数据点** (Phase 10 期间实测, 4 线程, Windows 11):

| 用例 | N | EMSuite 组装 (s) | EMSuite 求解 (s) | EMSuite 总计 (s) |
|------|---|----------------|----------------|----------------|
| Jet EFIE Direct | 14559 | 20.3 | 34.5 | 54.8 |
| Jet CFIE Direct | 14559 | 180.5 | 16.0 | 196.5 |
| Jet EFIE MLFMA | 14559 | 56.2 (setup) | 7.3 | 63.4 |
| Sphere CFIE MLFMA | 26424 | 131.0 (setup) | 7.6 | 138.6 |
| Sphere MFIE MLFMA | 26424 | 375.3 (setup) | 130.6 | 505.9 |

### Phase 9: 代码质量与发布 (剩余)
- [x] 测试套件清理: 138/138 全部通过
- [ ] JuliaFormatter.jl 统一代码风格
- [ ] 测试覆盖率统计与提升 (目标 > 80%)
- [ ] API 文档补全 (所有公共接口)
- [ ] 用户教程 (Quick Start, Advanced)
- [ ] 理论文档 (MoM, MLFMA, 积分方程推导)
- [ ] CHANGELOG.md 完善
- [ ] 发布到 Julia General Registry

---

## Legacy 因子对照表

| 项目 | Legacy | EMSuite | 说明 |
|------|--------|---------|------|
| EFIE Z | $1/16\pi$ | $1/4\pi$ + 显式 $l/2A$ | 数学等价 |
| FarField | $1/4\pi$ | $1/4\pi$ | 一致 |
| Translation | $1/16\\pi^2$ | $-jk/16\\pi^2$ | Legacy 对齐 |\n| SWG Disagg | N/A | $jk\\eta/(4\\pi)$ | VEFIE G 含 $1/(4\\pi)$ |
| 时间约定 | $e^{-j\omega t}$ | $e^{j\omega t}$ | Z 虚部符号相反 |
| Area/Length | 隐式包含在因子中 | 显式归一化 | 结果等价 |

---

## 已知问题

1. ~~**EMSuite vs Legacy ~3 dB 系統偏差**~~ — **已修复** (2026-02-28): 根因为 `edgev̂` 方向反转 + `calc_near_interaction!` 经验因子。修复后 RCS 偏差 < 0.3 dB RMSE。
2. ~~**MLFMA 远场 ×4 因子 + CFIE 符号**~~ — **已修复** (2026-03-02): EFIE MLFMA 系数误差 65.7% → 0.30%. CFIE MLFMA RMSE 3.45 dB → 0.003 dB.
3. **VEFIE Mie 偏差**: Legacy 和 EMSuite 均比 Mie 级数低 ~25dB — 属于 Legacy 算法固有问题, 标记为 "Legacy Parity", 物理修正为未来研究课题
4. **BiCGSTAB 收敛**: 需要预条件才能可靠收敛
5. **SCFIE 耦合项互易性**: 已修复 — Z_SV/κ = Z_VS^T 在机器精度成立 (2.99e-16)
6. ~~**EFIE 闭合体内部谐振**: EFIE 用于闭合导体时条件数差 (Direct vs MLFMA 系数差 62%)，应改用 CFIE~~ — **已确认**: 现在 CFIE MLFMA 正确工作 (RMSE 0.003 dB, 7 iterations)
7. **SWG MLFMA const_factor 符号**: `const_factor = jkη/(4π)` 可能应为 `-jkη/(4π)` (VEFIE `c1 = -jkηκ` 含负号). 需要 VEFIE MLFMA 精度测试验证.

---

## 更新日志

| 日期 | 更新内容 |
|------|----------|
| 2026-03-04 | **Phase 8.2 CFIE/MFIE 内核优化完成** — MFIE 零分配+循环重排+4点积分。CFIE 168.29→43.48s (-74%), CFIE/EFIE=2.31× ≤ 2.5× 达标, vs Legacy 1.19× |
| 2026-03-04 | **Phase 8.1 Z 组装去锁完成** — `Impedance.jl` 全局 SpinLock → Per-row SpinLock 数组 (N 把锁)。实测 Plate EFIE -54% (1.02→0.47s), Jet EFIE -12% (20.7→18.3s), EFIE 已与 Legacy 持平。138/138 测试通过 |
| 2026-03-04 | **Phase 8.0 基线完成** — EMSuite 7 用例 + Legacy 5 用例计时。见 `test_results/PERFORMANCE_BASELINE.md`。Legacy MLFMA 用例因 Julia 1.12 threadid() 兼容性失败 |
| 2026-03-04 | **Phase 10 完成** — MFIE MLFMA 支持 (近场+远场), 迭代求解器验证 (D2/E2 PASS), MLFMA 体积方程 (D3/E3 PASS), B2 MFIE MLFMA PASS. 全部 12/12 子测试 PASS (A2 跳过, MPI 延后) |
| 2026-03-03 | **Phase 8 性能优化计划** — 加入性能优化路线: 6 热点 (SpinLock去锁/CFIE合并/MLFMA Z_near/内存/SIMD/类型稳定), 8 步骤, 目标 ≤ Legacy 保底, ≤ 0.5× Legacy 挑战 |
| 2026-03-03 | **SCFIE Fss 边界修正** — 半基函数边界面积分修正。E1-VSEFIE RMSE 5.3→0.60 dB. D1-SWG VEFIE RMSE 0.95 dB. 138/138 测试通过 |
| 2026-03-02 | **MLFMA 因子修复×2** — (1) EFIE far-field ×4 因子: 系数误差 65.7%→0.30%, RMSE 3.1→0.028 dB; (2) CFIE MFIE 符号: ∇_{r'}G 给出 +jk k̂ (非 -jk k̂), RMSE 3.45→0.003 dB, GMRES 50→7 迭代 |
| 2026-03-01 | **Phase 10 计划** — 全方程全路径精度对齐设计完成: 5 方程 (S-EFIE/S-MFIE/S-CFIE/V-EFIE/VS-EFIE) × 4 路径 (Direct/Iterative/MLFMA/MPI), 全球面 1314 点采样, 共 17 子测试项 |
| 2026-02-28 | **精度效率报告 v2** — 全面基准测试: SEFIE Direct(RMSE 0.29dB), CFIE Direct, SEFIE MLFMA, SCFIE MLFMA(Sphere N=26424). 见 `test_results/ACCURACY_EFFICIENCY_REPORT.md` |
| 2026-02-27 | 修正验证脚本 benchmark/run_full_benchmark.jl API 错误 (MLFMAOperator 构造/排序透明性) |
| 2026-02-27 | 新增 MLFMA MatVec 快速测试脚本 (benchmark/quick_matvec_test.jl) |
| 2026-02-27 | 新增 Direct vs MLFMA 自一致性测试 (benchmark/self_consistency_test.jl) |
| 2026-02-28 | **测试套件清理完成** — 138/138 全部通过, 修复 6 个预存测试问题 |
| 2026-02-28 | 修复 `Vector{AbstractBasisFunction}` 类型派发 bug (Aggregation/Disaggregation) |
| 2026-02-28 | **SCFIE MLFMA 验证完成** — 修复 7+3 个 bug, 近场 1.58e-15, 远场 1.04% |
| 2026-02-28 | **SCFIE 耦合互易性修复** — Z_SV `(term1-term2)` → `(term1+term2)`, Z_VS c1 符号修正 |
| 2026-02-28 | 添加 `test/test_scfie.jl` 回归测试 (9 tests all pass) |
| 2026-02-28 | **~3 dB 偏差根因修复** — `BasisUtilities.jl` 边向量方向 + `EFIE.jl` 近交互面积归一化。Z 矩阵 ratio=1.0, RCS RMSE<0.3 dB |
| 2026-02-28 | **精度效率报告 v2** — 全面基准测试: SEFIE Direct(RMSE 0.29dB), CFIE Direct, SEFIE MLFMA, SCFIE MLFMA(Sphere N=26424). 见 `test_results/ACCURACY_EFFICIENCY_REPORT.md` |
| 2026-02-27 | 初始化进度文件, 从 copilot-instructions.md 迁移 |
| 2026-01-xx | VEFIE MLFMA 验证完成 (Rel Err 0.04%) |
| 2026-01-xx | SCFIE Direct 验证完成 (VSIE Plate) |
| 2026-01-xx | MPI/Threading 并行验证通过 |
| 2025-12-xx | Surface IE (EFIE/MFIE/CFIE) 全面验证完成 |
| 2025-12-xx | Legacy 因子对齐完成, 移除经验常数 |
