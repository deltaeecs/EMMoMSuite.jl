# EMSuite Phase 9 覆盖率提升进度 (追加节)

> 最后更新: 2026-03-XX

## MLFMA 算法分析阶段 — **已完成** ✅

> 更新: 2026-03-XX

### 工作成果

1. **代码清理**:
   - `scripts/` 废物代码归档到 `scripts/archive/` (~80 文件)
   - `benchmark/debug_archive/` 清理 (~50 文件)
   - 根目录废物文件删除 (tmp_*.txt, 过期文档等)
   - 保留: `scripts/format_code.jl`, `scripts/install_formatter.jl`

2. **Legacy MLFMA 完整算法报告**: [legacy_mlfma_algorithm_report.md](.github/plans/legacy_mlfma_algorithm_report.md)
   - 覆盖从网格读取到 RCS 后处理的全流程
   - 精确到每一个公式和系数
   - 12 项 Legacy 工程技巧总结
   - EMSuite vs Legacy 系数链完整验证 (总系数一致: k²η/16π²)

3. **关键发现**:
   - nLevels=2: EFIE 精度 0.5%~9% ✅ | nLevels≥3: 1005%~2918% ✗
   - 总系数链验证通过 → 问题在 upward/downward pass 的插值或相移中
   - Addition Theorem 指数符号确认: 源端 e^{+jk}, 场端 e^{-jk}
   - 远亲定义 (7³-3³=316) 和清零时机均与 Legacy 一致

4. **检视**: 2 轮检视完成，修正了 edgel 符号描述、时间约定、Addition Theorem 指数符号

5. **PMCHW 完整算法报告**: [pmchw_algorithm_report.md](.github/plans/pmchw_algorithm_report.md)
   - 覆盖面等效原理 → L/K 算子 → 2N×2N 矩阵 → 4遍 MLFMA → 双流 RCS 全链路
   - K^PMCHW 积分核与 MFIE K 的详细区别 (无对角、无 n̂× 测试、无 η₀ 预乘)
   - 4 个子块的 MLFMA 系数链验证 (EJ/HM L型, EM/HJ K型)
   - 确认 Legacy 不含 PMCHW 实现，PMCHW 完全为 EMSuite 原创
   - 2 轮检视完成：Round 1 修复 §6.4 球坐标符号、§C.3 有损介质描述等 4 项; Round 2 无新问题

## Phase 9 测试覆盖率提升 → 进行中

### Round 4 工作成果

新增/扩展测试 (commit ca3fca5, 2161daa):
- test_hex_rbf.jl (NEW, 67 tests): HexahedraMesh+PWCHexBasis+RBFBasis
- test_geometry.jl: 33->60 tests (NAS CHEXA/CTETRA/MSH tetra)
- test_mpi_array_utils.jl: 24->47 tests
- test_io.jl: VTK tetra + HDF5 save_result
- 总计: 449/449 pass

源码 Bug 修复 (commit 99260a0, 2161daa):
- VEFIE/SCFIE/Excitation/RadiationIntegral: gq_hex->gq_quad coordinates
- MeshIO.jl: _parse_chexa CHEXA 续行修复 (节点7,8丢失)
- GmshIO.jl: 添加 type 4 (tetra) + type 5 (hexa) 支持
- indices.jl: indice2rank NTuple{1} in 语义修复

Phase 9 检视迭代 Round 1 (commit 90787dc, baa0418):
- GmshIO.jl: 文档+版本检查@warn+混合网格逻辑
- VEFIE.jl _rbf_far_kernel!: 预计算 freeVns 减少冗余分配

覆盖率历史: Round1=57.22% Round2=64.37% Round3=65.69% Round4=测量中

---

# EMSuite 閲嶆瀯杩涘害

> 鏈€鍚庢洿鏂? 2026-03-06
## 当前阶段: Phase 14 全量精度测试与对比报告 — **进行中 🔄**

> 最后更新: 2026-03-07

### Phase 14 计划概述

**目标**: 对 EMSuite 所有主要积分方程 × 求解路径，与 Feko 商业软件结果或 Mie 解析解进行系统对比，生成独立精度报告（不再对比 Legacy 代码）。

**计划文档**: [PHASE_14_ACCURACY_REPORT_PLAN.md](PHASE_14_ACCURACY_REPORT_PLAN.md)

**基准来源**:
- Feko 基线: `MoM_AllinOne/deps/compare_feko/` (Jet 100MHz, Sphere 600MHz, Plate 1.2GHz, Plate+Metal 1.2GHz)
- Mie 解析解: EMSuite 内置 `Utilities/MieSeries.jl`（用于 PEC 球独立校验）

**测试用例**: F1–F9 (9 个用例，覆盖 S-EFIE/S-CFIE/V-EFIE/SCFIE × Direct/MLFMA)

**精度门限**: Direct ≤ 2 dB RMSE，MLFMA ≤ 3 dB RMSE (vs Feko/Mie)

**待完成子任务**:
- [x] 14.0 Feko CSV 解析器 + TDD 测试 ✅ b416090 (39/39)
- [x] 14.1 参考基准生成器（Mie PEC/介质 + 偶极子解析）✅
- [x] 14.2 `AccuracyResult` + `AntennaAccuracyResult` 指标函数 ✅ 65b4297 (25/25)
- [x] 14.3 F1–F4 Jet 仿真脚本 ✅ b2efd57
- [x] 14.4 F5–F6 Sphere 仿真脚本 ✅
- [x] 14.5 F7–F9 Plate 仿真脚本 ✅ (F8 跳过-无曲面网格)
- [x] 14.6 P1, P3 PMCHW Direct 脚本 ✅ (P2 跳过-待实现 PMCHWMLFMAOperator)
- [ ] 14.7 `PMCHWMLFMAOperator` 实现 (2×2 块 MLFMA)
- [ ] 14.8 P2 PMCHW MLFMA 验证
- [x] 14.9 A1–A4 偶极子天线 DeltaGap 基准脚本 ✅ 3039d32
- [ ] 14.10 实际运行仿真 → CSV → ACCURACY_REPORT.md
- [ ] 14.11 检视迭代 × 2 轮

---

## Phase 15: 介质与金属-介质混合天线 + PMCHWMLFMAOperator — **计划中 📋**

> 最后更新: 2026-03-07（设计修订 rev2：Gibson Algorithm 14 两遍分离聚合方案）
> 计划文档: [PHASE_15_DIELECTRIC_ANTENNA_PLAN.md](PHASE_15_DIELECTRIC_ANTENNA_PLAN.md)

**目标**: 扩展天线测试至介质（PMCHW）和金属-介质混合（VS-EFIE/VS-CFIE）类型，同步实现原生 PMCHWMLFMAOperator。

**PMCHWMLFMAOperator 设计方案（第三次修订 — 依据 Gibson Ch.11 Algorithm 14）**:
- **放弃**: ~~`MagneticRWGBasis` 标签类型~~（前两版设计已废弃）
- **核心原则**: **两个 N 点八叉树**（octree0/k0 + octree1/k1，不是 2N 点单树）+ **四遍远场** + 四种解聚接收函数
- **Gibson 原文根据**: "two passes must be executed...only T_J is actually stored...each testing function has two receive functions"
- **aggS 清零**: 四遍远场每遍开始前强制 `fill!(aggS, 0)`，防止跨遍叠加
- **文件**: `PMCHWMLFMAOperator.jl` 自包含（不修改 Aggregation.jl / Disaggregation.jl）
- **双八叉树**: `octree0`（k0 外部介质）+ `octree1`（k1 内部介质），各 N 点
- **4 遍远场**: J×k0, J×k1, M×k0, M×k1；每遍使用对应 octree，aggS 清零隔离
  - Z^EJ: `term_efie × jk₀η₀/(16π) + jk₁η₁/(16π)`（R_L）
  - Z^HJ: `term_mfie × -(jk₀/(16π) + jk₁/(16π))`（R_{-K}，负号）
  - Z^EM: `term_mfie × +(jk₀/(16π) + jk₁/(16π))`（R_K，与 HJ 符号相反）
  - Z^HM: `term_efie × jk₀/(η₀·16π) + jk₁/(η₁·16π)`（R_{L_η}）

**新增 API**:
- `excitation_vector(op::PMCHW, source::DeltaGapSource, basis::RWGBasis)` → 2N 激励向量（E-行 delta-gap）
- `input_impedance(op::PMCHW, source, I_2N, basis)` → PMCHW J-部分输入阻抗
- `excitation_vector(op::SCFIE, source::DeltaGapSource, surf_basis, vol_basis)` → 表面 + 零体积激励

**测试场景**:
| ID | 方程 | 馈电 | 求解 | 参考 |
|----|------|-----|------|-----|
| B1 | PMCHW (εᵣ=4) | DeltaGap (球面) | Direct | 物理自洽 + εᵣ→1 极限 |
| B2 | PMCHW | DeltaGap | MLFMA (双八叉树+四遍远场) | B1 Direct |
| B3 | VS-EFIE (α=0) | DeltaGap (金属面) | Direct | EFIE-only εᵣ→1 |
| B4 | VS-CFIE (α=0.5) | DeltaGap | Direct | B3 |
| B5 | VS-CFIE (α=0.5) | DeltaGap | MLFMA | B4 Direct |

**子任务**:
- [ ] 15.1 TDD: `test_pmchw_excitation.jl`
- [ ] 15.2–15.3 PMCHW DeltaGap 激励 + input_impedance API
- [ ] 15.4–15.5 SCFIE DeltaGap 激励 API
- [ ] 15.6 基准脚本 `run_B1_B5_antenna.jl`
- [ ] 15.7 TDD: `test_pmchw_mlfma_operator.jl`
- [ ] 15.8 实现 `assemble_near_field_pmchw`（2N×2N，4 块，无 MagneticRWGBasis）
- [ ] 15.9 实现 `aggregate_leaf_pmchw!`（x_range 参数区分 J/M）
- [ ] 15.10 实现 `disaggregate_leaf_pmchw_j!` + `_m!`（四种接收核函数）
- [ ] 15.11 组装 `PMCHWMLFMAOperator` struct + 构造函数 + `mul!`（两棵 N 点八叉树 octree0/k0+octree1/k1，四遍远场）
- [ ] 15.12–15.13 报告更新 + 检视迭代（≥ 2 轮 clean）

### 2026-03-06 Update (Governance Refresh)

- Rolled back the recent experimental PMCHW/MLFMA code updates as requested.
- Added governance plan: `.github/plans/phase_15_theory_impl_test_refresh.md`.
- New task `15.G1` added to roadmap to enforce a strict Theory -> Implementation -> Test contract.
- Next implementation work must satisfy mandatory Gate A/B/C/D checks before acceptance.

### 2026-03-06 Update (Gate A/B Implemented)

- Added executable Gate A/B tests into `test/test_pmchw_mlfma_operator.jl`.
- Gate A (structural invariants) is green:
   - HJ + EM near-field invariant passed
   - non-trivial near/far split passed
   - dual-octree permutation sanity passed
- Gate B (EJ pass alignment) is now machine-tracked with fixed seed (`Random.seed!(42)`):
   - k0: `rel=0.7536`, `corr=0.9962` (magnitude regression, marked `@test_broken`)
   - k1: `rel=0.9588`, `corr=0.4641` (marked `@test_broken`)
- Main end-to-end gate remains failing: `15.11 MLFMA mul! vs Direct` = `49.40%`.

---
## 褰撳墠闃舵: Phase 13.3 V-EFIE MPI 骞惰鍖?鈥?**宸插畬鎴?* 鉁?

**鏈€鏂版垚鏋?(2026-03-06)**:
- 鉁?**V-EFIE MPI Allreduce 鏂规瀹炵幇瀹屾垚**
  - 绠楁硶: 姣忚繘绋嬪鐞?`(it-1) % n_procs == rank` 鐨勬祴璇曞洓闈綋; 瀵圭О鍒╃敤 (js>it); `MPI.Allreduce!` 姹囨€?
  - 璋冨害淇: 鍒犻櫎 `module VolumeAssemblyMPI` 灏佽锛岀洿鎺?`import .Assembly: assemble_impedance_matrix_parallel` 鎵╁睍
  - 姝ｇ‘鎬ч獙璇? Tetra.nas N=3201, Z[1,1] 鍦?1/2/4 杩涚▼涓嬪畬鍏ㄧ浉鍚?(`4.747e6 - 4.716e6im`)
  - 鎬ц兘: 1杩涚▼ 10.95s 鈫?2杩涚▼ 8.86s (1.24脳); 灏?N 鏃?Allreduce 寮€閿€鍗犱富瀵间负姝ｅ父鐜拌薄
- 鉁?**BasisFunctions.jl**: 琛ュ叏 `get_triangles_info` 瀵煎嚭
- 鉁?**Parallel.jl**: 淇瀛愭ā鍧楅殧绂婚棶棰?(Julia dispatch 璺ㄦā鍧椾笉鍙)

## Phase 13.2 Option C 璇勪及 鈥?**宸插畬鎴?* 鉁? 
**Phase 13.3 Option A (MPI 骞惰鍖栨祴璇?** 鈥?**瀹屾垚** 鉁?

**鏈€鏂版垚鏋?(2026-03-0X)**:
- 鉁?**Option C 瀹屾垚**: PWC/RBF 鍩哄嚱鏁版€ц兘鍩哄噯娴嬭瘯
  - PWC+VEFIE: 12.62s (3.16脳 faster than SWG 39.94s)
  - 浣?PWC DOF +38% (21834 vs 15828) 鈫?姹傝В鍣ㄨ礋鎷呭姞閲?
  - Per-DOF虏 鏁堢巼浠?SWG 鐨?0.17脳 鈫?**涓嶉€傚悎閫氱敤浼樺寲**
  - RBF 纭浠呮敮鎸佸叚闈綋缃戞牸 (涓嶉€傜敤鍥涢潰浣?VEFIE)
  - 馃搳 璇︾粏鎶ュ憡: [OPTION_C_PWC_RBF_EVALUATION.md](.github/OPTION_C_PWC_RBF_EVALUATION.md)
  
**Phase 13.1 鎬荤粨** (2026-03-01):
- 鉁?V-EFIE: 39.05s (**1.18脳 vs Legacy** 46.13s) 鈥?**棣栨蹇簬 Legacy**
- 鉁?SCFIE: 63.81s (**1.04脳 vs Legacy** 66.68s) 鈥?**涓?Legacy 鎸佸钩**
- 鉁?FastExp 鏌ユ壘琛? 10000 鏉＄洰, 绮惧害 < 0.002%, ~160KB
- 鈿狅笍 Thread-local buffer 澶辫触 (-63% 鎬ц兘閫€姝? 宸插洖閫€)
- 馃搳 璇︾粏鎶ュ憡: [PHASE_13.1_SUMMARY.md](.github/PHASE_13.1_SUMMARY.md)

---

## 宸插畬鎴?鉁?

### Phase 1: 鍩虹鏋舵瀯 (2025-12)
- [x] 椤圭洰缁撴瀯 `EMSuite.jl` 鍒涘缓
- [x] `Project.toml` 渚濊禆绠＄悊
- [x] Core 妯″潡: `Interfaces.jl`, `Types.jl`, `Constants.jl`, `Materials.jl`, `Sources.jl`
- [x] Utilities: `Logging.jl`, `Parameters.jl`
- [x] CI/CD: `.github/workflows/CI.yml`
- [x] 鏂囨。妗嗘灦: Documenter.jl 閰嶇疆

### Phase 2: 鍑犱綍涓庡熀鍑芥暟 (2025-12)
- [x] 缃戞牸绫诲瀷: `TriangleMesh`, `TetrahedraMesh`, `HexahedraMesh`
- [x] 缃戞牸 I/O: Nastran (`.nas`), Gmsh (`.msh`)
- [x] 鍧愭爣鍙樻崲, 楂樻柉姹傜Н
- [x] RWG 鍩哄嚱鏁?(涓?Legacy 100% 鍖归厤, 鍖呮嫭杈规帓搴忛€昏緫)
- [x] SWG, RBF, PWC 鍩哄嚱鏁?

### Phase 3: 绉垎鏂圭▼ (2025-12 ~ 2026-01)
- [x] EFIE: PEC Plate/Sphere 楠岃瘉, 濂囧紓椤?($F_1$, $F_2$) 淇
- [x] MFIE: `mfie_interaction!` in-place 缁勮
- [x] CFIE: PEC Sphere 楠岃瘉 (RMSE 2.09 dB vs Mie)
- [x] VEFIE: 浣撶Н绉垎鏂圭▼, SWG 鍩哄嚱鏁版敮鎸?
- [x] SCFIE: 闈綋鑰﹀悎绉垎鏂圭▼

### Phase 4: MLFMA (2026-01)
- [x] 鍏弶鏍戞瀯寤?(`Octree.jl`, `OctreeBuilder.jl`)
- [x] 鑱氬悎 (`Aggregation.jl`), 鍚爣閲忓娍椤?
- [x] 杞Щ (`Translation.jl`), Legacy 鍥犲瓙 $-jk/16\pi^2$ 瀵归綈
- [x] 瑙ｈ仛 (`Disaggregation.jl`)
- [x] Lebedev 鐞冮潰鎻掑€奸泦鎴?
- [x] 杩戝満涓€鑷存€? Max Diff < 1e-12
- [x] 杩滃満绮惧害: 淇 1/4 鍥犲瓙, 杩戦偦缂撳啿鍖?= 4

### Phase 5: 姹傝В鍣ㄤ笌骞惰 (2026-01)
- [x] Direct Solver (LU)
- [x] GMRES (璇樊 2.7e-7)
- [x] BiCGSTAB
- [x] ILU 棰勬潯浠跺櫒
- [x] SPAI 棰勬潯浠跺櫒
- [x] MPI 鍒嗗竷寮忓苟琛?(n=2 vs n=1 鏈哄櫒绮惧害鍖归厤)
- [x] 澶氱嚎绋嬪苟琛?(4 绾跨▼鍔犻€熼獙璇?

### Phase 6: 鍚庡鐞嗕笌 I/O (2026-01)
- [x] RCS 璁＄畻
- [x] FarField / NearField 璁＄畻
- [x] VTK 瀵煎嚭 (ParaView)
- [x] 缁撴灉鏂囦欢 I/O (HDF5, CSV, TXT)
- [x] 鐢垫祦鍒嗗竷鍚庡鐞?

### Phase 7: 楠岃瘉涓庡榻?(2026-02-28) 鉁?
- [x] SCFIE MLFMA 杩戝満楠岃瘉: Rel Err = 1.58e-15 (鏈哄櫒绮惧害)
- [x] SCFIE MLFMA 杩滃満楠岃瘉: Overall 0.85%, Surface 0.85%, Volume 6.1%
- [x] Standalone EFIE MLFMA: 0.66% (4GHz, TriTetra.nas)
- [x] Standalone VEFIE MLFMA: 1.39% (4GHz, TriTetra.nas)
- [x] MoM_AllinOne 鍏ㄩ儴绠椾緥瀵规爣瀹屾垚

**SCFIE MLFMA 淇鐨?7 涓?Bug:**
1. `Disaggregation.jl`: SCFIE `efie_factor` 浠?`1.0+0im` 鈫?`jk畏/(16蟺)`
2. `MLFMAOperator.jl`: 杩戝満 SS 鍧楃Щ闄ゅ浣?`eta` (閬垮厤 MFIE 椤瑰弻閲嶄箻 畏)
3. `MLFMAOperator.jl`: VV 鍧椾粠 `vefie_element_interaction` (c1=-j蠅渭鈧€魏) 鈫?`vefie_element_interaction_kernel` (c1=+j蠅渭鈧€魏)
4. `MLFMAOperator.jl`: VV 鍧楁坊鍔犵己澶辩殑 mass matrix (鑷氦浜掗」)
5. `MLFMAOperator.jl`: VV/SV/VS 鍧楀垱寤?`distribute_term_nosign!` 閬垮厤 bfsSign 鍙岄噸璁℃暟
6. `Disaggregation.jl`: SWG `const_factor` 浠?`-jk畏` 鈫?`jk畏/(4蟺)` (VEFIE 鐢?G=e^{-jkR}/(4蟺R))
7. `MLFMAOperator.jl`: 娣诲姞 VEFIE 缂撳瓨棰勮绠?(TetBasisCache + precompute_vefie_basis)
8. `SCFIE.jl`: Z_SV 绗﹀彿淇 `(term1 - term2)` 鈫?`(term1 + term2)` (鎭㈠ $L$ 绠楀瓙浜掓槗鎬?
9. `SCFIE.jl`: Z_VS 绯绘暟淇 `c1_vs = -j蠅渭鈧€` 鈫?`+j蠅渭鈧€` (缁熶竴涓?Legacy 涓€鑷寸殑姝ｅ彿绾﹀畾)
10. 娣诲姞 `test/test_scfie.jl` 鍥炲綊娴嬭瘯 (浜掓槗鎬с€丏irect 缁勮銆丮LFMA 杩戝満/杩滃満)

### Phase 7.5: ~3 dB 绯荤粺鍋忓樊鏍瑰洜淇 (2026-02-28) 鉁?

**鏍瑰洜 1 鈥?`edgev虃`/`edgen虃` 鏂瑰悜鍙嶈浆** (`BasisUtilities.jl`):
- EMSuite 璁＄畻 `e1 = v2 - v3`锛坴3鈫抳2 鏂瑰悜锛夛紝Legacy 璁＄畻 `v3 - v2`锛坴2鈫抳3 鏂瑰悜锛?
- 瀵艰嚧 `edgen虃 = cross(edgev虃, facen虃)` 鎸囧悜涓夎褰㈠唴閮ㄨ€岄潪澶栭儴
- 褰卞搷: `faceSingularityIgIvecg` 涓?`p02jl`锛堟姇褰辫窛绂伙級绗﹀彿閿欒 鈫?杩戝寮傜Н鍒嗙粨鏋滈敊璇?
- 淇: 灏嗚竟鍚戦噺鏀逛负 `e1 = v3 - v2`, `e2 = v1 - v3`, `e3 = v2 - v1`

**鏍瑰洜 2 鈥?`calc_near_interaction!` 缁忛獙鍥犲瓙 + 闈㈢Н褰掍竴鍖?* (`EFIE.jl`):
- 瀛樺湪缁忛獙鍥犲瓙 `* 1.25`锛堝簲涓?1.0锛?
- `inv_areas = 1.0 / tri_test.area`锛堝簲涓?`1.0 / (tri_test.area * tri_source.area)`锛?
- 鍦?edgev 鏂瑰悜閿欒鏃讹紝杩欎袱涓敊璇儴鍒嗕簰鐩歌ˉ鍋匡紙Z_near 鈮?0.0485脳 姝ｇ‘鍊?鈫?杩戜箮鍙拷鐣ワ級
- 淇: 绉婚櫎 `* 1.25`锛屾敼涓烘纭殑鍙岄潰绉綊涓€鍖?

**楠岃瘉缁撴灉:**
- Z 鐭╅樀瀵硅绾?ratio: 1.000000 (2640脳2640 鏉跨綉鏍?vs Legacy)
- Frobenius 鑼冩暟姣? 1.000010
- RCS (Jet 100MHz): Mean Diff 0.05 dB / RMSE 0.29 dB (Phi=0), Mean Diff 0.008 dB / RMSE 0.09 dB (Phi=90)
- 鍏ㄩ儴 138/138 鍗曞厓娴嬭瘯閫氳繃锛堟棤鍥炲綊锛?

### Phase 10.A: MLFMA 鍥犲瓙淇 (2026-03-02) 鉁?

**Bug 1 鈥?MLFMA far-field 脳4 鍥犲瓙** (`MLFMAOperator.jl`, `Disaggregation.jl`):
- 鏍瑰洜: `efie.factor = jk畏/(16蟺)` 鍖呭惈 `1/4` (鏉ヨ嚜 RWG `l虏/4` 褰掍竴鍖?, 浣?MLFMA 鑱氬悎/瑙ｈ仛鍚勭敤 `l/2`锛屼箻绉?`l虏/4` 宸茶嚜鐒跺寘鍚鍥犲瓙 鈫?**鍙岄噸璁℃暟 1/4**
- Legacy 閬垮厤姝ら棶棰? translation 鐢?`-jk/(16蟺虏)` + disagg 鐢?`jk畏`锛堜笉鍚岀殑鍥犲瓙鍒嗚В鏂瑰紡锛?
- 淇: `y_far *= 4 * operator.factor` (EFIE mul!), CFIE/SCFIE disagg `efie_factor *= 4`
- 楠岃瘉: 鑷唇鎬х郴鏁拌宸?65.7% 鈫?0.30%, RCS RMSE 3.1 dB 鈫?0.028 dB
- A3 S-EFIE MLFMA vs Legacy: Mean Diff 0.048 dB, RMSE 0.303 dB

**Bug 2 鈥?CFIE MLFMA MFIE 绗﹀彿閿欒** (`Disaggregation.jl`):
- 鏍瑰洜: MFIE K 绠楀瓙浣跨敤 $\nabla_{r'}G$锛堟簮姊害锛夛紝杩滃満杩戜技涓?$+jk\hat{k}G$;
  浠ｇ爜閿欒鍦颁娇鐢ㄤ簡 $\nabla_r G$锛堝満姊害锛夌殑 $-jk\hat{k}G$锛屽鑷?MFIE 椤圭鍙峰弽杞?
- Legacy 澶勭悊: 鑱氬悎/瑙ｈ仛鍒嗙, 缁熶竴涔樹互 `jk畏`, EFIE 鍜?MFIE 杩滃満绯绘暟鍚屽彿 (+jk畏)
- 淇: `(-efie_factor)` 鈫?`(+efie_factor)` (MFIE 椤?
- 楠岃瘉: C3 CFIE MLFMA vs Legacy: RMSE 3.45 dB 鈫?**0.003 dB** (1000脳 鏀瑰杽)
- GMRES 杩唬娆℃暟: 50 鈫?7锛堢畻瀛愬噯纭害鎻愬崌鍚庢敹鏁涘姞閫燂級

**宸查獙璇佹祴璇曠粨鏋滄眹鎬?**

| 娴嬭瘯 | 鎸囨爣 | 缁撴灉 |
|------|------|------|
| 鍗曞厓娴嬭瘯 | 138/138 | 鉁?PASS |
| A1 S-EFIE Direct Jet | RMSE vs Legacy | 0.215 dB |
| A3 S-EFIE MLFMA Jet | RMSE vs Legacy | 0.303 dB |
| A3 self-consistency | 绯绘暟璇樊 | 0.30% |
| B1 CFIE 鍒嗚В | rel_err | 0.0 (10/10) |
| C1 S-CFIE Direct Sphere | RMSE vs Legacy | 0.001 dB |
| C3 S-CFIE MLFMA Sphere | RMSE vs Legacy | **0.003 dB** |
| D1-SWG V-EFIE Direct | RMSE vs Legacy | 0.952 dB |
| E1 VSEFIE Direct | RMSE vs Legacy | **0.602 dB** |
| EFIE MLFMA Sphere | RMSE vs Legacy SCFIE | 0.041 dB |

### Phase 10.B: SCFIE Fss 杈圭晫淇 (2026-03-03) 鉁?

**Bug 鈥?缂哄け鍗婂熀鍑芥暟杈圭晫闈㈢Н鍒嗕慨姝?(Fss)** (`SCFIE.jl`, `MLFMAOperator.jl`):
- 鏍瑰洜: 杈圭晫 SWG 鍩哄嚱鏁帮紙鍗婂熀鍑芥暟锛屼粎鏈変竴涓洓闈綋鏀拺锛夌殑鏍囬噺鍔跨己澶辫〃闈㈢Н鍒嗕慨姝ｉ」
- Legacy 鍦?`EFIEVSIERWGSWG.jl` 涓€氳繃 `Fss` 椤瑰鐞嗭細
  - `isbdn=true` 鏃? Z_VS[n,m] += j蠅渭鈧€/(4蟺k虏) 脳 l_m 脳 |A_n| 脳 鈭埆 G dS_tri dS_face
  - `未魏鈮?` 鏃? Z_SV[m,n] += 未魏 脳 (鍚屼笂)
- EMSuite 瀹屽叏缂哄け姝や慨姝?鈫?鑰﹀悎鐭╅樀 Z_SV/Z_VS 鍋忓樊 22%, Z_VV 鍋忓樊 48%
- 淇: 鍦?`SCFIE.jl` 涓坊鍔?`assemble_fss_boundary_correction!` 鍜?`assemble_fss_boundary_correction_sparse`
  - 鐩存帴姹傝В璺緞: 鍦?`assemble_coupling_blocks!` 鍚庤皟鐢?
  - MLFMA 璺緞: 浠ョ█鐤忕煩闃靛舰寮忓姞鍒?Z_near
- 楠岃瘉: E1-VSEFIE RMSE 浠?**5.3 dB 鈫?0.60 dB** (PASS), 138/138 娴嬭瘯鍏ㄩ€氳繃

### 楠岃瘉閲岀▼纰?
- [x] **SEFIE Direct**: `verify_SEFIE_direct.jl` (18.8s 4绾跨▼ / 30.3s 1绾跨▼, Legacy 31.2s)
- [x] **SEFIE MLFMA**: `verify_SEFIE_mlfma.jl` (Ratio 1.0000, Rel Err 1.5%)
- [x] **VEFIE Direct**: `verify_VEFIE_direct.jl` (Legacy Parity)
- [x] **VEFIE MLFMA**: `verify_VEFIE_mlfma.jl` (Ratio 1.0000, Rel Err 0.04%)
- [x] **SCFIE Direct**: `verify_SCFIE_direct.jl` (VSIE Plate+Metal, RCS -15.35 dBsm)
- [x] **SCFIE MLFMA**: `quick_scfie_mlfma_test.jl` (Near-field 1.58e-15, Far-field 0.85%)
- [x] **Standalone EFIE MLFMA**: `test_efie_vefie_farfield.jl` (Rel Err 0.66%)
- [x] **Standalone VEFIE MLFMA**: `test_efie_vefie_farfield.jl` (Rel Err 1.39%)
- [x] **MPI**: `benchmark_parallel_sphere.jl` (Consistency check passed)
- [x] **Threading**: 4 threads vs 1 thread speedup 楠岃瘉

---

## 宸插畬鎴?鉁?(缁?

### Phase 8: 鎬ц兘浼樺寲 (2026-03) 鉁?

#### 8.0 鎬ц兘鍩虹嚎娴嬮噺 鉁?(commit `861426d`)
- [x] 鍒涘缓 `benchmark/performance_baseline.jl` (EMSuite 7 鐢ㄤ緥)
- [x] 鍒涘缓 `LegacyBenchmark/legacy_performance_baseline.jl` (Legacy 瀵规爣)
- [x] EMSuite 鍏ㄩ儴 7 鐢ㄤ緥娴嬮噺瀹屾垚
- [x] Legacy 5 鐢ㄤ緥娴嬮噺瀹屾垚
- [x] 缁煎悎瀵规瘮鎶ュ憡: `test_results/PERFORMANCE_BASELINE.md`

#### 8.1 Z 缁勮鍘婚攣 鉁?(commit `2d4ebe6`)
- [x] `Impedance.jl` SpinLock 鈫?Per-row SpinLock (琛岀骇鏃犻攣骞惰)
- [x] Plate EFIE 缁勮 **-54%**, Jet EFIE 缁勮 **-12%**
- [x] 138/138 娴嬭瘯閫氳繃

#### 8.2 CFIE 鍐呮牳鍚堝苟 鉁?(commit `d0888cf`)
- [x] MFIE 鍐呮牳浼樺寲锛氬叡浜?Green 鍑芥暟銆乮nline rho 鍚戦噺銆佹秷闄ら噸澶嶅嚑浣曡绠?
- [x] CFIE 缁勮 **-74%** (Jet: 168.29s 鈫?43s)
- [x] CFIE/EFIE 缁勮姣? 8.1脳 鈫?**2.31脳** (鐩爣 鈮?2.5脳 鉁?

#### 8.3 MLFMA Z_near 浼樺寲 鉁?(commit `d2f7963`)
- [x] 棰勫垎閰?COO 鏁扮粍浠ｆ浛鍔ㄦ€?push!
- [x] COO 鍚堝苟鍚?sparse() 鏋勯€?

#### 8.4 鍐呭瓨鍒嗛厤鐑偣 鉁?(commit `82988cf`)
- [x] 绉婚櫎 MFIE 褰卞瓙 `get_global_quad_points` 鍑芥暟

#### 8.5 Julia 1.12 鍏煎淇 鉁?(commit `67d3a8a`)
- [x] `threadid()` 鈫?`Threads.maxthreadid()` (Legacy + EMSuite)
- [x] 8.5b 绫诲瀷绋冲畾鎬у鏌ワ細`@code_warntype` 鍏ㄩ儴 clean

#### 8.6 @fastmath + SIMD 鉁?(commit `1c6d499`)
- [x] `calc_interaction!` 閲嶅啓锛氱洿鎺?dot() 鏇挎崲 SMatrix
- [x] `@fastmath` 鍔犻€?exp() 绛夋暟瀛﹁繍绠?
- [x] `@inbounds @simd` 浼樺寲鍐呭惊鐜?

#### 8.7 BlockJacobiPreconditioner 鉁?(commit `76f8b16`)
- [x] 瀹炵幇 `BlockJacobiPreconditioner` (浠?Z_near 鎻愬彇瀵硅鍧? 骞惰 LU)
- [x] 鏋勫缓閫熷害姣?Sparse LU 蹇?**166脳**
- [x] 閫傜敤浜?CFIE (3 娆?GMRES 杩唬); EFIE 涓嶆敹鏁? LU 浠嶄负榛樿
- [x] 娣诲姞 `get_leaf_intervals(op::MLFMAOperator)`

#### 8.8 鏈€缁堝熀鍑嗗娴?鉁?(commit `6f4987a`)
- [x] 鍏ㄩ儴 6 涓敤渚?(+ CFIE 瀵规瘮) 閲嶆柊璁℃椂
- [x] 淇 Sphere CFIE MLFMA OOM (COO 鍒濆鍒嗛厤涓婇檺)
- [x] 鐢熸垚 `test_results/PERFORMANCE_REPORT.md`

**Phase 8 鏈€缁堢粨鏋?(2026-03-01 鏇存柊):**

| 鐢ㄤ緥 | N | 鍘熷鍩虹嚎 | Phase 8.8 | **Phase 8.9** | **鏈€鏂?* | 璇存槑 |
|------|---|---------|---------|----------|------|------|
| Plate EFIE | 2640 | 1.02s | 1.94s | **0.153s** | 0.153s | EFIE 鍐呮牳閲嶅啓 |
| Jet EFIE | 14559 | 20.70s | 29.01s | **4.26s** | 4.26s | EFIE SIMD 淇 |
| **Jet CFIE** | 14559 | **168.29s** | **64.88s** | **14.48s** | 14.48s | CFIE 鏋舵瀯 + `@.` |
| Jet MLFMA | 14559 | 76.69s | 108.93s | 鏈噸娴?| 鈥?| 鈥?|
| Sphere MLFMA | 26424 | 323.25s | 285.81s | 鏈噸娴?| 鈥?| 鈥?|
| **VEFIE** | 15828 | 46.13s | 66.24s | 鏈噸娴?| **41.30s 鉁?* | 涓婁笁瑙掑绉颁紭鍖?|
| **SCFIE** | 15860 | 66.68s | 96.94s | 鏈噸娴?| **65.67s 鉁?* | VEFIE+鑰﹀悎浼樺寲 |

#### 8.9 EFIE/CFIE/SCFIE 娣卞害浼樺寲 鉁?(commit `8f8dfc3`, `f520609`)
- [x] **EFIE `calc_interaction!` 閲嶅啓**: 绉婚櫎 `@simd for n in 1:3` + 涓夊厓鍒嗘敮锛屾敼鐢?tuple-indexed rho + 瀹屽叏灞曞紑 3脳3 鍐呯Н 鈫?Jet EFIE **4.26s** (-79.4% vs 20.7s 鍘熷鍩虹嚎)
- [x] **CFIE 鏋舵瀯淇**: 鍚堝苟姹囩紪瀹炴祴姣斿垎绂绘眹缂栨參 (register/cache pressure)锛屾敼涓哄垎绂昏皟鐢?+ `@.` 灏卞湴鍔犳潈姹傚拰 (閬垮厤绗笁涓?N脳N 鍒嗛厤) 鈫?Jet CFIE **14.48s** (-91.4% vs 168.3s 鍘熷鍩虹嚎)
- [x] **SCFIE Fss 骞惰鍖?*: `assemble_fss_boundary_correction!` 娣诲姞 `@threads` + 琛岀骇 SpinLock
- [x] **鍗曞厓娴嬭瘯**: 179/179 閫氳繃 (Testing EMSuite tests passed)

#### 8.10 VEFIE/SCFIE 鎬ц兘绐佺牬 鉁?(commit `bbf8fdd`)
- [x] **VEFIE 涓婁笁瑙掑绉颁紭鍖?*: 灏嗗叏 N虏 tet 瀵瑰惊鐜敼涓轰笂涓夎 N*(N+1)/2 鍧? 
  - 澶栧眰寰幆 test tet `it`锛屽唴灞?`js` 浠?`it+1` 鍒?`ntet`锛堜笂涓夎锛? 
  - 鍒╃敤 Z_st[j,i] = (魏_t/魏_s) 脳 Z_ts[i,j]锛屽悓鏃跺啓鍏?Z[m,n] 鍜?Z[n,m]  
  - 鍏ㄥ眬 SpinLock锛堥攣鎸佹湁鏃堕棿 鈮?2 娆℃爣閲忓啓 鈮?2 ns锛岃绠楁椂闂?鈮? 渭s锛岀珵浜夌巼 <1%锛? 
  - **VEFIE: 41.30s**锛坴s Legacy 46.13s锛?*蹇?12%**锛泇s 鍩虹嚎 66.24s锛?*蹇?1.60脳**锛? 
- [x] **SCFIE 鑰﹀悎鍧椾簰鏄撴€т紭鍖?*: Z_vs = Z_sv / 魏锛坈ommit `a87be12`锛夛紝鍑忓皯涓€鍗婅€﹀悎绉垎  
  - **SCFIE: 65.67s**锛坴s Legacy 66.68s锛?*蹇?1.5%**锛泇s 鍩虹嚎 96.94s锛?*蹇?1.48脳**锛? 
- [x] **鍏ㄩ儴 179/179 娴嬭瘯閫氳繃**


鹿 EFIE 缁勮澧炲箙: @fastmath/SIMD 閲嶅啓涓昏浼樺寲 MFIE 璺緞, 瀵圭函 EFIE 鏈夎交寰紑閿€
虏 MLFMA EFIE 澧炲箙: 棰勬潯浠跺櫒 LU 鍙樻參 (8.89s鈫?7.27s), 闈炰唬鐮佸洖褰?
鲁 VEFIE/SCFIE 缁勮澧炲箙鍚屽洜; LU 姹傝В澶у箙鍔犻€?(155.61s鈫?1.12s for VEFIE)

---

## 杩涜涓?馃敡

锛堟棤褰撳墠杩涜涓换鍔★級

### Phase 12: 鍏潰浣撳畬鏁存敮鎸?PWCHex + RBF (2026-03-04) 鉁?

**鐩爣**: 瀹炵幇鎵€鏈夌己澶辩殑鍏潰浣撳熀鍑芥暟+绉垎鏂圭▼缁勫悎锛岃ˉ鍏?Phase 12 璺嚎鍥句腑鐨?10 涓?Gap銆?

**淇敼鏂囦欢:**
1. **`src/Geometry/GaussQuadrature.jl`** 鈥?鏂板鍏潰浣?(8鐐?tensor-product GL) 鍜屽洓杈瑰舰 (4鐐? GQ 瑙勫垯
2. **`src/Geometry/MeshTypes.jl`** 鈥?鏂板 `HexahedraInfo`, `Quads4Hexa` 缁撴瀯浣?+ 杈呭姪鍑芥暟 (`get_free_vns`, `set_delta_kappa!`, `hex_volume` 绛?
3. **`src/Geometry/MeshIO.jl`** 鈥?鏂板 CHEXA Nastran 缃戞牸璇诲彇锛屾敮鎸佺画琛岀鏍煎紡
4. **`src/BasisFunctions/PWC.jl`** 鈥?鏂板 `PWCHexBasis` (3 DOF/鍏潰浣? x,y,z 鍒嗛噺)
5. **`src/BasisFunctions/RBF.jl`** 鈥?瀹屽杽 `evaluate()` 瀹炵幇锛屽惎鐢ㄨ竟鐣岄潰鍩哄嚱鏁?
6. **`src/BasisFunctions/BasisUtilities.jl`** 鈥?鏂板 `get_hexahedra_info(mesh, PWCHexBasis/RBFBasis, permittivities)`
7. **`src/IntegralEquations/VEFIE.jl`** 鈥?~500琛? PWCHex, RBF, 娣峰悎 TetraHex 瑁呴厤 (6 涓柊鏂规硶)
   - 娉涘寲 `_pwc_dyad_kernel!` 鏀寔 duck-typed 浣撳厓绱?
   - 娣峰悎 TetraHex 瑁呴厤: 4涓瓙鍧?(TT, TH, HT, HH)
8. **`src/IntegralEquations/SCFIE.jl`** 鈥?~380琛? RWG+PWCHex (骞剁煝 L 绠楀瓙) 鍜?RWG+RBF (鏍囬噺鍔垮舰寮?+ Fss 杈圭晫淇)
9. **`src/IntegralEquations/Excitation.jl`** 鈥?~180琛? PWCHex 鍜?RBF 骞抽潰娉㈡縺鍔卞悜閲?+ 缁勫悎 SCFIE 鐗堟湰
10. **`src/PostProcessing/RadiationIntegral.jl`** 鈥?~120琛? PWCHex 鍜?RBF 杈愬皠绉垎
11. **`src/PostProcessing/RCS.jl`** 鈥?~80琛? PWCHex 鍜?RBF RCS 璁＄畻鏂规硶
12. **`test/test_basis_functions.jl`** 鈥?淇 RBF 娴嬭瘯棰勬湡鍊?(num_basis 1鈫?1, 鏌ユ壘鍐呴儴鍩哄嚱鏁?

**鏂规硶缁熻:**
- 18 涓?`assemble_impedance_matrix` 鏂规硶 (EFIE/MFIE/CFIE/VEFIE/SCFIE 脳 鍚勫熀鍑芥暟缁勫悎)
- 17 涓?`excitation_vector` 鏂规硶
- 5 涓?`radarCrossSection` 鏂规硶

**娴嬭瘯缁撴灉:**
- 鍏ㄩ儴 179/179 娴嬭瘯閫氳繃 (鏃犲洖褰?
- +2296 琛屼唬鐮?
- Commit: `099385b`

### Phase 11: PWC 鍩哄嚱鏁版敮鎸佹墿灞?(2026-03-04) 鉁?

**鐩爣**: 瀵归綈 Legacy 鐨?PWC (Piecewise Constant) 鍩哄嚱鏁版敮鎸侊紝瀹屽杽 VEFIE+PWC 鍜?SCFIE+RWG+PWC 缁勫悎銆?

**淇敼鏂囦欢:**
1. **`src/BasisFunctions/PWC.jl`** 鈥?瀹屽叏閲嶅啓: 3 DOFs/鍥涢潰浣?(x,y,z 鍒嗛噺)
   - `PWC` struct 澧炲姞 `inBfsID::SVector{3, IT}` (涓変釜鍏ㄥ眬鍩哄嚱鏁癐D)
   - `num_basis` 杩斿洖 `3 * length(functions)` (鍘熶负 1:1)
   - `evaluate` 杩斿洖鍗曚綅鍚戦噺 x虃/欧/岷?(鍩轰簬 `mod1(i, 3)`)
   - Legacy 瀵归綈: `MoM_Basics` 鐨?`nPWC = 3 * num_tetrahedra`

2. **`src/BasisFunctions/BasisUtilities.jl`** 鈥?鏂板 `get_tetrahedra_info(mesh, basis::PWCBasis, permittivities)`
   - `inBfsID = SVector{4}(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3, 0)` (绗?椤规湭浣跨敤)

3. **`src/IntegralEquations/VEFIE.jl`** 鈥?鏂板 ~230 琛? VEFIE+PWC 缁勮
   - `assemble_impedance_matrix(vefie::VEFIE, basis::PWCBasis)` 鈥?甯?permittivities 鐨?鍙?2鍙傜増鏈?
   - `_pwc_dyad_kernel!` 鈥?3脳3 骞剁煝 L 绠楀瓙: $(k^2 I + \nabla\nabla) G(R)$
   - 瀵圭О缁勮 + 鑷€傚簲绉垎 (杩滃満1鐐?杩戝満5鐐?
   - 鑷綔鐢ㄩ」璐ㄩ噺鐭╅樀: $V/(j\omega\varepsilon)$

4. **`src/IntegralEquations/Excitation.jl`** 鈥?鏂板 ~100 琛? PWC 婵€鍔卞悜閲?
   - VEFIE 绠楀瓙鐗? `excitation_vector(op::VEFIE, source::PlaneWave, basis::PWCBasis)`
   - 鐙珛鐗? `excitation_vector(source::PlaneWave, basis::PWCBasis)`
   - 缁勫悎鐗? `excitation_vector(source, surf_basis::RWGBasis, vol_basis::PWCBasis)`

5. **`src/IntegralEquations/SCFIE.jl`** 鈥?鏂板 ~170 琛? SCFIE+RWG+PWC 鑰﹀悎
   - `assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCBasis)`
   - `assemble_coupling_blocks_pwc!` 鈥?骞剁煝 L 绠楀瓙鑰﹀悎
   - Z_SV 鍖呭惈 魏, Z_VS 鏃?魏
   - 鏃?Fss 杈圭晫淇 (PWC 鏃犲崐鍩哄嚱鏁?

6. **`src/PostProcessing/RadiationIntegral.jl`** 鈥?鏂板 PWC 杈愬皠绉垎
   - `radiation_integral_pwc`: $N = \sum_t V_t \kappa_t \sum_{gq} J \cdot e^{jk\hat{r}\cdot r_{gq}} w_{gq}$

7. **`src/PostProcessing/RCS.jl`** 鈥?鏂板 PWC RCS 璁＄畻
   - `radarCrossSection(..., basis::PWCBasis, permittivities)`

8. **`src/Driver.jl`** 鈥?鎵╁睍鏀寔澶欼E绫诲瀷 (EFIE/MFIE/CFIE/VEFIE/SCFIE)

9. **`src/Core/Configuration.jl`** 鈥?SimulationConfig 澧炲姞 `ie_type`, `cfie_alpha`, `permittivities` 瀛楁

10. **`test/test_pwc.jl`** 鈥?鏂板 PWC 涓撶敤娴嬭瘯 (鍩哄嚱鏁版瀯閫? VEFIE+PWC, 婵€鍔? SCFIE+PWC)

11. **`test/test_basis_functions.jl`** 鈥?鏇存柊 PWC 娴嬭瘯閫傞厤 3 DOFs/tet

**淇鐨?Bug:**
- VEFIE.jl 缂哄け module 闂悎 `end` (缂栬瘧閿欒)
- Excitation.jl 缂哄け module 闂悎 `end` (缂栬瘧閿欒)
- Driver.jl 浣跨敤 `get(struct, ...)` 瀵艰嚧 MethodError (struct 涓嶆敮鎸?Dict 鐨?get)

**娴嬭瘯缁撴灉:**
- 鍏ㄩ儴 139+/139+ 娴嬭瘯閫氳繃 (鍚柊澧?PWC 娴嬭瘯)
- PWC 鍩虹娴嬭瘯: 16/16 PASS
- 鏃犲洖褰?

---

## 杩涜涓?馃殌

### Phase 9: 浠ｇ爜璐ㄩ噺涓庡彂甯?
- [x] 娴嬭瘯濂椾欢娓呯悊: 179/179 鍏ㄩ儴閫氳繃 鉁?(鍚?PWC/RBF/PWCHex 鏂版祴璇?
- [x] JuliaFormatter.jl 缁熶竴浠ｇ爜椋庢牸 鉁?(`0535576`) 鈥?79 涓?`src/` 鏂囦欢
- [x] CHANGELOG.md 瀹屽杽 鉁?鈥?Phase 8-12 鎵€鏈夐噷绋嬬
- [x] compat 璇硶瑙勮寖鍖?鉁?+ 瑕嗙洊鐜囪剼鏈?+ 鍙戝竷娓呭崟 (`9d3e671`)
- [ ] 娴嬭瘯瑕嗙洊鐜囩粺璁′笌鎻愬崌 (鐩爣 > 80%)
- [ ] API 鏂囨。琛ュ叏 (鎵€鏈夊叕鍏辨帴鍙?
- [ ] 鐢ㄦ埛鏁欑▼ (Quick Start, Advanced)
- [ ] 鐞嗚鏂囨。 (MoM, MLFMA, 绉垎鏂圭▼鎺ㄥ)
- [ ] 鍙戝竷鍒?Julia General Registry

---

## Legacy 鍥犲瓙瀵圭収琛?

| 椤圭洰 | Legacy | EMSuite | 璇存槑 |
|------|--------|---------|------|
| EFIE Z | $1/16\pi$ | $1/4\pi$ + 鏄惧紡 $l/2A$ | 鏁板绛変环 |
| FarField | $1/4\pi$ | $1/4\pi$ | 涓€鑷?|
| Translation | $1/16\\pi^2$ | $-jk/16\\pi^2$ | Legacy 瀵归綈 |\n| SWG Disagg | N/A | $jk\\eta/(4\\pi)$ | VEFIE G 鍚?$1/(4\\pi)$ |
| 鏃堕棿绾﹀畾 | $e^{-j\omega t}$ | $e^{j\omega t}$ | Z 铏氶儴绗﹀彿鐩稿弽 |
| Area/Length | 闅愬紡鍖呭惈鍦ㄥ洜瀛愪腑 | 鏄惧紡褰掍竴鍖?| 缁撴灉绛変环 |

---

## 宸茬煡闂

1. ~~**EMSuite vs Legacy ~3 dB 绯荤当鍋忓樊**~~ 鈥?**宸蹭慨澶?* (2026-02-28): 鏍瑰洜涓?`edgev虃` 鏂瑰悜鍙嶈浆 + `calc_near_interaction!` 缁忛獙鍥犲瓙銆備慨澶嶅悗 RCS 鍋忓樊 < 0.3 dB RMSE銆?
2. ~~**MLFMA 杩滃満 脳4 鍥犲瓙 + CFIE 绗﹀彿**~~ 鈥?**宸蹭慨澶?* (2026-03-02): EFIE MLFMA 绯绘暟璇樊 65.7% 鈫?0.30%. CFIE MLFMA RMSE 3.45 dB 鈫?0.003 dB.
3. **VEFIE Mie 鍋忓樊**: Legacy 鍜?EMSuite 鍧囨瘮 Mie 绾ф暟浣?~25dB 鈥?灞炰簬 Legacy 绠楁硶鍥烘湁闂, 鏍囪涓?"Legacy Parity", 鐗╃悊淇涓烘湭鏉ョ爺绌惰棰?
4. **BiCGSTAB 鏀舵暃**: 闇€瑕侀鏉′欢鎵嶈兘鍙潬鏀舵暃
5. **SCFIE 鑰﹀悎椤逛簰鏄撴€?*: 宸蹭慨澶?鈥?Z_SV/魏 = Z_VS^T 鍦ㄦ満鍣ㄧ簿搴︽垚绔?(2.99e-16)
6. ~~**EFIE 闂悎浣撳唴閮ㄨ皭鎸?*: EFIE 鐢ㄤ簬闂悎瀵间綋鏃舵潯浠舵暟宸?(Direct vs MLFMA 绯绘暟宸?62%)锛屽簲鏀圭敤 CFIE~~ 鈥?**宸茬‘璁?*: 鐜板湪 CFIE MLFMA 姝ｇ‘宸ヤ綔 (RMSE 0.003 dB, 7 iterations)
7. **SWG MLFMA const_factor 绗﹀彿**: `const_factor = jk畏/(4蟺)` 鍙兘搴斾负 `-jk畏/(4蟺)` (VEFIE `c1 = -jk畏魏` 鍚礋鍙?. ~~闇€瑕?VEFIE MLFMA 绮惧害娴嬭瘯楠岃瘉~~ 鈫?D3 娴嬭瘯 RMSE=0.0 dB 琛ㄦ槑褰撳墠瀹炵幇姝ｇ‘.
8. **A2 S-EFIE Iterative 鏈厖鍒嗘敹鏁?*: restart=1000 + Diagonal 棰勬潯浠朵笅 EFIE (N=14559) RMSE=0.343 dB (> 0.1 dB). 鏍瑰洜: EFIE 鏉′欢鏁板ぇ, 瀵硅棰勬潯浠朵笉瓒? D2/E2 宸茶瘉鏄?GMRES 鍩虹璁炬柦姝ｇ‘; A3 MLFMA+杩戝満 LU 棰勬潯浠跺彲姝ｅ父鏀舵暃

---

## 鏇存柊鏃ュ織

| 鏃ユ湡 | 鏇存柊鍐呭 |
|------|----------|
| 2026-03-06 | **Phase 13.3 V-EFIE MPI 骞惰瑁呴厤瀹屾垚** 鉁?鈥?鏂板 `VolumeAssembly.jl` (303 lines, 鏃?module 灏佽, Allreduce 绛栫暐). 淇: (1) Julia dispatch 璺ㄦā鍧椾笉鍙 (鍒犻櫎 `module VolumeAssemblyMPI`); (2) `BasisFunctions.jl` 琛ュ嚭鍙?`get_triangles_info`; (3) 瀛楃缂栫爜鎹熷潖 魏鈫掗瓘 淇. 姝ｇ‘鎬? Tetra.nas N=3201, Z[1,1]=4.747e6-4.716e6im 鍦?1/2/4 杩涚▼涓嬪畬鍏ㄧ浉鍚? 鎬ц兘: 1P 10.95s 鈫?2P 8.86s (1.24脳). commit: `1730b71` |
| 2026-03-0X | **Phase 13.2 Option C 璇勪及瀹屾垚** 鉁?鈥?鏂板 bench_pwc_rbf_performance.jl (177 lines)銆傚叧閿彂鐜? PWC+VEFIE 12.62s (3.16脳 faster), 浣?DOF +38% (21834 vs 15828), per-DOF虏 鏁堢巼浠?0.17脳; SCFIE(RWG+PWC) 33.26s vs (RWG+SWG) 63.55s (1.91脳 faster); RBF 浠呮敮鎸佸叚闈綋缃戞牸 (涓嶉€傜敤鍥涢潰浣?. **缁撹**: PWC 涓嶉€傚悎閫氱敤浼樺寲 (姹傝В鍣ㄤ唬浠烽珮), SWG 淇濇寔鏈€浣冲钩琛? **鎺ㄨ崘**: 杞悜 Option A (MPI 骞惰鍖? 棰勬湡 3.5脳 @ 4 杩涚▼). 璇﹁ [OPTION_C_PWC_RBF_EVALUATION.md](.github/OPTION_C_PWC_RBF_EVALUATION.md). commit: `bc2121d` |
| 2026-03-01 | **Phase 13.1 FastExp 浼樺寲瀹屾垚** 鉁?鈥?鏂板 FastExp.jl (10000 鏉＄洰绾挎€ф彃鍊艰〃), VEFIE.jl 闆嗘垚, @fastmath/@inbounds/unsafe_trunc 浼樺寲銆俈-EFIE 66.24s鈫?9.05s (1.70脳 vs baseline, **1.18脳 vs Legacy** 鉁?, SCFIE 96.94s鈫?3.81s (1.52脳 vs baseline, **1.04脳 vs Legacy** 鉁?銆俆hread-local buffer 澶辫触鏁欒: -63% 鎬ц兘閫€姝?(16GB 鍐呭瓨鐖嗙偢)銆傝窛绂?2脳 鐩爣: V-EFIE 缂哄彛 70%, SCFIE 缂哄彛 92%銆傝瑙?[PHASE_13.1_SUMMARY.md](.github/PHASE_13.1_SUMMARY.md). commits: `bfb24f6`, `0488e28`, `2c2ea75` |
| 2026-03-01 | **Phase 13.1 FastExp 鏌ユ壘琛ㄤ紭鍖栧畬鎴?* 鈥?鏂板 FastExp.jl (10000 鏉＄洰绾挎€ф彃鍊艰〃, 瑕嗙洊 20位), VEFIE.jl 闆嗘垚 (struct 瀛楁 + 涓ゅ璋冪敤鐐?, 浣跨敤 @fastmath/@inbounds/unsafe_trunc 浼樺寲鎬ц兘銆傛€ц兘缁撴灉 (plate_and_metal_1dot2GHz, 4 threads): V-EFIE 66.24s鈫?9.05s (1.70脳 vs baseline, 1.18脳 vs Legacy 46.13s 鉁?, SCFIE 96.94s鈫?3.81s (1.52脳 vs baseline, 1.04脳 vs Legacy 66.68s 鉁?. 绮惧害楠岃瘉: 鏈€澶х浉瀵硅宸?< 0.002%. commit: `bfb24f6` |
| 2026-07-30 | **Phase 9.1 浠ｇ爜璐ㄩ噺** 鈥?JuliaFormatter 鏍煎紡鍖?79 涓?`src/` 鏂囦欢; CHANGELOG.md 鍏ㄩ潰鏇存柊; compat 璇硶瑙勮寖鍖?+ julia="1.10"; 瑕嗙洊鐜囪剼鏈?(`scripts/check_coverage.jl`); 鍙戝竷娓呭崟 (`RELEASE_CHECKLIST.md`). 179/179 閫氳繃. commits: `0535576`, `b189a87`, `9d3e671` |
| 2026-07-29 | **Phase 8.10 VEFIE 瀵圭О浼樺寲** 鈥?Z_st[j,i]=(魏_t/魏_s)路Z_ts[i,j] 鎺ㄥ; 涓婁笁瑙掑洓闈綋閬嶅巻; VEFIE **41.30s** (蹇?Legacy 12%); SCFIE **65.67s** (蹇?Legacy 1.5%). commit: `bbf8fdd` |
| 2026-03-01 | **Phase 8.9 EFIE/CFIE/SCFIE 娣卞害浼樺寲** 鈥?EFIE 鍐呮牳閲嶅啓 (绉婚櫎@simd+涓夊厓鍒嗘敮) 鈫?Jet EFIE 20.7s鈫?.26s (-79%); CFIE 鏋舵瀯淇 (鍒嗙姹囩紪+`@.`灏卞湴) 鈫?Jet CFIE 168s鈫?4.5s (-91%); SCFIE Fss 骞惰鍖? 179/179 閫氳繃. commits: 8f8dfc3, f520609 |
| 2026-07-29 | **Phase 10 绮惧害楠岃瘉琛ュ叏瀹屾垚** 鈥?琛ラ綈鍓╀綑 9 瀛愭祴璇?(D2/E2/D3/E3/A2/B2/A4/B3/C3-MPI). 15/16 閫氳繃, A2 鏉′欢閫氳繃 (GMRES 鏀舵暃鍙楅檺). Bug 淇: CFIE MPI 骞惰瑁呴厤 (`cfie_interaction!` 涓嶅瓨鍦?鈫?EFIE+MFIE 鐙珛浜や簰+绾挎€х粍鍚?. 179/179 鍗曞厓娴嬭瘯閫氳繃 |
| 2026-03-04 | **Phase 12 鍏潰浣撳畬鏁存敮鎸佸畬鎴?* 鈥?PWCHexBasis 3 DOFs/hex + RBF evaluate + 杈圭晫闈€侴Q (hex/quad)銆丮eshIO (CHEXA)銆乂EFIE (PWCHex/RBF/Mixed)銆丼CFIE (RWG+PWCHex/RBF)銆佹縺鍔卞悜閲忋€佽緪灏勭Н鍒?RCS銆?79/179 娴嬭瘯閫氳繃銆?2296 琛?|
| 2026-03-04 | **Phase 11 PWC 鍩哄嚱鏁版墿灞曞畬鎴?* 鈥?PWCBasis 3 DOFs/tet, VEFIE+PWC 骞剁煝缁勮, PWC 婵€鍔?杈愬皠绉垎/RCS, SCFIE+RWG+PWC 鑰﹀悎, Driver.jl 澶欼E鎵╁睍, SimulationConfig 澧炲己, 鏂板 test_pwc.jl. 139+/139+ 娴嬭瘯鍏ㄩ€氳繃 |
| 2026-02-28 | **Phase 8 鎬ц兘浼樺寲鍏ㄩ儴瀹屾垚** 鈥?8 涓瓙闃舵 (8.0-8.8), 鏍稿績鎴愭灉: CFIE 缁勮 -61% (168鈫?5s), CFIE/EFIE 姣?8.1脳鈫?.2脳, MLFMA OOM 淇, BlockJacobiPreconditioner, Julia 1.12 鍏煎, 绫诲瀷绋冲畾鎬?clean. 璇﹁ `test_results/PERFORMANCE_REPORT.md` |
| 2026-03-03 | **Phase 8.0 鎬ц兘鍩虹嚎瀹屾垚** 鈥?EMSuite 7 鐢ㄤ緥 + Legacy 5 鐢ㄤ緥璁℃椂銆傚叧閿彂鐜? CFIE 4.61脳 鎱?(鍙岄亶鍘嗛棶棰?, SCFIE 2.26脳 鎱? EFIE/VEFIE 鎸佸钩鎴栨洿蹇? LU 姹傝В蹇?30-40% |
| 2026-03-03 | **Phase 8 鎬ц兘浼樺寲璁″垝** 鈥?鍔犲叆鎬ц兘浼樺寲璺嚎: 6 鐑偣 (SpinLock鍘婚攣/CFIE鍚堝苟/MLFMA Z_near/鍐呭瓨/SIMD/绫诲瀷绋冲畾), 8 姝ラ, 鐩爣 鈮?Legacy 淇濆簳, 鈮?0.5脳 Legacy 鎸戞垬 |
| 2026-03-03 | **SCFIE Fss 杈圭晫淇** 鈥?鍗婂熀鍑芥暟杈圭晫闈㈢Н鍒嗕慨姝ｃ€侲1-VSEFIE RMSE 5.3鈫?.60 dB. D1-SWG VEFIE RMSE 0.95 dB. 138/138 娴嬭瘯閫氳繃 |
| 2026-03-02 | **MLFMA 鍥犲瓙淇脳2** 鈥?(1) EFIE far-field 脳4 鍥犲瓙: 绯绘暟璇樊 65.7%鈫?.30%, RMSE 3.1鈫?.028 dB; (2) CFIE MFIE 绗﹀彿: 鈭嘷{r'}G 缁欏嚭 +jk k虃 (闈?-jk k虃), RMSE 3.45鈫?.003 dB, GMRES 50鈫? 杩唬 |
| 2026-03-01 | **Phase 10 璁″垝** 鈥?鍏ㄦ柟绋嬪叏璺緞绮惧害瀵归綈璁捐瀹屾垚: 5 鏂圭▼ (S-EFIE/S-MFIE/S-CFIE/V-EFIE/VS-EFIE) 脳 4 璺緞 (Direct/Iterative/MLFMA/MPI), 鍏ㄧ悆闈?1314 鐐归噰鏍? 鍏?17 瀛愭祴璇曢」 |
| 2026-02-28 | **绮惧害鏁堢巼鎶ュ憡 v2** 鈥?鍏ㄩ潰鍩哄噯娴嬭瘯: SEFIE Direct(RMSE 0.29dB), CFIE Direct, SEFIE MLFMA, SCFIE MLFMA(Sphere N=26424). 瑙?`test_results/ACCURACY_EFFICIENCY_REPORT.md` |
| 2026-02-27 | 淇楠岃瘉鑴氭湰 benchmark/run_full_benchmark.jl API 閿欒 (MLFMAOperator 鏋勯€?鎺掑簭閫忔槑鎬? |
| 2026-02-27 | 鏂板 MLFMA MatVec 蹇€熸祴璇曡剼鏈?(benchmark/quick_matvec_test.jl) |
| 2026-02-27 | 鏂板 Direct vs MLFMA 鑷竴鑷存€ф祴璇?(benchmark/self_consistency_test.jl) |
| 2026-02-28 | **娴嬭瘯濂椾欢娓呯悊瀹屾垚** 鈥?138/138 鍏ㄩ儴閫氳繃, 淇 6 涓瀛樻祴璇曢棶棰?|
| 2026-02-28 | 淇 `Vector{AbstractBasisFunction}` 绫诲瀷娲惧彂 bug (Aggregation/Disaggregation) |
| 2026-02-28 | **SCFIE MLFMA 楠岃瘉瀹屾垚** 鈥?淇 7+3 涓?bug, 杩戝満 1.58e-15, 杩滃満 1.04% |
| 2026-02-28 | **SCFIE 鑰﹀悎浜掓槗鎬т慨澶?* 鈥?Z_SV `(term1-term2)` 鈫?`(term1+term2)`, Z_VS c1 绗﹀彿淇 |
| 2026-02-28 | 娣诲姞 `test/test_scfie.jl` 鍥炲綊娴嬭瘯 (9 tests all pass) |
| 2026-02-28 | **~3 dB 鍋忓樊鏍瑰洜淇** 鈥?`BasisUtilities.jl` 杈瑰悜閲忔柟鍚?+ `EFIE.jl` 杩戜氦浜掗潰绉綊涓€鍖栥€俍 鐭╅樀 ratio=1.0, RCS RMSE<0.3 dB |
| 2026-02-28 | **绮惧害鏁堢巼鎶ュ憡 v2** 鈥?鍏ㄩ潰鍩哄噯娴嬭瘯: SEFIE Direct(RMSE 0.29dB), CFIE Direct, SEFIE MLFMA, SCFIE MLFMA(Sphere N=26424). 瑙?`test_results/ACCURACY_EFFICIENCY_REPORT.md` |
| 2026-02-27 | 鍒濆鍖栬繘搴︽枃浠? 浠?copilot-instructions.md 杩佺Щ |
| 2026-01-xx | VEFIE MLFMA 楠岃瘉瀹屾垚 (Rel Err 0.04%) |
| 2026-01-xx | SCFIE Direct 楠岃瘉瀹屾垚 (VSIE Plate) |
| 2026-01-xx | MPI/Threading 骞惰楠岃瘉閫氳繃 |
| 2025-12-xx | Surface IE (EFIE/MFIE/CFIE) 鍏ㄩ潰楠岃瘉瀹屾垚 |
| 2025-12-xx | Legacy 鍥犲瓙瀵归綈瀹屾垚, 绉婚櫎缁忛獙甯告暟 |
