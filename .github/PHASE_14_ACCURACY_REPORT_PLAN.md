# Phase 14: 全量精度测试与对比报告

> 创建日期: 2026-03-04  
> 状态: 计划中

---

## 0. 开发原则遵循声明

本 Phase 适用以下核心原则（见 `copilot-instructions.md`）：

- **原则 1 (TDD)**: 先写解析器/读取器的单元测试，再实现
- **原则 5 (Git 提交规范)**: 每个功能模块（Feko 解析器、Mie 参考、报告生成器）独立提交
- **原则 6 (终端输出)**: 不重定向到本地文件；用 Julia 内部 `open(...) do io end` 写报告
- **原则 7 (Phase 结束检视)**: 完成后至少 2 轮 clean 检视再关闭
- **原则 9 (计划文档规范)**: 本文档已包含所有必要节

---

## 1. 目标

对 EMSuite 所有主要积分方程 × 主要求解路径，与 **Feko 商业软件结果** 或 **Mie 解析解** 进行系统对比，输出量化精度报告。

**与 Phase 10 的区别**：

| 维度 | Phase 10 | Phase 14 |
|------|---------|---------|
| 对比基准 | Legacy 代码（Python 重写版 MoM_AllinOne） | Feko 基线 + Mie 解析解 |
| 报告形式 | 内嵌于 Roadmap 表格 | 独立 Markdown + CSV 报告文件 |
| 测试覆盖 | EFIE/CFIE/VEFIE/SCFIE × Direct/MLFMA | 同上 + 新增 PWC/PWCHex/RBF 介质 |
| 精度判据 | vs Legacy RMSE < 1 dB | vs Feko/Mie RMSE < 2 dB (物理意义上的"正确") |

---

## 2. Legacy 对齐基准（数据来源）

### 2.1 Feko 基线文件

路径: `C:\Users\12253\OneDrive\MoM\MoM_AllinOne\deps\compare_feko\`

| 文件名 | 几何体 | 频率 | 适用方程 | 采样规格 |
|--------|--------|------|----------|---------|
| `jet_100MHzRCS.csv` | Jet 战机 (PEC) | 100 MHz | S-EFIE, S-CFIE | θ ∈ [-180°, 180°], 0.5°步长 × 2 个 φ 切面 |
| `sphere_600MHzRCS.csv` | PEC 球 (r≈0.3m) | 600 MHz | S-CFIE, S-MFIE, S-EFIE | 同上 |
| `plate_1dot2GHzRCS.csv` | 介质板 | 1.2 GHz | V-EFIE (SWG), SCFIE | 同上 |
| `plate_metal_1dot2GHzRCS.csv` | 介质+金属板 | 1.2 GHz | SCFIE (RWG+SWG) | 同上 |

**Feko 文件格式**（固定宽度文本）：
```
   THETA    PHI      magn.    phase     magn.    phase        in m*m  ...
 -180.00    0.00   5.222E-01  -36.23  5.946E-04   69.30      3.42739E+00  ...
```
- 列 1: θ (deg)
- 列 2: φ (deg)
- 列 3-4: Eθ 分量 (幅值/相位)
- 列 5-6: Eφ 分量 (幅值/相位)
- 列 7: **RCS (m²)**  ← 主要对比列
- 后续列: 极化率、方向等（不使用）

每文件: 2 行标题 + 1441 行数据（721 点 × 2 个 φ 切面，φ=0° 和 φ=90°）

### 2.2 Mie 解析解（球体）

- EMSuite 已有 `Utilities/MieSeries.jl`（实现了 `mie_rcs_sphere`）
- 用于 PEC 球 600 MHz 的独立精度核查（补充 Feko 对比）
- 球半径从网格文件 `sphere_600MHz.nas` 自动提取

### 2.3 网格文件路径

```
C:\Users\12253\OneDrive\MoM\MoM_AllinOne\meshfiles\
├── jet_100MHz.nas           # Jet PEC 表面三角网格
├── sphere_600MHz.nas        # PEC 球面三角网格
├── plate_1dot2GHz.nas       # 介质板四面体+三角混合网格
└── plate_and_metal_1dot2GHz.nas  # 介质+金属板混合网格
```

---

## 3. 测试矩阵

### 3.1 主测试用例

| ID | 几何体 | 频率 | 方程 | 基函数 | 求解路径 | 基准 | 阈值 |
|----|--------|------|------|--------|----------|------|------|
| **F1** | Jet (PEC) | 100 MHz | S-EFIE | RWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F2** | Jet (PEC) | 100 MHz | S-CFIE (α=0.5) | RWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F3** | Jet (PEC) | 100 MHz | S-EFIE | RWG | MLFMA+GMRES | Feko | RMSE < 3 dB |
| **F4** | Jet (PEC) | 100 MHz | S-CFIE (α=0.5) | RWG | MLFMA+GMRES | Feko | RMSE < 3 dB |
| **F5** | Sphere (PEC) | 600 MHz | S-CFIE (α=0.5) | RWG | Direct (LU) | Feko + Mie | RMSE < 2 dB |
| **F6** | Sphere (PEC) | 600 MHz | S-CFIE (α=0.5) | RWG | MLFMA+GMRES | Feko + Mie | RMSE < 3 dB |
| **F7** | Plate (介质) | 1.2 GHz | V-EFIE | SWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F8** | Plate (介质) | 1.2 GHz | VS-EFIE (SCFIE) | RWG+SWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F9** | Plate+Metal | 1.2 GHz | VS-EFIE (SCFIE) | RWG+SWG | Direct (LU) | Feko | RMSE < 2 dB |

### 3.2 补充测试用例（可选，若时间允许）

| ID | 几何体 | 方程 | 基函数 | 基准 | 说明 |
|----|--------|------|--------|------|------|
| **X1** | Sphere | S-EFIE | RWG | Mie | 单独验证 EFIE 精度 |
| **X2** | Plate (介质) | V-EFIE | PWC | Feko | 与 F7 对比 SWG vs PWC 精度 |
| **X3** | Plate (介质) | VS-EFIE | RWG+PWC | Feko | 与 F8 对比 |

---

## 4. 实施计划

### 4.1 F0 — Feko 数据解析器（独立模块）

**位置**: `benchmark/accuracy/feko_reader.jl`

**接口设计**:
```julia
"""
    read_feko_rcs(filepath) -> (theta_deg, phi_deg, rcs_sqm, rcs_dBsm)

解析 MoM_AllinOne `compare_feko/*.csv` 格式的 Feko 输出文件。
返回:
- theta_deg  :: Vector{Float64}  # θ 角度，单位度
- phi_deg    :: Vector{Float64}  # φ 角度，单位度
- rcs_sqm    :: Vector{Float64}  # RCS，单位 m²
- rcs_dBsm   :: Vector{Float64}  # RCS，单位 dBsm = 10*log10(rcs_sqm)
"""
function read_feko_rcs(filepath::String)
```

**解析逻辑**:
- 跳过标题行（行首含 "THETA" 或空白行）
- 按固定宽度或空格分割，取 col[1]=θ, col[2]=φ, col[7]=RCS(m²)
- 返回完整数组（包含所有 φ 切面的点）

**分离 φ 切面**:
```julia
"""
    split_phi_cuts(theta, phi, rcs_dBsm) -> Dict{Float64, NamedTuple}

按 φ 值分组，返回 φ→(theta, rcs) 的字典。
"""
```

**TDD**: 先在 `test/test_feko_reader.jl` 写测试，验证：
- 文件行数 = 1441（4 个文件一致）
- θ 范围 [-180, 180]，步长 0.5°
- RCS(m²) > 0（物理约束）
- dBsm = 10*log10(sqm) 转换正确

### 4.2 F1 — Mie 解析解参考生成器

**位置**: `benchmark/accuracy/mie_reference.jl`

**接口**:
```julia
"""
    generate_mie_rcs_dBsm(mesh_file, freq_hz, theta_deg_vec) -> (rcs_phi0, rcs_phi90)

从网格文件提取球半径，计算 Mie 解析解 RCS。
对 PEC 球，RCS 不依赖 φ（轴对称），两切面相同。
"""
function generate_mie_rcs_dBsm(mesh_file::String, freq_hz::Real, theta_deg_vec::Vector)
```

**实现**:
- 从 Nastran `.nas` 文件提取节点坐标，计算球半径（节点到原点平均距离）
- 调用现有 `EMSuite.Utilities.mie_rcs_sphere(a, k, theta_vec)`

### 4.3 F2 — 精度指标计算

**位置**: `benchmark/accuracy/accuracy_metrics.jl`

**指标集**:
```julia
struct AccuracyResult
    label     :: String
    n_points  :: Int
    rmse_dB   :: Float64      # √(mean((A-B)²)) in dB
    max_err_dB :: Float64     # max(|A-B|) in dB
    mean_bias_dB :: Float64   # mean(A-B) in dB
    backscatter_err_dB :: Float64  # θ=180° 处的单点误差
    pass      :: Bool         # rmse_dB < threshold
    threshold_dB :: Float64
end
```

```julia
function compute_accuracy(rcs_emsuite_dBsm, rcs_ref_dBsm, label; threshold=2.0)
```

### 4.4 F3 — 各用例仿真脚本

**位置**: `benchmark/accuracy/`

| 脚本 | 覆盖测试 ID |
|------|------------|
| `run_F1_F4_jet.jl` | F1, F2, F3, F4 |
| `run_F5_F6_sphere.jl` | F5, F6 (含 Mie 对比) |
| `run_F7_F9_plate.jl` | F7, F8, F9 |

**每个脚本结构**:
```julia
# 1. 加载网格
# 2. 组装阻抗矩阵 / 构建 MLFMAOperator
# 3. 计算激励向量 (平面波入射，θ_inc=0, φ_inc=0)
# 4. 求解
# 5. 计算 RCS（全球面角度扫描，与 Feko 相同采样点）
# 6. 读取 Feko 参考
# 7. 计算 AccuracyResult
# 8. 输出 CSV + 打印摘要
```

**RCS 计算角度**: 与 Feko 保持一致
- θ ∈ [-180°, 180°]，步长 0.5°（721 点）
- φ 切面: φ=0° 和 φ=90°

### 4.5 F4 — 报告生成器

**位置**: `benchmark/accuracy/generate_report.jl`

**输入**: 各用例仿真脚本产生的 CSV 文件（位于 `test_results/accuracy/`）

**输出**:
- `test_results/accuracy/ACCURACY_REPORT.md` — Markdown 汇总表
- `test_results/accuracy/<ID>_rcs_phi0.csv` — 各用例 φ=0° 切面数据
- `test_results/accuracy/<ID>_rcs_phi90.csv` — 各用例 φ=90° 切面数据
- `test_results/accuracy/<ID>_summary.txt` — 单用例指标摘要

**报告格式**（`ACCURACY_REPORT.md`）:
```markdown
# EMSuite 精度对比报告 (vs Feko / Mie 解析解)

生成时间: 2026-XX-XX  
EMSuite 版本: vX.X.X (commit: xxxxxxx)

## 汇总

| ID | 几何 | 方程 | 求解器 | 基准 | RMSE (dB) | Max Err (dB) | Bias (dB) | 后向散射误差 (dB) | PASS? |
|----|------|------|--------|------|-----------|-------------|-----------|------------------|-------|
| F1 | Jet 100MHz | S-EFIE | Direct | Feko | X.XX | X.XX | ±X.XX | X.XX | ✓/✗ |
...

## 结论

- 通过/失败用例数: X/9
- 最大误差用例: FX (X.XX dB RMSE)
- 所有 Direct 求解路径 RMSE < 2 dB: ✓/✗
- 所有 MLFMA 求解路径 RMSE < 3 dB: ✓/✗
```

---

## 5. 目录结构

```
benchmark/
└── accuracy/
    ├── feko_reader.jl          # F0: Feko CSV 解析器
    ├── mie_reference.jl        # F1: Mie 解析解参考
    ├── accuracy_metrics.jl     # F2: AccuracyResult 结构体 + 指标函数
    ├── run_F1_F4_jet.jl        # F3a: Jet 用例仿真
    ├── run_F5_F6_sphere.jl     # F3b: Sphere 用例仿真
    ├── run_F7_F9_plate.jl      # F3c: Plate 用例仿真
    └── generate_report.jl      # F4: 汇总报告生成

test/
└── test_feko_reader.jl         # F0 单元测试 (TDD)

test_results/
└── accuracy/
    ├── ACCURACY_REPORT.md      # 最终汇总报告
    ├── F1_SEFIE_Jet_Direct_phi0.csv
    ├── F1_SEFIE_Jet_Direct_phi90.csv
    ├── F1_summary.txt
    ├── ...
    └── F9_summary.txt
```

---

## 6. 精度验收门限

| 求解路径 | 几何/方程 | Feko/Mie RMSE 门限 | 最大误差门限 |
|---------|----------|-------------------|------------|
| Direct | S-EFIE, S-CFIE | ≤ 2.0 dB | ≤ 5.0 dB |
| Direct | V-EFIE, SCFIE | ≤ 2.0 dB | ≤ 5.0 dB |
| MLFMA+GMRES | S-EFIE, S-CFIE | ≤ 3.0 dB | ≤ 6.0 dB |
| 任意 | 后向散射 (θ=180°) | 对应门限内 | — |

> **说明**: Feko 本身与 Mie 解析解存在 ~0.5–1 dB 的网格离散误差，因此对比 Feko 的门限适度放宽（相比 Phase 10 的 Legacy 对比）。

---

## 7. 实施顺序与时间估算

| 步骤 | 工作内容 | 预估工时 |
|------|---------|---------|
| **Step 1** | F0: Feko 解析器 + TDD 测试 | 0.5 天 |
| **Step 2** | F1: Mie 参考生成器 | 0.5 天 |
| **Step 3** | F2: AccuracyResult 指标 | 0.5 天 |
| **Step 4** | F3a: run_F1_F4_jet.jl 运行 + 调试 | 1 天 |
| **Step 5** | F3b: run_F5_F6_sphere.jl 运行 + 调试 | 1 天 |
| **Step 6** | F3c: run_F7_F9_plate.jl 运行 + 调试 | 1 天 |
| **Step 7** | F4: 报告生成器 + ACCURACY_REPORT.md | 0.5 天 |
| **Step 8** | 检视迭代 (≥ 2 轮) | 1 天 |
| **合计** | | ~6 天 |

---

## 8. DoD（完成定义）

- [ ] `test/test_feko_reader.jl` 全部测试通过
- [ ] F1–F9 全部用例仿真运行无报错
- [ ] F1–F9 中至少 7/9 通过精度门限
- [ ] `ACCURACY_REPORT.md` 已生成，包含所有 9 个用例的指标表
- [ ] 所有 CSV 数据文件已保存到 `test_results/accuracy/`
- [ ] 检视迭代 ≥ 2 轮 clean（无新问题）
- [ ] `REFACTORING_ROADMAP.md` 中 Phase 14 全部子项勾选
- [ ] `REFACTORING_PROGRESS.md` 更新日志追加 Phase 14 完成记录

---

## 9. 检视迭代计划

Phase 14 完成仿真后，必须进行 ≥ 2 轮检视：

**检视重点**：
1. **算法**: Feko 激励方向是否与 EMSuite 一致（平面波入射角、极化）？
2. **数据对齐**: φ 切面 0° 对应 Legacy 的哪一个切面？角度偏移？
3. **RCS 公式**: EMSuite 的 RCS 定义 (σ = 4π|F|²) 是否与 Feko 一致（面积归一化/距离归一化）？
4. **报告完整性**: 是否覆盖所有基函数（SWG/RWG/PWC/RBF）？
5. **软件工程**: 解析器是否有充分的错误处理？文件路径是否可配置？

---

## 10. 待机链接

- 路线图更新位置: [REFACTORING_ROADMAP.md](REFACTORING_ROADMAP.md) → Phase 14 节
- 进度追踪位置: [REFACTORING_PROGRESS.md](REFACTORING_PROGRESS.md) → Phase 14 条目
- Feko 基线路径: `C:\Users\12253\OneDrive\MoM\MoM_AllinOne\deps\compare_feko\`
- 网格文件路径: `C:\Users\12253\OneDrive\MoM\MoM_AllinOne\meshfiles\`
