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

## 1b. 新增范围（本次更新）

除原 RCS 散射用例外，新增：

| 类别 | 测试 | 参考基准 | MLFMA | 说明 |
|------|------|---------|-------|------|
| **PMCHW** | 介质球散射 RCS | Mie 介质级数（已实现） | 需新实现 | P1–P2 |
| **天线（集总端口）** | 半波偶极子辐射方向图、输入阻抗 | 解析公式 | 需新实现 | A1–A2 |
| **天线（差分/同轴端口）** | 完整 S 参数 | 解析公式 | 可选 | A3–A4 |

> **MLFMA PMCHW 与天线 MLFMA** 是新功能实现任务，在 Phase 14 测试框架搭建后作为子任务完成（见第 4.6 节）。

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

`Utilities/MieSeries.jl` 已实现两类：

| 函数 | 用途 |
|------|------|
| `calculate_mie_rcs_pec_sphere(radius, freq, theta)` | PEC 球 RCS → 用于 F5/F6/X1 |
| `calculate_mie_rcs_dielectric_sphere(radius, freq, theta, eps_r, mu_r)` | 介质球 RCS → 用于 P1/P2（PMCHW 验证） |

球半径从网格文件坐标节点到原点距离计算。

### 2.3 天线解析参考基准

天线问题**不存在 Feko 基线**，使用经典理论解析值作为精度门限：

| 天线类型 | 频率 | 解析参考 | 精度指标 |
|---------|------|---------|--------|
| 半波偶极子（0.5λ, PEC, 自由空间） | 300 MHz (λ=1m, L=0.5m) | 输入阻抗 Z_in = 73.1 + j42.5 Ω | `\|Z_in_EMSuite - Z_in_analytic\|/\|Z_in_analytic\|` < 5% |
| 半波偶极子 | 同上 | 最大方向性 D_max = 1.64 (2.15 dBi) | < 1 dBi 误差 |
| 半波偶极子 | 同上 | 辐射方向图 E 面 RMSE vs 解析 sinc | < 1 dB |
| 谐振偶极子（0.47λ，无虚部阻抗） | 318 MHz (L=0.47m) | Im(Z_in) ≈ 0，Re(Z_in) ≈ 73 Ω | Im(Z_in) < 5 Ω |

> 注：半波偶极子 Z_in 的解析值来自 Balanis《Antenna Theory》表 8.1。

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

### 3.1 散射问题（RCS 精度，vs Feko/Mie PEC）

| ID | 几何体 | 频率 | 方程 | 基函数 | 求解路径 | 基准 | 阈值 |
|----|--------|------|------|--------|----------|------|------|
| **F1** | Jet (PEC) | 100 MHz | S-EFIE | RWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F2** | Jet (PEC) | 100 MHz | S-CFIE (α=0.5) | RWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F3** | Jet (PEC) | 100 MHz | S-EFIE | RWG | MLFMA+GMRES | Feko | RMSE < 3 dB |
| **F4** | Jet (PEC) | 100 MHz | S-CFIE (α=0.5) | RWG | MLFMA+GMRES | Feko | RMSE < 3 dB |
| **F5** | Sphere (PEC) | 600 MHz | S-CFIE (α=0.5) | RWG | Direct (LU) | Feko + Mie PEC | RMSE < 2 dB |
| **F6** | Sphere (PEC) | 600 MHz | S-CFIE (α=0.5) | RWG | MLFMA+GMRES | Feko + Mie PEC | RMSE < 3 dB |
| **F7** | Plate (介质) | 1.2 GHz | V-EFIE | SWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F8** | Plate (介质) | 1.2 GHz | VS-EFIE (SCFIE) | RWG+SWG | Direct (LU) | Feko | RMSE < 2 dB |
| **F9** | Plate+Metal | 1.2 GHz | VS-EFIE (SCFIE) | RWG+SWG | Direct (LU) | Feko | RMSE < 2 dB |

### 3.2 PMCHW 均匀介质体散射（vs Mie 介质级数）

| ID | 几何体 | 频率 | 方程 | 基函数 | 求解路径 | 基准 | 阈值 |
|----|--------|------|------|--------|----------|------|------|
| **P1** | 介质球 (εᵣ=4, 无损, r=0.15m) | 600 MHz | PMCHW | RWG (2N DOF) | Direct (LU) | Mie 介质级数 | RMSE < 2 dB |
| **P2** | 介质球 (εᵣ=4, 无损, r=0.15m) | 600 MHz | PMCHW | RWG (2N DOF) | MLFMA+GMRES | Mie 介质级数 | RMSE < 3 dB |
| **P3** | 介质球 (εᵣ=2.2-j0.1，有损) | 300 MHz | PMCHW | RWG | Direct (LU) | Mie 介质级数 | RMSE < 2 dB |

> **P2 前置条件**: 需先实现 PMCHW 的 MLFMA 算子（2×2 块结构，见 4.6 节）。

### 3.3 天线端口辐射（vs 解析公式）

| ID | 几何体 | 频率 | 方程 | 端口类型 | 求解路径 | 基准 | 阈值 |
|----|--------|------|------|---------|----------|------|------|
| **A1** | 半波偶极子 (L=0.5m, PEC) | 300 MHz | S-EFIE | LumpedPort（delta-gap，单边） | Direct (LU) | 解析 Z_in/方向图 | \|ΔZ_in\|/\|Z_in\| < 5%，D_max 误差 < 1 dBi |
| **A2** | 半波偶极子 (L=0.5m, PEC) | 300 MHz | S-EFIE | LumpedPort | MLFMA+GMRES | vs A1 结果 | RMSE < 0.5 dB |
| **A3** | 谐振偶极子 (L=0.47m，无虚部) | 319 MHz | S-EFIE | LumpedPort | Direct (LU) | Im(Z_in) ≈ 0 | \|Im(Z_in)\| < 5 Ω |
| **A4** | 半波偶极子 + 50Ω 匹配 | 300 MHz | S-EFIE | LumpedPort (负载端口) | Direct (LU) | S11 解析值 | \|S11_dB\| < 0.5 dB误差 |

> **网格**: 偶极子需从 Gmsh 或 Nastran 生成细长线段三角面（模拟导线表面），或使用已有线天线网格。

### 3.4 补充测试（可选）

| ID | 几何体 | 方程 | 基准 | 说明 |
|----|--------|------|------|------|
| **X1** | Sphere (PEC) | S-EFIE Direct | Mie PEC | 单独验证 EFIE |
| **X2** | Plate (介质) | V-EFIE, PWC | Feko | SWG vs PWC 精度对比 |
| **X3** | 介质球 | PMCHW, εᵣ=10 | Mie 介质 | 高对比度介质测试 |

---

## 4. 实施计划

### 4.1 F0 — Feko 数据解析器（独立模块）

**位置**: `benchmark/accuracy/feko_reader.jl`

**接口设计**:
```julia
"""
    read_feko_rcs(filepath) -> (theta_deg, phi_deg, rcs_sqm, rcs_dBsm)

解析 MoM_AllinOne `compare_feko/*.csv` 格式的 Feko 输出文件。
"""
function read_feko_rcs(filepath::String)
```

**解析逻辑**:
- 跳过标题行（行首含 "THETA" 或空白行）
- 按固定宽度或空格分割，取 col[1]=θ, col[2]=φ, col[7]=RCS(m²)
- 返回完整数组（包含所有 φ 切面的点）

**TDD**: 先写 `test/test_feko_reader.jl`，验证：
- 文件行数 = 1441（4 个文件一致）
- θ 范围 [-180, 180]，步长 0.5°
- RCS(m²) > 0
- dBsm = 10*log10(sqm) 转换正确

### 4.2 F1 — 参考基准生成器

**位置**: `benchmark/accuracy/reference_data.jl`

```julia
# Mie PEC 球
function mie_pec_rcs_dBsm(mesh_file, freq_hz, theta_deg_vec)

# Mie 介质球（用于 PMCHW）
function mie_dielectric_rcs_dBsm(radius_m, freq_hz, eps_r, mu_r, theta_deg_vec)

# 半波偶极子解析输入阻抗
function dipole_halfwave_Zin_analytic()  # → 73.1 + j42.5 Ω (精确半波)
function dipole_resonant_Zin_analytic()  # → ~73.0 + j0 Ω (0.47λ)

# 半波偶极子解析方向图（E面，sin(θ)加权）
function dipole_halfwave_farfield_analytic(theta_vec)  # → E_theta(θ)，归一化
```

### 4.3 F2 — 精度指标计算

**位置**: `benchmark/accuracy/accuracy_metrics.jl`

```julia
struct AccuracyResult
    label          :: String
    n_points       :: Int
    rmse_dB        :: Float64
    max_err_dB     :: Float64
    mean_bias_dB   :: Float64
    backscatter_err_dB :: Float64
    pass           :: Bool
    threshold_dB   :: Float64
end

struct AntennaAccuracyResult
    label              :: String
    Zin_rel_err        :: Float64    # |ΔZ|/|Z_ref|
    D_max_err_dBi      :: Float64    # |D_max_EMSuite - D_max_analytic| in dBi
    pattern_rmse_dB    :: Float64    # 方向图 RMSE (dB)
    S11_err_dB         :: Float64    # |S11_EMSuite - S11_analytic| in dB
    Im_Zin             :: Float64    # Im(Z_in)，用于谐振测试
    pass               :: Bool
end
```

### 4.4 F3 — 仿真脚本

**位置**: `benchmark/accuracy/`

| 脚本 | 覆盖 ID |
|------|--------|
| `run_F1_F4_jet.jl` | F1–F4（Jet S-EFIE/S-CFIE） |
| `run_F5_F6_sphere.jl` | F5–F6（Sphere S-CFIE + Mie PEC） |
| `run_F7_F9_plate.jl` | F7–F9（Plate V-EFIE/SCFIE） |
| `run_P1_P3_pmchw.jl` | P1–P3（PMCHW 介质球 + Mie 介质） |
| `run_A1_A4_antenna.jl` | A1–A4（半波偶极子 + LumpedPort） |

**天线仿真流程**（`run_A1_A4_antenna.jl`）：
```julia
# 1. 生成偶极子网格（Gmsh API 或程序生成）：
#    L = 0.5m，直径 d = 0.01m，沿 z 轴，表面三角剖分
#    f = 300 MHz → λ = 1m，单元尺寸 ≈ λ/20 = 0.05m
# 2. 找到中央边作为 delta-gap 端口边
# 3. 构造 LumpedPort(id=1, edge_idx=central_edge, impedance=50Ω, type=:voltage_source, voltage=1V)
# 4. 装配 EFIE 阻抗矩阵
# 5. 添加集总端口阻抗贡献: assemble_lumped_port_impedance!(Z, port)
# 6. 构造激励向量: add_port_excitation!(V, port)
# 7. 求解: I = Z \ V
# 8. 计算输入阻抗: Z_in = port_voltage(I, port) / port_current(I, port)
# 9. 计算方向图: antenna_directivity(θs, ϕs, I, basis)
# 10. 与解析值对比
```

**PMCHW 仿真流程**（`run_P1_P3_pmchw.jl`）：
```julia
# 1. 加载球面网格（sphere_600MHz.nas，r=0.15m，重用或生成新网格）
# 2. 构造 PMCHW(freq, eps_r=4.0, mu_r=1.0)
# 3. 装配 2N×2N 阻抗矩阵
# 4. 构造平面波激励（Excitation 模块的 PMCHW 变体）
# 5. 求解 2N×2N 系统
# 6. 后处理：从 J（前 N）+ M（后 N）计算散射场/RCS
# 7. 与 Mie 介质球对比
```

### 4.5 F4 — 报告生成器

**位置**: `benchmark/accuracy/generate_report.jl`

**输出目录**: `test_results/accuracy/`

**报告结构**（`ACCURACY_REPORT.md`）：
```
# EMSuite 精度对比报告 (vs Feko / Mie / 解析解)
生成时间 / 版本 / commit

## Part I: RCS 散射精度 (F1–F9, vs Feko)
| ID | 几何 | 方程 | 求解器 | RMSE | Max Err | Bias | 后向散射误差 | PASS? |

## Part II: PMCHW 介质体散射 (P1–P3, vs Mie 介质级数)
| ID | 几何 | ε_r | 求解器 | RMSE | Max Err | PASS? |

## Part III: 天线辐射（A1–A4, vs 解析）
| ID | 几何 | 端口 | 求解器 | ΔZ_in (%) | D_max 误差 (dBi) | Pattern RMSE | PASS? |

## 结论
通过/失败统计, 最大误差用例
```

### 4.6 MLFMA PMCHW 实现（新功能，Phase 14 核心任务）

**当前状态**: `PMCHW.jl` 仅支持 Direct (LU)。  
**目标**: 为 PMCHW 的 2N×2N 系统构造 `MLFMAOperator` 或 2×2 块 MLFMA 结构。

**实现方案**：

```
PMCHW MLFMA 算子结构：
              
  [EFIE_far(k0) + EFIE_far(k1)        K_far(k0) + K_far(k1)      ]   [J]
  [-(K_far(k0) + K_far(k1))           EFIE_far_inv(k0) + EFIE_far_inv(k1)]   [M]
```

**Block MLFMA 设计**：
- 每个媒质 (k0, k1) 独立构建 MLFMA 树（共用几何，但波数不同）
- `mul!(y, op::PMCHWMLFMAOperator, x)` 接受长度 2N 向量，拆分为 J/M 两部分，分别调用各子算子
- 近场 Z_near 仍用 Dense 直接计算（2N×2N 子块）
- 预条件器：2×2 Block-Jacobi，子块分别为 Z_EJ_near 和 Z_HM_near

**实施步骤**：
1. 定义 `PMCHWMLFMAOperator` 结构体（包含 2 个 EFIE MLFMAOperator + 2 个 K 远场算子）
2. 实现 `build_pmchw_mlfma_operator(pmchw, basis, ...)`
3. 实现 `LinearAlgebra.mul!(y, op, x)` — 4 个远场 mul 的线性组合
4. 加入 GMRES 迭代求解路径
5. 实现 P2 测试用例验证

**天线 MLFMA**（A2）：
- 标准 S-EFIE MLFMAOperator（已有），直接复用
- LumpedPort 作为激励向量修正，不影响 MLFMA 算子本身
- 仅需在 GMRES 迭代路径中正确处理 Z_port 的近场贡献

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
| **Step 2** | F1: 参考基准生成器（Mie PEC/介质 + 偶极子解析） | 0.5 天 |
| **Step 3** | F2: `AccuracyResult` + `AntennaAccuracyResult` | 0.5 天 |
| **Step 4** | F3a-c: Jet/Sphere/Plate 用例运行 + 调试 | 2 天 |
| **Step 5** | F3d: PMCHW Direct 用例 + Mie 介质对比 (P1, P3) | 1 天 |
| **Step 6** | `PMCHWMLFMAOperator` 实现 + 单元测试 | 2 天 |
| **Step 7** | P2: PMCHW MLFMA 验证 | 0.5 天 |
| **Step 8** | F3e: 偶极子天线 + LumpedPort 用例 (A1–A4) | 1.5 天 |
| **Step 9** | F4: 报告生成器 + `ACCURACY_REPORT.md` | 0.5 天 |
| **Step 10** | 检视迭代 (≥ 2 轮) | 1 天 |
| **合计** | | **~10 天** |

---

## 8. DoD（完成定义）

**基础设施**:
- [ ] `test/test_feko_reader.jl` 全部测试通过
- [ ] `benchmark/accuracy/` 目录结构完整（包括 reference_data.jl）

**F1–F9 散射精度**:
- [ ] F1–F9 全部用例仿真运行无报错
- [ ] F1–F9 至少 7/9 通过精度门限

**PMCHW 介质精度**:
- [ ] P1、P3 Direct 仿真完成，与 Mie 介质 RMSE < 2.5 dB
- [ ] `PMCHWMLFMAOperator` 实现，单元测试通过
- [ ] P2 MLFMA 验证，与 P1 Direct RMSE < 1 dB（一致性）

**天线端口**:
- [ ] 偶极子网格已生成 (Gmsh 或程序自动)
- [ ] A1 输入阻抗 \|ΔZ_in\|/\|Z_analytic\| < 5%
- [ ] A1 最大方向性 D_max 误差 < 1 dBi
- [ ] A3 调振 Im(Z_in) < 5 Ω

**报告**:
- [ ] `ACCURACY_REPORT.md` 已生成，包含 F/P/A 三类用例指标表
- [ ] 所有 CSV 数据文件保存到 `test_results/accuracy/`

**质量门控**:
- [ ] 检视迭代 ≥ 2 轮 clean（无新问题）
- [ ] `REFACTORING_ROADMAP.md` 中 Phase 14 全部子项勾选
- [ ] `REFACTORING_PROGRESS.md` 更新日志追加 Phase 14 完成记录

---

## 9. 检视迭代计划

Phase 14 完成仿真后，必须进行 ≥ 2 轮检视：

**检视重点**：
1. **激励对齐**: Feko 平面波入射方向是否与 EMSuite 一致（入射角、极化）？
2. **数据对齐**: φ 切面 0° 对应 Feko 哪一个切面？角度偏移？
3. **RCS 公式**: EMSuite 的 RCS 定义 (σ = 4π|F|²) 是否与 Feko 一致？
4. **PMCHW 验证**: J/M 向量分割正确？散射场 RCS 后处理公式是否包含 J 和 M 两者贡献？
5. **天线验证**: delta-gap 端口边索引是否正确？端口电压/电流提取方式？
6. **PMCHW MLFMA**: 2x2 块结构的近场节点占位是否正确？预条件器是否覆盖 J/M 全部分量？
7. **软件工程**: 解析器是否有充分错误处理？文件路径是否可配置？

---

## 10. 待机链接

- 路线图更新位置: [REFACTORING_ROADMAP.md](REFACTORING_ROADMAP.md) → Phase 14 节
- 进度追踪位置: [REFACTORING_PROGRESS.md](REFACTORING_PROGRESS.md) → Phase 14 条目
- Feko 基线路径: `C:\Users\12253\OneDrive\MoM\MoM_AllinOne\deps\compare_feko\`
- 网格文件路径: `C:\Users\12253\OneDrive\MoM\MoM_AllinOne\meshfiles\`
