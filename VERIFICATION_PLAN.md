# EMSuite Verification Plan (Updated Jan 2026)

This plan outlines the remaining verification steps for `EMSuite`.
**CRITICAL INSTRUCTION**: The verification focus is **STRICT LEGACY PARITY**.
**DO NOT** compare against Mie Series or other theoretical benchmarks for debugging. The Mie series implementation itself may be incorrect.
The **ONLY** source of truth is the Legacy Code (`MoM_Kernels`, `MoM_Basics`, `MoM_AllinOne`).

## 1. Completed Verification
- [x] **Surface Integral Equations (PEC Sphere)**
    - **EFIE**: Verified (Internal Resonance observed as expected).
    - **CFIE**: Verified (Resonance suppressed).
    - **Direct Solver**: Verified.
    - **MLFMA**: Verified (Forward/Backscatter match Theory).
    - **Matrix Accuracy**: Confirmed `EMSuite` uses higher-precision analytical integrals than Legacy.

## 2. Upcoming Verification Tasks

### 2.1 Volume Integral Equation (VEFIE)
**Goal**: Verify scattering from a dielectric object (e.g., Dielectric Sphere).
- **Target**: Dielectric Sphere ($r=0.5m, \epsilon_r=2.0, f=300MHz$).
- **Metrics**:
    - **RCS Accuracy**: **Compare against Legacy Code (`MoM_AllinOne`)**.
    - **Planes**: E-plane ($\phi=0^\circ$) AND H-plane ($\phi=90^\circ$).
    - **Tolerance**: Error $< 1.0 \text{ dB}$ vs Legacy.
- **Steps**:
    1.  **Direct Solver**: Verify `verify_VEFIE_direct.jl` against Legacy output.
    2.  **MLFMA**: Verify `verify_VEFIE_mlfma.jl` against Direct Solver.
    3.  **Debug Flow**: If RCS fails -> Check Matrix Elements ($Z_{vv}$) vs Legacy -> Check RHS ($V$) vs Legacy.

### 2.2 Surface-Volume Integral Equation (SVIE/SCFIE)
**Goal**: Verify scattering from a composite object (PEC + Dielectric).
- **Target**: Coated Sphere or Composite Geometry.
- **Metrics**:
    - **RCS Accuracy**: **Compare against Legacy Code**.
    - **Planes**: E-plane AND H-plane.
    - **Tolerance**: Error $< 1.0 \text{ dB}$ vs Legacy.
- **Steps**:
    1.  **Direct Solver**: Verify `verify_SCFIE_direct.jl`.
    2.  **MLFMA**: Verify `verify_SCFIE_mlfma.jl`.

### 2.3 Parallel Verification (MPI)
**Goal**: Ensure MLFMA works correctly across multiple processes.
- **Target**: PEC Sphere (Large Mesh).
- **Steps**:
    1.  Run `benchmark_parallel_sphere.jl` with `mpiexec -n 4`.
    2.  Compare RCS with Serial result.
    3.  Tolerance: Exact match (Machine Precision).

## 3. Debugging Protocol
If RCS verification fails ($> 1 \text{ dB}$ error vs Legacy):
1.  **Check Geometry**: Mesh quality, normals, material properties.
2.  **Check Matrix ($Z$)**:
    - Extract small sub-matrix (e.g., 2 elements).
    - Compare **Legacy vs EMSuite** element-by-element.
3.  **Check Excitation ($V$)**:
    - Verify incident field phase and polarization against Legacy.
4.  **Check Solver**:
    - Verify convergence history (Residual $< 10^{-3}$).

## 4. Success Criteria
- **RCS**: E-plane and H-plane match **Legacy** within 1 dB.
- **Stability**: No internal resonances for CFIE/PMCHWT.
- **Performance**: MLFMA scales as $O(N \log N)$.
- [ ] **验证**: 确保结果与串行版本一致（无竞争冒险）。

---

## 执行计划 (Action Plan)

1.  **停止 Mie 对比**: 废弃 `verify_VEFIE_mie.jl`，不再使用理论解作为判断依据。
2.  **创建 Legacy 对比脚本**: 编写 `scripts/verification/verify_VEFIE_legacy.jl`。
    - 运行 `MoM_AllinOne` 的 VEFIE 算例，保存 RCS 结果。
    - 运行 `EMSuite` 的 VEFIE 算例，保存 RCS 结果。
    - 对比两者差异。
3.  **矩阵元素对比**: 如果 RCS 差异大，编写脚本对比 Z 矩阵元素。
4.  **全流程测试**: 确保所有算子 (EFIE/MFIE/CFIE/VEFIE) 与 Legacy 保持一致。
