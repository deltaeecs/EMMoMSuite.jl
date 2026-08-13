# EMMoMSuite Test Configuration

## MPI/混合并行测试

以下测试需要 MPI 启动器（MSMPI），不能并入常规套件：

```powershell
mpiexec -n 2 julia -t 2 --project=. test/test_hybrid_mlfma.jl   # MLFMA MPI×线程精度门
mpiexec -n 2 julia -t 2 --project=. test/test_hybrid_aim.jl     # AIM  MPI×线程精度门
mpiexec -n 2 julia -t 2 --project=. test/test_hybrid_pmchw.jl   # PMCHW MPI×线程精度门
mpiexec -n 4 julia -t 2 --project=. test/test_distributed_lu.jl # 分布式稠密 LU 精度门
mpiexec -n 2 julia -t 2 --project=. test/test_scalapack_lu.jl   # ScaLAPACK（本机 MinGW/MSMPI）分布式稠密 LU 精度门
mpiexec -n 4 julia -t 2 --project=. test/test_hybrid_preconditioner.jl # MPI 分布式预条件门（BlockJacobi/Diagonal）
```

ScaLAPACK 精度门需要本机 ScaLAPACK 动态库。库路径自动探测（常见 MSYS2
mingw64/ucrt64/clang64 路径与 PATH），也可用环境变量 `SCALAPACK_LIB_PATH`
显式指定（最高优先）；未找到时 `scalapack_lu_solve` 会给出安装指引而不是崩溃。
安装示例（MSYS2，MSMPI 版）：
```powershell
pacman -S mingw-w64-x86_64-scalapack      # mingw64（本机验证 2.2.2）
# 或 mingw-w64-ucrt-x86_64-scalapack      # ucrt64
```
无 ScaLAPACK 环境时可改用自研 MPI LU（`mpi_lu!`/`mpi_lu_solve!`，见
`test/test_distributed_lu.jl`）。

混合效率基准（秩内 FFTW 线程在脚本内设置）：

```powershell
mpiexec -n 2 julia -t 4 --project=. benchmark/benchmark_hybrid_mlfma.jl 2.0 24 48 40
mpiexec -n 2 julia -t 4 --project=. benchmark/benchmark_hybrid_aim.jl   2.0 24 48 40
mpiexec -n 4 julia -t 2 --project=. benchmark/benchmark_scalapack_lu.jl 9600   # 分布式稠密 LU（ScaLAPACK）效率基准
```

## Test Tiers

EMMoMSuite 测试套件分为三个层次，以平衡速度和覆盖率：

### 🚀 Fast Tests (`runtests_fast.jl`)
**运行时间**: < 5 分钟  
**用途**: 本地开发快速验证，CI/CD 主流程

**包含**:
- 核心单元测试（Materials, Geometry, BasisFunctions）
- 积分方程基础测试（IntegralEquations, FastExp）
- 覆盖缺口测试（Phase 21 回归测试）
- 工具函数测试（IO, Lebedev, Ports）

**运行方式**:
```bash
julia --project=. test/runtests_fast.jl
```

### ⚙️ Medium Tests (`runtests_medium.jl`)
**运行时间**: < 30 分钟  
**用途**: PR 验证，夜间构建

**包含**:
- 求解器测试（Solvers, Preconditioners）
- MLFMA 测试
- 并行测试（Threads, MPI）
- 中等规模数值验证
- 体积方法测试（SCFIE, PWC）
- 后处理测试

**运行方式**:
```bash
julia --project=. --threads=4 test/runtests_medium.jl
```

### 🔬 Full Tests (`runtests_full.jl`)
**运行时间**: < 2 小时  
**用途**: 发布前验证，完整覆盖率分析

**包含**:
- Fast + Medium 所有测试
- 大规模集成测试
- PMCHW/N-Müller 完整验证
- Legacy 对齐测试
- 性能基准测试
- 特殊场景测试

**运行方式**:
```bash
julia --project=. --threads=8 test/runtests_full.jl
```

## 运行特定测试

```bash
# 单个测试文件
julia --project=. test/test_coverage_gaps.jl

# 使用 Test 模块
julia --project=. -e 'using Test; include("test/test_coverage_gaps.jl")'

# 带覆盖率
julia --project=. --code-coverage=user test/runtests_fast.jl
```

## CI/CD 集成

GitHub Actions 自动运行：
- **每次 push/PR**: Fast Tests（所有平台）
- **PR 通过后**: Medium Tests（Ubuntu）
- **定期/手动**: Full Tests + 覆盖率报告

## 本地开发建议

**提交前**:
```bash
julia --project=. test/runtests_fast.jl
```

**重要修改后**:
```bash
julia --project=. --threads=4 test/runtests_medium.jl
```

**发布前**:
```bash
julia --project=. --threads=8 test/runtests_full.jl
```

## 覆盖率目标

| 模块 | 目标覆盖率 | 当前状态 |
|------|-----------|---------|
| Core | ≥90% | ~92% |
| IntegralEquations | ≥85% | ~88% |
| FastAlgorithms | ≥80% | ~82% |
| Solvers | ≥85% | ~87% |
| Utilities | ≥90% | ~91% |
| **Overall** | **≥85%** | **~87%** |

## 添加新测试

1. 单元测试 → 添加到相应的 `test_*.jl` 或创建新文件
2. 快速测试 → 添加到 `runtests_fast.jl`
3. 中等测试 → 添加到 `runtests_medium.jl`
4. 慢速测试 → 添加到 `runtests_full.jl`

## Phase 21 测试增强

- [x] Phase 21.1 Round 1: 覆盖缺口测试（`test_coverage_gaps.jl`）
- [x] Phase 21.5: 测试分层结构（fast/medium/full）
- [x] Phase 21.5: CI/CD 配置更新
- [ ] Phase 21.2: 集成测试增强（vs Legacy）
- [ ] Phase 21.3: 数值精度验证
- [ ] Phase 21.4: 并行测试扩展
