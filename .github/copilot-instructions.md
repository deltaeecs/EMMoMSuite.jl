# EMSuite.jl Copilot Instructions

## 计划与进度文件

- **重构路线图**: [REFACTORING_ROADMAP.md](.github/REFACTORING_ROADMAP.md)
- **重构进度**: [REFACTORING_PROGRESS.md](.github/REFACTORING_PROGRESS.md)

> **原则**: 每次完成具有实质意义的进展后，必须同步更新上述两个文件。

---

## 核心开发原则

### 1. TDD 工作流 (必须遵循)
1. **RED**: 先写失败测试，定义期望接口和行为
2. **GREEN**: 写最简代码使测试通过
3. **REFACTOR**: 改善结构，保持测试绿色

### 2. 严格 Legacy 对齐
- **唯一真相源**: Legacy 代码 (`MoM_Kernels`, `MoM_Basics`, `MoM_AllinOne`)
- **禁止**: 使用 Mie 级数或其他解析解作为调试基准
- **禁止**: 添加经验常数 (`1/23π`, `22.0/k²` 等) 来校准结果
- 结果不匹配时，算法有误，必须在 Legacy 代码中找到根源

### 3. 差异排查流程
当结果与 Legacy 不一致时，逐项比较:
1. **几何**: 顶点坐标、边长、基函数定义 (符号、支撑)
2. **常数**: $k$, $\eta$, $1/4\pi$ vs $1/16\pi$
3. **积分**:
   - 矢量势项: $\int \mathbf{f} \cdot \mathbf{f}' G$
   - 标量势项: $\int (\nabla \cdot \mathbf{f}) (\nabla \cdot \mathbf{f}') G$
   - 奇异项: $F_1$ (1/R), $F_2$ (Rho·Rho/R)
4. **组装**: 矩阵元素幅度和相位

### 4. 进度同步原则
- 每次有实质进展 (功能实现、验证通过、Bug 修复) 后，**必须**更新:
  - `REFACTORING_ROADMAP.md` 中对应任务的勾选状态
  - `REFACTORING_PROGRESS.md` 中的详细状态和更新日志
- 不在本提示词文件中直接记录进度

### 5. Git 提交规范 (必须遵循)
- **每次对源码的修改确认成功后，必须立即提交代码**
- 提交粒度: 一个逻辑完整的修改 = 一次提交 (Bug 修复、功能新增、测试添加、文档更新等)
- 提交信息格式: `<type>: <简要描述>`
  - `fix:` Bug 修复
  - `feat:` 新功能
  - `test:` 测试添加/修改
  - `docs:` 文档更新
  - `refactor:` 重构
  - `bench:` 基准测试脚本
  - `chore:` 杂项 (配置、脚本等)
- 禁止积累大量未提交的修改；禁止将不相关的修改混入同一次提交

### 6. 其他约定
- **跳过 Lebedev 测试**: 除非关键问题，假设 Lebedev 路径正确
- **量级检查**: 用数量级检查因子 ($4\pi$, $k^2$) 而非完整解析推导
- **代码规范**: PascalCase 类型, snake_case 函数, UPPER_SNAKE_CASE 常量
