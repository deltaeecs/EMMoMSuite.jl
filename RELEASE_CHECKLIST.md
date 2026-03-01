# EMSuite.jl 发布清单 (Julia General Registry)

## 前置条件

在发布到 Julia General Registry 之前，需要完成以下步骤：

---

## 1. 补全 Project.toml 元数据

打开 `Project.toml`，填写实际值：

```toml
name = "EMSuite"
uuid = "<运行 `import UUIDs; UUIDs.uuid4()` 生成唯一 UUID>"
version = "0.1.0"
authors = ["Your Real Name <your.email@example.com>"]
```

生成 UUID 的方法：
```julia
julia -e "import UUIDs; println(UUIDs.uuid4())"
```

---

## 2. 本地依赖处理

当前 `EMSuite` 依赖以下本地（未发布）包：
- `MoM_AllinOne`
- `MoM_Basics`
- `MoM_Kernels`
- `MoM_Visualizing`

选项：
- **A. 先发布子包**: 将这些包独立发布到 Registry，然后在 `[deps]` 中保留标准引用。
- **B. 内联代码**: 将这些包的代码移入 `EMSuite/src/` 并删除外部依赖。
- **C. 私有 Registry**: 使用组织内部的 LocalRegistry.jl。

---

## 3. GitHub Repository 设置

1. 在 GitHub 创建与包名一致的仓库：`EMSuite.jl`
2. 推送代码：
   ```powershell
   git remote add origin https://github.com/<your-org>/EMSuite.jl.git
   git push -u origin master
   ```
3. 确保仓库是公开的（General Registry 要求）。

---

## 4. 打版本标签

```powershell
git tag v0.1.0
git push origin v0.1.0
```

---

## 5. 注册包

方法 A — GitHub App（推荐）：
1. 安装 [JuliaRegistrator](https://github.com/JuliaRegistries/Registrator.jl#registering-a-package-in-the-general-registry) GitHub App。
2. 在 GitHub 上对 `v0.1.0` 打标签的 commit 发评论：
   ```
   @JuliaRegistrator register
   ```
3. 等待自动 PR 到 General Registry。

方法 B — 本地（LocalRegistry.jl）：
```julia
using LocalRegistry
create_registry("EMSuiteRegistry", "git@github.com:<your-org>/EMSuiteRegistry.git")
register("EMSuite")
```

---

## 6. CI/CD 建议 (`.github/workflows/`)

创建 `CI.yml`：
```yaml
name: CI
on:
  push:
    branches: [master]
  pull_request:
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
        with:
          version: '1.10'
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
```

---

## 7. 发布后

- [ ] 在 `CHANGELOG.md` 中将 `[Unreleased]` 替换为 `[0.1.0] — YYYY-MM-DD`
- [ ] 更新 `docs/make.jl` 中的 `deploydocs` 配置（取消注释并填写仓库 URL）
- [ ] 开始开发 `v0.2.0`（功能分支）

---

## 当前状态检查

| 项目 | 状态 |
|------|------|
| uuid (真实值) | ⚠️ 需替换占位符 |
| authors (真实值) | ⚠️ 需替换占位符 |
| 本地依赖 (MoM_*) | ⚠️ 需要处理 |
| GitHub 仓库 | ⚠️ 需要创建 |
| 测试套件 (179/179) | ✅ |
| JuliaFormatter 格式化 | ✅ |
| CHANGELOG.md | ✅ |
| compat 语法 | ✅ |
| docs 构建 | ✅ |
