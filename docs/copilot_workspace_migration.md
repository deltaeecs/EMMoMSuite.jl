# Copilot Workspace Migration

本文件记录一个尽量保持 GitHub Copilot 会话连续性的工作区迁移方案。

## 目标

当工程从 OneDrive 迁移到本地磁盘时，尽量不要继续依赖“直接打开工程文件夹”的方式，而是改为始终打开固定路径的 workspace 文件。

当前固定 workspace 文件位置：

- `C:\Users\12253\AppData\Local\VSCodeWorkspaces\EMMoMSuite.code-workspace`

## 当前状态

- 固定 workspace 文件已经创建。
- 它当前指向 `C:\Users\12253\OneDrive\MoM\EMMoMSuite`。
- 后续真正迁移到本地后，只需要修改这个 workspace 文件中的 `folders[0].path`。

## 建议迁移步骤

1. 先使用 `C:\Users\12253\AppData\Local\VSCodeWorkspaces\EMMoMSuite.code-workspace` 打开工程，而不是直接打开 OneDrive 下的文件夹。
2. 将工程完整复制到本地非 OneDrive 路径，例如 `D:\Dev\EMMoMSuite`。
3. 编辑 workspace 文件，把 `folders[0].path` 从 OneDrive 路径改成新的本地路径。
4. 之后固定通过这个 workspace 文件打开工程。

## 原因

- VS Code 文档明确区分“打开文件夹”和“打开 `.code-workspace` 文件”两种 workspace 形态。
- 当前 Copilot 会话数据位于 `workspaceStorage` 下，实际行为更接近“与 workspace 身份绑定”而不是“与仓库逻辑名称绑定”。
- 官方文档没有保证“改文件夹路径后自动继承同一会话历史”，因此固定 workspace 文件路径是更稳妥的工程化做法。

## 限制

- 该方案是为了“尽量保持”历史连续，不是 GitHub Copilot 官方承诺的强保证。
- 如果未来需要再次搬迁工程目录，优先只改 workspace 文件里的路径，不要更换 workspace 文件本身的位置。