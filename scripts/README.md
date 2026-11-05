# scripts

包装脚本与冒烟。保持**薄**：只调用服务入口 / 测试入口，不承载业务逻辑。

| 脚本 | 作用 |
| :--- | :--- |
| `sync-cursor-skills.ps1` | 同步 `.cursor/skills/` 的 Matt Pocock skill 副本 |
| `list_mcp_tools.py` | 列出已实现服务的 MCP 工具面（实现后） |
| `smoke_*.py` | stdio 端到端冒烟 |

命名：Python 脚本用 `snake_case`；Windows 与 Linux 各一份包装时行为保持一致。
