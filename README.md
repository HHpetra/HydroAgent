# HydroAgent

HydroAgent 是一个**纯 MCP 服务**项目：面向水库防洪调度，提供水库资料查询与静库容调洪演算的标准 MCP 工具面。项目使用 **Python** 技术栈，不依赖 Node.js / TypeScript。

## 状态

**设计阶段（Design）** `src/hydro_agent/servers/` 尚未实现。`hydro_agent.core` 已提供 SCS-CN 产流。

## 服务边界

| MCP server | 职责 | 传输 |
| :--- | :--- | :--- |
| `reservoir-info` | 只读资料：特征水位、库容曲线、泄流曲线或闸门几何、入库 / 降雨 / 水位 / 出库 / 开度过程 | stdio |
| `floodroute` | 静库容调洪：无闸 q(Z) 或有闸堰流 / 孔流，调度方案与特征值 | stdio |

两个服务都使用**场景时钟 `as_of`**；工具必须使用资料返回的 `as_of`，不得使用现实当前时间。

## 技术栈

| 依赖 | 版本 | 说明 |
| :--- | :--- | :--- |
| Python | >= 3.11 | `pyproject.toml` `requires-python` |
| uv | 0.12.5（本机验证） | 包管理与虚拟环境；锁文件为 `uv.lock` |
| MCP SDK | `mcp>=2.3,<3` | Python SDK；精确版本由 `uv.lock` 锁定 |
| pytest | >= 8.0 | 单元与集成测试 |
| ruff / mypy | dev 依赖 | 格式、lint 与类型检查 |

## 安装与测试

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run mypy src
```

MCP 服务实现后，入口约定为：

```bash
uv run python -m hydro_agent.servers.reservoir_info
uv run python -m hydro_agent.servers.floodroute
```

宿主通过 stdio 启动（按实际仓库路径替换）：

```json
{
  "mcpServers": {
    "reservoir-info": {
      "command": "uv",
      "args": [
        "run",
        "--directory",
        "C:/Agent/Agent05-HydroAgent",
        "python",
        "-m",
        "hydro_agent.servers.reservoir_info"
      ]
    }
  }
}
```

## 目录说明

| 路径 | 用途 |
| :--- | :--- |
| `src/hydro_agent/` | Python 包 |
| `src/hydro_agent/core/` | 共享领域逻辑：时序、曲线、调洪算法 |
| `src/hydro_agent/servers/` | MCP 服务入口：reservoir-info / floodroute |
| `tests/` | pytest 测试与 `fixtures/` 样例数据 |
| `docs/` | 架构说明 + `adr/` + `agents/`（长期） |
| `docs/research/` | 调研与可行性 |
| `scripts/` | 包装脚本与冒烟 |
| `tools/` | 本地 MCP 调试桥 |
| `GLOSSARY.md` | 领域语言 |
| `AGENTS.md` | 协作者与 Agent 的工程约定 |

## 延伸阅读

- [`docs/README.md`](docs/README.md) — 架构说明
- [`docs/agents/domain.md`](docs/agents/domain.md) — 领域文档规则
- [`GLOSSARY.md`](GLOSSARY.md) — 领域语言
- [`AGENTS.md`](AGENTS.md) — 工程约定

## 许可

MIT. Copyright (c) 2026 HydroAgent Contributors.
