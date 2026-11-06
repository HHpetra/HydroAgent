# 架构说明

HydroAgent 是纯 Python MCP 服务项目：不包含桌面壳、launcher 或 dsh profile。MCP 宿主通过 stdio 加载服务。

- 只写架构和关键原理，不写实现细节与过程性 spec。
- 具体实现以阅读代码为准。

## 分层

| 层 | 位置 | 职责 |
| :--- | :--- | :--- |
| MCP 服务 | `src/hydro_agent/servers/reservoir_info.py` | 只读水库资料与过程序列 |
| MCP 服务 | `src/hydro_agent/servers/floodroute.py` | 静库容调洪与调度方案计算 |
| 领域核心 | `src/hydro_agent/core/` | SCS-CN 产流、时序、曲线、水力与调洪算法（共享） |
| 数据 | `tests/fixtures/`、`HYDRO_AGENT_DATA_DIR` | 可复现样例与运行资料 |
| 宿主 | 任意 MCP 客户端 | stdio 启动、工具调用与结果展示 |

## 关键数据流

1. 宿主先调用 `reservoir-info` 的 `list_reservoirs` / `get_flood_routing_inputs`。
2. 宿主把返回的 `as_of`、`z0`、入库过程、库容曲线、出流侧（`discharge` 或 `gates`）与实测过程原样传给 `floodroute`。
3. `floodroute` 返回特征值与过程序列；过程图、调度建议与最终汇报由宿主组织。

## 场景时钟

所有工具以 `as_of` 为时间基准，不使用系统时间；`as_of` 从资料工具返回并原样传递。

## SCS-CN 产流

流域净雨在 `hydro_agent.core` 中按一场雨的累积降雨计算：先得到累积净雨，再差分出逐时段净雨。初损在整场雨里只消耗一次。按场景时刻截取余雨时，先算完整场雨，再取该时刻之后的增量。直接相连的不透水面积默认不扣损。

## 边界

- stdout 只写 MCP 协议消息，日志写 stderr。
- 两个服务共享 `hydro_agent.core`，不互相 import。
- 不 fork 官方 MCP Python SDK；不引入 dsh、Node.js 或 TypeScript 工程。

## 同目录其他文档

| 路径 | 用途 |
| :--- | :--- |
| `adr/` | 架构决策记录：一决策一文件，编号递增（首次写决策时创建） |
| `agents/` | Agent 消费的配置：`domain.md` / `issue-tracker.md` / `triage-labels.md` |
| `research/` | 可行性研究与技术调研 |

> 短期工作（spec、票据）不落在仓库里，走 GitHub Issues；见 [`agents/issue-tracker.md`](agents/issue-tracker.md)。
