# HydroAgent - Agent 指南与项目规范

## 1. 项目简介

- **项目名称**：HydroAgent
- **定位**：纯 MCP 服务项目，不是 dsh 插件，也不包含桌面 / 浏览器产品壳。
- **愿景**：让 Agent 通过标准 MCP 工具完成水库资料查询、调洪演算与防洪情势分析。

---

## 2. 当前项目阶段

- **当前阶段**：设计阶段 (Design)。仓库目前只有 Python 项目骨架，`src/hydro_agent/servers/` 尚未实现。
- **核心目标**：
  1. 定义 reservoir-info / floodroute 两个 MCP 服务的工具面与数据契约。
  2. 实现静库容调洪核心算法与可复现的样例数据。
  3. 保持 stdio 独立运行，可被任意 MCP 宿主加载。

---

## 3. 目录与文档规范

项目严格遵循以下目录规范：

| 目录路径 | 用途说明 | 保留策略 |
| :--- | :--- | :--- |
| `GLOSSARY.md` | **领域语言**：本项目特有概念的唯一定义 | **长期保留**；只放词义，不放实现细节与决策 |
| `docs/` | 面向用户与 Agent 的**架构与关键原理**介绍 | **长期保留**；只写架构与原理，不写实现细节与过程性 spec |
| `docs/adr/` | **架构决策记录（ADR）**：一决策一文件 | **长期保留**；只记关键决定 |
| `docs/agents/` | Agent 消费的配置（领域文档规则 / issue tracker / 分诊标签） | **长期保留** |
| `docs/research/` | 可行性研究与技术调研文档 | 调研阶段产物；结论沉淀进 ADR / `docs/` 后可不再当日常真源 |
| `src/hydro_agent/` | Python 包；`core/` 放共享领域逻辑，`servers/` 放 MCP 入口 | 长期 |
| `tests/` | pytest 测试与 `fixtures/` 样例数据 | 长期 |
| `scripts/` | 包装脚本、冒烟、开发工具 | 长期 |
| `tools/` | 本地 MCP / 调试桥等开发工具 | 长期 |

**短期工作不进仓库**：spec、票据、任务全部走 GitHub Issues，见 `docs/agents/issue-tracker.md`。

### 3.1 ADR（长期）

- 路径：`docs/adr/`；**一决策一文件**：`0001-slug.md`、`0002-slug.md`。
- **目录懒创建**：第一次真要写 ADR 时再建 `docs/adr/`；没有决策就不建空目录。
- 编号扫现存最大值 +1；允许缺号，**不重新编号**。
- **三条门槛同时满足才写**：难以回退 / 没有上下文会让后人困惑 / 真有过取舍。缺一条就跳过。
- 正文 = 标题 + 1–3 句（背景、决定、为什么）；`后果`、`被否决的方案` 只在真有价值时展开。
- 取代：新增一条并写 `**状态**：superseded by ADR-NNNN`；旧文件正文不改，也不另建 `deprecated.md`。
- **不记录**：实现细节、局部重构、显然只有一个选项的事。

### 3.2 Spec 与票据（短期）

- **不发文件**：spec 由 `/to-spec` 发布成 GitHub Issue，票据由 `/to-tickets` 拆成带阻塞边的子 issue。
- 仓库里不设 `spec/` 目录；`.scratch/` 只在没有 remote 时作临时容器。
- Agent 不应把过期 issue / 旧文档当作现行规范；以 `GLOSSARY.md`、ADR、`docs/` 与代码为准。

### 3.3 `docs/`（架构导读）

- 介绍项目主要功能与分层，方便快速理解架构。
- **只记录架构和关键原理**；具体实现以**阅读代码**为准。

### 3.4 调研文档（`docs/research/`）

- 建议包含：背景、核心诉求、方案对比、可行性、推荐结论、下一步。
- 已拍板结论同步写入 `docs/adr/`；面向读者的架构叙述写入 `docs/`。

---

## 4. 技术栈

| 项 | 选型 | 说明 |
| :--- | :--- | :--- |
| 运行时 | Python >= 3.11 | 根 `pyproject.toml` `requires-python` |
| 包管理 | uv | **只允许 uv**；锁文件 `uv.lock`；不要混用 pip / poetry / conda |
| 语言 | Python 3.11+ | 使用类型注解；Node / TypeScript 不属于本项目技术栈 |
| MCP | `mcp` Python SDK（实现时在 `pyproject.toml` 钉版本） | 不 fork 上游 |
| 传输 | stdio 优先 | stdout 只走 MCP 协议 |
| 测试 | pytest | 与实现同目录或 `tests/` 下 `test_*.py` |
| 静态检查 | ruff + mypy | 提交前运行 |

---

## 5. 构建与测试

```bash
uv sync --extra dev   # 创建 .venv 并安装依赖
uv run pytest         # 单元测试 + 集成测试
uv run ruff check .   # lint
uv run mypy src       # 类型检查
```

服务实现后，入口约定为 `python -m hydro_agent.servers.<server>`。

| 环境变量 | 作用 |
| :--- | :--- |
| `HYDRO_AGENT_DATA_DIR` | 水库资料目录（实现时使用，禁止写死绝对路径） |

---

## 6. MCP 服务约定

1. **服务拆分**：`reservoir_info` 只读资料；`floodroute` 只做计算。共享代码放 `hydro_agent.core`，不复制。
2. **stdio 传输**：stdout 只写 MCP 协议消息，日志用 `logging` 写 stderr，避免污染协议流。
3. **场景时钟**：所有工具显式接收并返回 `as_of`；禁止用现实当前时间。
4. **单位固定**：水位 m，库容 m³，流量 m³/s，降雨 mm；时间 RFC3339。
5. **工具面**：参数用类型化模型（MCP tool schema / pydantic）；结果用 JSON 文本；错误以工具错误返回，不把异常抛给宿主。
6. **确定性**：样例数据放 `tests/fixtures/`，不隐式访问网络。
7. **不 fork 上游**：直接依赖官方 MCP Python SDK；不引入 dsh、Node.js 或 TypeScript 工程。
8. **单一包管理器**：只用 uv；不在子目录单独创建虚拟环境或混用包管理器。

---

## 7. Agent 协作指引

1. **文档落点**：关键决定 → `docs/adr/`；当前开发 → GitHub Issues（`docs/agents/issue-tracker.md`）；架构导读 → `docs/`；调研 → `docs/research/`；术语 → `GLOSSARY.md`。
2. **真源优先级**：长期行为以 **`GLOSSARY.md` + ADR + `docs/` + 代码** 为准。
3. **不做长期 spec 库**：不要把完整产品 spec 堆进仓库长期保存。
4. **代码检索**：优先使用语言服务器 / IDE 符号跳转；命令行检索用 `rg` 并限定路径与文件类型，不要一上来全局关键词扫仓。
5. **动手前先读**：进入某区域前先读 `GLOSSARY.md` 与相关 ADR；与既有 ADR 冲突时显式标注（"与 ADR-NNNN 冲突，但值得重开，因为…"），不静默覆盖。

---

## 8. 代码检索

- 找符号定义 / 调用关系：优先用 IDE、Pyright / mypy 语言服务或 `rg -t py`。
- 跨目录检索时限定路径与 glob，避免把 `.venv/`、缓存和样例数据当作源码。
- 不要引入 Node.js 代码检索工具作为默认工作流。

---

## 9. 远程与协作

| 项 | 值 |
| :--- | :--- |
| 默认分支 | `main` |
| 探索分支 | `invest`（调研草稿，不作为现行真源） |
| 远程 | https://github.com/HHpetra/HydroAgent.git |

提交前建议先同步，避免多机冲突。

---

## 10. Agent skills（Matt Pocock 流程）

由 `setup-matt-pocock-skills` 生成；本模板已预置默认值，可直接改。

### Issue tracker

本地 markdown 只作**无 remote 时的临时容器**；正常情况 spec 与票据都是 GitHub Issue，阻塞关系用原生依赖。见 `docs/agents/issue-tracker.md`。

### Triage labels

`needs-triage` / `needs-info` / `ready-for-agent` / `ready-for-human` / `wontfix`，只用于外部进来的 issue。见 `docs/agents/triage-labels.md`。

### Domain docs

单上下文：根 `GLOSSARY.md` + `docs/adr/`；动手前先读，冲突显式标注。见 `docs/agents/domain.md`。

---

## 11. 许可

MIT. Copyright (c) 2026 HydroAgent Contributors.
