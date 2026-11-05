# tools

本地开发工具（调试桥、本地 MCP 等）。

示例结构（按需创建）：

```text
tools/
  mcp_bridge/
    __init__.py
    __main__.py
    README.md
```

约定：

- 工具通过 `.cursor/mcp.json` 或宿主自己的 MCP 配置注册（stdio）。
- 路径与环境变量用参数 / `HYDRO_AGENT_DATA_DIR` 注入，避免依赖当前工作目录。
