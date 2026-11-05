# tests

pytest 测试与样例数据。

- 单元测试优先与实现同目录 `test_*.py`；跨包 e2e / 冒烟放本目录。
- `fixtures/` 放样例水库、过程线与期望结果。
- 运行：`uv run pytest`。
- fixtures 不要放敏感数据或超大二进制。
