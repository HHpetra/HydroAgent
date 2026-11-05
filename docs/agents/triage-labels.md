# 分诊标签（Triage Labels）

`triage` skill 用这五个角色把**外来** issue（bug 报告、功能请求）整理成 agent 可接手的票。标签字符串与角色同名，不做本地化。

首次在仓库启用时建一次（已存在的会报错，忽略即可）：

```powershell
'needs-triage','needs-info','ready-for-agent','ready-for-human','wontfix' |
  ForEach-Object { gh label create $_ 2>$null }
```

| 标签 | 含义 |
| :--- | :--- |
| `needs-triage` | 刚进来，还没判断 |
| `needs-info` | 缺关键信息，等提出者补充 |
| `ready-for-agent` | 范围清晰，agent 可直接开工 |
| `ready-for-human` | 需要人做决定或人手操作（例如发版、第三方账号配置） |
| `wontfix` | 不做；理由与替代方案写进 issue |

只用于外部进来的 issue。`to-tickets` 拆出来的票本身已 agent-ready，**不要**再过一遍分诊。

> 未安装 `triage` skill 时可忽略本文件；`AGENTS.md` 的 Agent skills 块里也不用写 Triage labels 一段。
