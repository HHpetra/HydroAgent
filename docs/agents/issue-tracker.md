# Issue Tracker

本仓用 **GitHub Issues** 作为**唯一**的短期工作载体：spec、票据、bug、任务都在 issue 里。仓库里**不设** `spec/` 目录——短期的东西一旦落成文件，就会在几周后变成"过期文档被当成规范"。

## 约定

- **spec**：由 `/to-spec` 从对话直接发布成 issue，用固定模板（Problem / Solution / User Stories / Implementation Decisions / Testing Decisions / Out of Scope），打 `ready-for-agent` 标签。**不要**落成仓库文件。
- **tickets**：`/to-tickets` 把 spec 拆成 tracer-bullet 子 issue；**阻塞关系要写成原生依赖**（sub-issue / blocked-by），不要只在正文里写"依赖 #12"。
- **取票**：只取阻塞项都已关闭的票（frontier）；一票一个会话，只带那张 issue 进上下文。
- **收口**：实现 + review 后关闭 issue；值得长期保留的结论进 `docs/adr/`，issue 不是决策的最终载体。

## 常用命令

| 动作 | 命令 |
| :--- | :--- |
| 建 issue | `gh issue create --title "…" --body-file body.md --label ready-for-agent` |
| 列队列 | `gh issue list --state open --label ready-for-agent` |
| 看单条 | `gh issue view <n> --comments` |
| 关闭 | `gh issue close <n> --comment "已合入 #<pr>"` |
| 改标签 | `gh issue edit <n> --add-label needs-info` |

阻塞关系：`gh` 没有原生的 sub-issue / blocked-by 子命令，用 `gh api graphql` 或网页设置；团队已有约定时以约定为准。

## 还没有远端时

刚复制出来的项目通常还没建 GitHub 仓库，先建了再用：

```powershell
git init
gh repo create <owner>/<name> --private --source . --remote origin --push
```

确实没法建远端（离线 / 内网）时，临时退到本地 markdown：`.scratch/<feature>/issues/<NN>-<slug>.md`，一票一文件，`Blocked by:` 写编号。它只是**临时容器**，有 remote 后迁回 issue，不要长期双轨。
