# vareffect agent harness

The **main engineer** (the session agent) owns planning, implementation, tests, documentation and integrating findings, and is the only editor. Three read-only agents support it; skills carry the workflow and domain knowledge. Everything here is native Claude Code configuration — edit files directly, there is no generator. Adapted from the LIA harness (2026-10-07): LIA's clinical, security and UI reviewers, its ledgers and its wiki rules do not apply to a library and were left out.

| Path | Contents |
| --- | --- |
| `INVARIANTS.md` | Binding, maintainer-owned charter (V-1 … V-9, plus its known conflicts with current code). Agents quote and flag conflicts; only an authorized human revision changes it. |
| `agents/` | `vareffect-reviewer` (engineering), `vareffect-bioinformatics` (annotation correctness, pinned to Opus), `vareffect-researcher` (targeted lookups). Shared obligations and the report format are preloaded from `vareffect-findings`. |
| `skills/` | Workflow: `vareffect-ground`, `vareffect-review`, `vareffect-sequence`. Support: `vareffect-findings`, `vareffect-verification`, `concordance-measurement`. Conventions: `idiomatic-rust` (the full Rust rule set; `CLAUDE.md` keeps only the hard rules). Reference: `citation-verification`. |
| `hooks/pin-reviewer-models.py` | Blocks spawning a reviewer below its model floor. |
| `settings.json` | Project model/effort defaults, charter edit deny, secret-file read deny, the hook. |

`scripts/vareffect_harness.py` routes changes (`route`), freezes reviewed bytes (`snapshot`) and checks this structure (`check`); `make harness-check` runs the check and its tests. Task records live in gitignored `.vareffect-tasks/`.

## Workflow

`vareffect-ground` scopes a task, assigns a risk tier (T0 trivial → T3 annotation semantics) and, from T2 up, gets an independent plan review. The main engineer implements and verifies (`vareffect-verification`). Anything that can move annotation output is measured before and after on a pinned corpus (`concordance-measurement`). `vareffect-review` freezes the reviewed bytes, writes one review packet and dispatches fresh reviewers in parallel. Non-protected findings are re-checked by continuing the reviewer that raised them; protected findings need a fresh `vareffect-bioinformatics` instance to confirm their closure. Maintainer decisions — accepting a VEP divergence, transcript-source policy, breaking changes, releases — go to the user under "Decisions for you". There is no findings ledger: pre-existing issues are reported, and become GitHub issues only at the user's direction.

Limits: two plan revision rounds; two repair attempts per finding (the second needs new evidence) and two repair rounds per task; stop on no progress. Unavailable review or checks are reported as incomplete, never as a pass. Commits and releases each need their own authorization.

## Model policy

| Agent | Model / effort | Why |
| --- | --- | --- |
| Main engineer | Opus / high (project default) | Owns design and integration. |
| `vareffect-bioinformatics` | **Opus / high, pinned** | Annotation correctness feeds clinical interpretation downstream; it runs only when routing or semantics call for it. |
| `vareffect-reviewer` | Sonnet / high; spawn with `model: opus` for T2/T3 | Routine engineering review is well served by Sonnet; protected or cross-cutting work gets Opus. |
| `vareffect-researcher` | Sonnet / medium | Bounded lookups. |

A spawn-time `model` overrides an agent's frontmatter, and environment settings can override or remap both, so the pin is enforced by the `PreToolUse` hook: `vareffect-bioinformatics` accepts only Opus-class models (the `opus`/`fable` aliases or a full Opus/Fable id, including Bedrock and Vertex forms), and the hook also blocks it when `CLAUDE_CODE_SUBAGENT_MODEL_FORCE` is set without an Opus `CLAUDE_CODE_SUBAGENT_MODEL`, or `ANTHROPIC_DEFAULT_OPUS_MODEL` maps the alias to a non-Opus id. No reviewer may run on Haiku. Any failure to run the hook — a missing `python3`, a crash, no hook found — exits 2 and blocks every Agent call, fail closed by design. Hooks are fixed when a session starts: after changing the hook or its command, start a new session. Never substitute a cheaper model when a required one is unavailable — report the limit.

## Read-only controls

Every agent allows exactly `Read, Grep, Glob, WebSearch, WebFetch` and denies `Bash, Edit, Write, NotebookEdit, Agent, Skill, mcp__*`. Skills are preloaded through `skills:`. Do not add `memory` (it grants Write and Edit), `hooks`, delegated tools or MCP servers to an agent; the tool inventory, not the words "read-only", is the boundary.

The charter is protected by `Edit(/.claude/INVARIANTS.md)`, which covers every built-in editing tool including Write; it anchors at the directory the session starts in, so start sessions at the repository root, and it does not stop subprocess writes. Credential files (`.env` variants, SSH keys, cargo, cloud and Docker credentials, certificate and key files) are read-denied with absolute `//**/` rules; `snapshot` also refuses them and sequencing data files. These are governance aids, not a security boundary against someone who changes the configuration. `make harness-check` verifies structure statically; it proves nothing about a live session.

## Smoke procedure (after a Claude Code upgrade or harness change)

1. `make harness-check`; record `claude --version`.
2. `claude -p --agent vareffect-reviewer --permission-mode dontAsk --no-session-persistence --verbose --output-format stream-json` with a harmless prompt asking for a write to a fresh temporary path outside the repository. In the init event confirm the three agents and all skills are registered once and the tools are exactly the five read-only ones; confirm no file was created.
3. In a fresh `claude -p --output-format stream-json --verbose` session (with session persistence), make two Agent calls: `vareffect-bioinformatics` with `model: sonnet` (the hook must block it) and `vareffect-reviewer` without a model, asking it to quote its report format's verdict line from context without tools (proves preloading). Continue that reviewer with `SendMessage` and ask it to list its tools.
4. Record each result separately; never summarize a partial run as a pass.

Not yet run for this repository.
