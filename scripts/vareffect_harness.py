#!/usr/bin/env python3
"""Harness structure checking, conservative Git change routing and review snapshots (Python 3.11+).

Adapted from the LIA harness. `route` gives the minimum reviewers, verification categories and tier floor
for a comparison; `snapshot` freezes the reviewed worktree bytes as one Git tree id; `check` validates the
agents, skills, settings and hook statically.
"""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote

# Agent -> (model, effort, preloaded skills). The hook's Opus pins must agree with this table.
AGENT_POLICY = {
    "vareffect-reviewer": ("sonnet", "high", {"vareffect-findings", "idiomatic-rust"}),
    "vareffect-bioinformatics": ("opus", "high", {"vareffect-findings", "citation-verification", "concordance-measurement"}),
    "vareffect-researcher": ("sonnet", "medium", {"citation-verification"}),
}
AGENT_KEYS = {"name", "description", "tools", "disallowedTools", "model", "effort", "skills"}
REQUIRED_SKILLS = {
    "vareffect-ground", "vareffect-review", "vareffect-sequence", "vareffect-findings", "vareffect-verification",
    "concordance-measurement", "citation-verification", "idiomatic-rust",
}
AGENT_TOOLS = {"Read", "Grep", "Glob", "WebSearch", "WebFetch"}
AGENT_DENY = {"Bash", "Edit", "Write", "NotebookEdit", "Agent", "Skill", "mcp__*"}
SPECIALIST = "vareffect-bioinformatics"
ENGINEERING = "vareffect-reviewer"
HOOK = ".claude/hooks/pin-reviewer-models.py"
# Project copy first; the working tree's only when the launched project has none and both are checkouts of the
# same repository (a session that entered a worktree before the hook was merged). Any failure to run it, including
# a missing python3 or no hook found, blocks the spawn.
HOOK_COMMAND = (
    f'f="$CLAUDE_PROJECT_DIR/{HOOK}"; '
    'if [ ! -f "$f" ]; then '
    'a="$(git rev-parse --path-format=absolute --git-common-dir 2>/dev/null)"; '
    'b="$(git -C "$CLAUDE_PROJECT_DIR" rev-parse --path-format=absolute --git-common-dir 2>/dev/null)"; '
    'if [ "${a#/}" != "$a" ] && [ "$a" = "$b" ]; '
    f'then f="$(git rev-parse --show-toplevel 2>/dev/null)/{HOOK}"; fi; fi; '
    'python3 "$f" || exit 2'
)
# Absolute forms: relative rules do not reach files outside the project directory, and permission rules have
# no negation, so `.env.example` stays readable by not being listed.
SECRET_READ_DENIES = {
    f"Read(//**/{name})"
    for name in (
        ".env", "*.env", ".env.local", ".env.*.local", ".env.production", ".env.development", ".env.staging",
        ".env.test", ".env.prod", ".env.dev", ".env.ci", ".env.bak", ".env.backup", ".env.old", ".env.qa",
        ".env.uat", ".env.docker", ".envrc", ".netrc", ".pgpass", ".git-credentials", ".npmrc", ".pypirc",
        ".ssh/**", "id_rsa", "id_ed25519", "id_ecdsa", "id_dsa", ".doppler/**", ".cargo/credentials*", ".kube/config", ".config/gcloud/**", ".config/gh/hosts.yml",
        ".aws/credentials", ".docker/config.json", "*credentials*.json", "*.pem", "*.key", "*.p12", "*.pfx", "*.jks",
    )
}
TASK_DIR = ".vareffect-tasks/"
# Machine-local paths a harness document may name although they are absent from a fresh clone.
MACHINE_LOCAL = (".claude/settings.local.json", ".claude/plans/", ".claude/worktrees/", TASK_DIR)
MAX_SKILL_DESCRIPTION = 300
ROOT_INSTRUCTIONS_LIMIT = 16 * 1024
CHAIN_LIMIT = 24 * 1024
PATH_PREFIXES = (".claude/", "vareffect/", "vareffect-cli/", "scripts/", ".github/")

# Everything the library does is annotation, and almost every CLI module feeds or emits it; the exempt
# modules only scaffold or check configuration. Tests and fixtures hold recorded ground truth.
ANNOTATION_PATHS = ("vareffect/src/", "vareffect/tests/", "vareffect-cli/src/", "vareffect-cli/vareffect_build.toml")
ANNOTATION_EXEMPT = ("vareffect-cli/src/init.rs", "vareffect-cli/src/check.rs")
DIVERGENCES = "vareffect/vep_divergences.md"
# Code whose change means rebuilding data files before the data-gated suites say anything.
STORE_PATHS = (
    "vareffect-cli/src/builders/", "vareffect-cli/src/setup.rs", "vareffect-cli/src/common.rs",
    "vareffect-cli/vareffect_build.toml", "vareffect/src/types.rs", "vareffect/src/transcript.rs", "vareffect/src/fasta.rs",
    "vareffect/src/chrom.rs",
)
# Public contract surfaces (V-6): every library module re-exports or feeds public types, and the CSQ layout is
# the CLI's output contract.
CONTRACT_PATHS = ("vareffect/src/", "vareffect-cli/src/csq.rs")
# Harness files carry review authority. Only pure style references that decide nothing are exempt.
HARNESS_REFERENCE_ONLY = (".claude/skills/idiomatic-rust/",)
BUILD_CONFIG = {"makefile", "cargo.toml", "cargo.lock", "deny.toml", "rust-toolchain.toml", ".gitignore"}
UNSAFE = re.compile(r"\bunsafe\s*(\{|\(|fn\b|impl\b|trait\b|extern\b)")
# Files snapshot refuses to store or write into a review diff unless each path is explicitly allowed.
SENSITIVE_FILE = re.compile(
    r"(^|/)(\.env(\.(?!example$)[^/]*)?|[^/]*\.env|\.envrc|\.netrc|\.pgpass|\.git-credentials|\.npmrc|\.pypirc"
    r"|id_(rsa|ed25519|ecdsa|dsa)|settings\.local\.json|[^/]*credentials[^/]*\.json)$"
    r"|(^|/)\.(ssh|doppler)/|(^|/)\.kube/config$|(^|/)\.cargo/credentials|(^|/)\.config/(gcloud/|gh/hosts\.yml$)"
    r"|(^|/)\.aws/credentials$|(^|/)\.docker/config\.json$"
    r"|\.(vcf|gvcf|bcf|bam|cram|sam|fastq|fq|pem|key|p12|pfx|jks)(\.(gz|bgz|bz2|xz|zst))?$",
    re.I,
)


class HarnessError(Exception):
    """Invalid repository, comparison, or harness input."""


# Inherited pathspec modes would turn the snapshot's exclude pathspecs into literal names.
PATHSPEC_ENV = ("GIT_LITERAL_PATHSPECS", "GIT_GLOB_PATHSPECS", "GIT_NOGLOB_PATHSPECS", "GIT_ICASE_PATHSPECS")


def git(root, *args, env=None):
    """Run a Git command without external diff drivers, optional locks or inherited pathspec modes."""
    env = {k: v for k, v in (env or os.environ).items() if k not in PATHSPEC_ENV}
    env["GIT_OPTIONAL_LOCKS"] = "0"
    result = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, env=env, check=False
    )
    if result.returncode:
        raise HarnessError(result.stderr.decode(errors="replace").strip())
    return result.stdout


def commit(root, ref):
    if not ref or ref.startswith("-"):
        raise HarnessError(f"invalid revision: {ref!r}")
    return git(root, "rev-parse", "--verify", "--end-of-options", ref + "^{commit}").decode().strip()


def comparison(root, expression):
    """Resolve a base, two-dot endpoint comparison, or three-dot merge-base comparison."""
    if "..." in expression:
        parts = expression.split("...")
        mode = "merge-base"
    elif ".." in expression:
        parts = expression.split("..")
        mode = "endpoints"
    else:
        parts = [expression, "HEAD"]
        mode = "base-to-HEAD"
    if len(parts) != 2 or not all(parts):
        raise HarnessError("comparison needs explicit valid endpoints")
    left, right = (commit(root, ref) for ref in parts)
    start = left
    if mode == "merge-base":
        bases = git(root, "merge-base", "--all", left, right).decode().splitlines()
        if len(bases) != 1:
            raise HarnessError("comparison must have exactly one merge base")
        start = bases[0]
    return {"expression": expression, "mode": mode, "left": left, "right": right, "start": start}


def diff_paths(raw):
    """Parse NUL-delimited name-status output, retaining both sides of renames/copies."""
    fields = raw.split(b"\0")
    if fields[-1] == b"":
        fields.pop()
    paths = []
    index = 0
    while index < len(fields):
        status = fields[index].decode("ascii")
        count = 2 if status.startswith(("R", "C")) else 1
        if not status or index + count >= len(fields):
            raise HarnessError("malformed Git name-status output")
        for name in fields[index + 1:index + count + 1]:
            paths.append((os.fsdecode(name), status))
        index += count + 1
    return paths


def collect_changes(root, expression):
    comp = comparison(root, expression)
    args = ("diff", "--no-ext-diff", "--no-textconv", "--name-status", "-z", "--find-renames")
    selections = [
        ("comparison", (*args, comp["start"], comp["right"], "--")),
        ("staged", (*args, "--cached", "HEAD", "--")),
        ("unstaged", (*args, "--")),
    ]
    paths = {}
    for source, command in selections:
        for path, status in diff_paths(git(root, *command)):
            paths.setdefault(path, []).append({"source": source, "status": status})
    for raw in git(root, "ls-files", "--others", "--exclude-standard", "-z").split(b"\0"):
        if raw:
            paths.setdefault(os.fsdecode(raw), []).append({"source": "untracked", "status": "?"})
    return comp, paths


def route_path(path):
    """Canonical minimum path rules; semantic risk must add to these, never subtract."""
    p = path.lower()
    name = Path(p).name
    reviewers = {ENGINEERING}
    verification = set()
    reasons = ["every changed file needs engineering review"]
    by_role = {}
    floor = "T1"

    def add(category, reason, role=None, tier=None):
        nonlocal floor
        verification.add(category)
        reasons.append(reason)
        if role:
            reviewers.add(role)
            by_role.setdefault(role, []).append(reason)
        if tier and tier > floor:
            floor = tier

    harness = p.startswith((".claude/", "scripts/vareffect_harness", "scripts/test_vareffect_harness")) or name in {"agents.md", "claude.md"}
    if harness:
        verification.add("harness")
        if not p.startswith(HARNESS_REFERENCE_ONLY):
            add("harness", "harness authority: charter, review rules, routing, tool grants, hooks or settings", SPECIALIST, "T2")
    code = p.endswith(".rs") or p.endswith(".toml")
    if p.startswith(ANNOTATION_PATHS) and not p.startswith(ANNOTATION_EXEMPT):
        add("concordance" if code or p.startswith("vareffect/tests/") else "documentation", "annotation logic, its inputs, outputs or recorded ground truth", SPECIALIST, "T2")
    if p == DIVERGENCES:
        add("documentation", "documented VEP divergences (V-4)", SPECIALIST, "T2")
    if p.startswith(STORE_PATHS):
        add("store-rebuild", "transcript store or genome build and load: rebuild the data before the data-gated suites")
    if p.startswith(CONTRACT_PATHS):
        add("public-contract", "public API, JSON or CSQ contract surface (V-6): check semver and CHANGELOG", tier="T2")
    if p == ".github/workflows/release.yml":
        add("release", "release workflow: publishing and its credentials", tier="T2")
    elif p.startswith(".github/workflows/"):
        add("build-config", "CI workflow: what gates every change", tier="T2")
    if name in {"cargo.toml", "cargo.lock"}:
        # Dependencies sit in the annotation path (coitrees overlap order, rmp-serde store decode, serde_json
        # output, mmap, decompression) and the release profile keeps unwinding for catch_unwind.
        add("store-rebuild", "dependency or profile change: can change store compatibility, overlap order or panic survival", SPECIALIST, "T2")
        verification.add("concordance")
    if p.endswith((".md", ".mdx")):
        verification.add("documentation")
    if p.endswith(".rs") or name in {"cargo.toml", "cargo.lock"}:
        verification.add("rust")
    if name in {"cargo.lock", "deny.toml"}:
        verification.add("supply-chain")
    if p.startswith((".github/", "scripts/")) or name in BUILD_CONFIG:
        verification.add("build-config")
    if not verification:
        verification.add("scope-specific")
    return {"reviewers": sorted(reviewers), "verification": sorted(verification), "reasons": reasons, "by_role": by_role, "tier_floor": floor}


def content_texts(root, path, revisions):
    """The file's text in the worktree and at each revision that has it; (texts, unreadable)."""
    texts = []
    unreadable = False
    for revision in dict.fromkeys(revisions):
        listed = git(root, "ls-tree", "-z", "--name-only", revision, "--", path).split(b"\0")
        if os.fsencode(path) not in listed:
            continue
        try:
            texts.append(git(root, "show", f"{revision}:{path}").decode(errors="replace"))
        except HarnessError:
            unreadable = True
    if (root / path).exists():
        try:
            texts.append((root / path).read_text(errors="replace"))
        except OSError:
            unreadable = True
    return texts, unreadable


def route_content(root, path, rules, revisions=()):
    """Add content-based minima a path cannot show: Rust that adds, removes or sits beside `unsafe` code."""
    if not path.endswith(".rs"):
        return rules
    texts, unreadable = content_texts(root, path, revisions)
    if unreadable or any(UNSAFE.search(text) for text in texts):
        rules["verification"] = sorted({*rules["verification"], "unsafe"})
        rules["reasons"].append("file could not be read for its unsafe check" if unreadable else "file contains unsafe code: document and re-check its invariants")
        rules["tier_floor"] = max(rules["tier_floor"], "T2")
    return rules


def route(root, expression, verbose=False):
    """Route changes to minimum reviewers; compact by default, per-file detail when verbose."""
    comp, changes = collect_changes(root, expression)
    head = commit(root, "HEAD")
    revision = None if comp["right"] == head else comp["right"]
    if revision:
        # The snapshot diff runs from the comparison start to the worktree, so commits after `right` count too.
        args = ("diff", "--no-ext-diff", "--no-textconv", "--name-status", "-z", "--find-renames", revision, head, "--")
        for path, status in diff_paths(git(root, *args)):
            changes.setdefault(path, []).append({"source": "after-comparison", "status": status})
    revisions = (comp["start"], comp["right"], head)
    files = [{"path": path, "changes": changes[path], **route_content(root, path, route_path(path), revisions)} for path in sorted(changes)]
    why = {}
    for f in files:
        for role in f["reviewers"]:
            entry = why.setdefault(role, {"files": 0, "examples": [], "reasons": set()})
            entry["files"] += 1
            if len(entry["examples"]) < 3:
                entry["examples"].append(f["path"])
            entry["reasons"].update(f["by_role"].get(role, [f["reasons"][0]]))
    result = {
        "comparison": comp,
        "file_count": len(files),
        "minimum_reviewers": sorted({r for f in files for r in f["reviewers"]}),
        "tier_floor": max((f["tier_floor"] for f in files), default="none"),
        "verification_categories": sorted({v for f in files for v in f["verification"]}),
        "why": {role: {**entry, "reasons": sorted(entry["reasons"])} for role, entry in sorted(why.items())},
        "limits": "Path and content minima only. Add reviewers and raise the tier for semantic effects (T3 when output is meant to move); this is not review or approval.",
    }
    if verbose:
        result["files"] = files
    return result


def snapshot(root, base=None, diff_out=None, allow=()):
    """Identify the exact worktree bytes (tracked + untracked, non-ignored) as one Git tree id.

    Works on a copy of the index, so the real index, refs and files are untouched. Untracked
    files become unreferenced blobs that garbage collection later removes, so paths that look
    like secrets or sequencing data are refused (untracked, modified or intent-to-add before storing; any
    change but a deletion after) unless each is listed in `allow`. With `diff_out`, writes the binary diff
    from the comparison start (or HEAD) to the tree: the durable record of what was reviewed. A deleted
    sensitive file appears there without its old bytes, so removing committed real data stays reviewable.
    """
    allow = set(allow)
    head = commit(root, "HEAD")
    head_tree = git(root, "rev-parse", "HEAD^{tree}").decode().strip()

    def refuse(paths, when):
        sensitive = sorted(p for p in paths if SENSITIVE_FILE.search(p) and p not in allow)
        if sensitive:
            raise HarnessError(f"{when} files look like secrets or sequencing data; remove or ignore them, or pass --allow <path> for each confirmed-synthetic fixture: {sensitive}")

    def names(*args, env=None):
        return [os.fsdecode(p) for p in git(root, *args, env=env).split(b"\0") if p]

    refuse(names("ls-files", "--others", "--exclude-standard", "-z"), "untracked")
    index = Path(git(root, "rev-parse", "--path-format=absolute", "--git-path", "index").decode().strip())
    sparse = git(root, "config", "--bool", "--default", "false", "core.sparseCheckout").decode().strip() == "true"
    with tempfile.TemporaryDirectory() as tmp:
        temp_index = Path(tmp) / "index"
        env = dict(os.environ, GIT_INDEX_FILE=str(temp_index))
        if index.is_file():
            # Keep the index's mtime: Git re-reads entries modified as late as the index write ("racy git"),
            # and a fresh mtime would make such an edit look clean and drop it from the snapshot. Bytes and
            # mtime come from one open file, since Git replaces the index by rename.
            with open(index, "rb") as source, open(temp_index, "wb") as copy:
                shutil.copyfileobj(source, copy)
                stat = os.fstat(source.fileno())
            os.utime(temp_index, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        else:
            git(root, "read-tree", head, env=env)
        # Assume-unchanged and skip-worktree bits would hide worktree bytes, modes and deletions from `add`:
        # clear them in the copy, except for files a sparse checkout legitimately leaves absent.
        entries = [(line[:1], line[2:]) for line in git(root, "ls-files", "-v", "-z", env=env).decode(errors="replace").split("\0") if line]
        assumed = [p for tag, p in entries if tag.islower()]
        skipped = [p for tag, p in entries if tag in "Ss" and not (sparse and not os.path.lexists(root / p))]
        if assumed:
            git(root, "update-index", "--no-assume-unchanged", "--", *assumed, env=env)
        if skipped:
            git(root, "update-index", "--no-skip-worktree", "--", *skipped, env=env)
        # Refuse before `add` stores new bytes: modified, retyped and intent-to-add files (deletions expose nothing).
        refuse(names("diff", "--name-only", "-z", "--diff-filter=d", env=env), "modified or intent-to-add")
        git(root, "add", "--all", "--", ".", env=env)
        tree = git(root, "write-tree", env=env).decode().strip()
    start = comparison(root, base)["start"] if base else head
    changed = names("diff", "--name-only", "-z", "--no-renames", "--diff-filter=d", start, tree, "--")
    refuse(changed, "added or modified")
    # A rename carries the old content under the new name, so a sensitive source is refused like an addition.
    fields = names("diff", "--name-status", "-z", "--find-renames", start, tree, "--")
    renamed_from, index = [], 0
    while index < len(fields):
        if fields[index][:1] in "RC":
            renamed_from.append(fields[index + 1])
        index += 3 if fields[index][:1] in "RC" else 2
    refuse(renamed_from, "renamed")
    deleted = names("diff", "--name-only", "-z", "--no-renames", "--diff-filter=D", start, tree, "--")
    withheld = sorted(p for p in deleted if SENSITIVE_FILE.search(p) and p not in allow)
    result = {"head": head, "tree": tree, "matches_head": tree == head_tree}
    if diff_out:
        # Withheld deletions are kept out of rename pairing (whose hunks would quote old lines) and written
        # separately without their old bytes.
        args = ("diff", "--binary", "--no-ext-diff", "--no-textconv")
        pathspec = ("--", ".", *(f":(exclude,literal){p}" for p in withheld))
        if set(withheld) & set(names("diff", "--name-only", "-z", "--no-renames", start, tree, *pathspec)):
            raise HarnessError("the exclude pathspec did not withhold a sensitive deletion; no diff written")
        diff = git(root, *args, "--find-renames", start, tree, *pathspec)
        if withheld:
            diff += git(root, *args, "--no-renames", "--irreversible-delete", start, tree, "--", *(f":(literal){p}" for p in withheld))
        Path(diff_out).parent.mkdir(parents=True, exist_ok=True)
        Path(diff_out).write_bytes(diff)
        result["diff"] = {"from": start, "to": tree, "path": str(diff_out), "bytes": len(diff)}
    result["withheld_deletions"] = withheld
    result["allowed_sensitive"] = sorted(allow & set(changed + deleted))
    result["limits"] = "Ignored files (e.g. .env, stores, build output) are excluded; unreferenced blobs persist until git gc."
    return result


def frontmatter(path):
    """Parse the scalar, folded-scalar and string-list YAML used by agent and skill headers."""
    content = path.read_text()
    if not content.startswith("---\n") or "\n---\n" not in content[4:]:
        raise HarnessError(f"missing frontmatter: {path}")
    header, body = content[4:].split("\n---\n", 1)
    result = {}
    key = None
    for line in header.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if line.startswith((" ", "\t")):
            if key is None:
                raise HarnessError(f"invalid frontmatter: {path}")
            item = line.strip()
            if isinstance(result[key], list) and item.startswith("- "):
                result[key].append(item[2:].strip().strip('"\''))
            elif isinstance(result[key], str):
                result[key] = f"{result[key]} {item}".strip()
            else:
                raise HarnessError(f"invalid frontmatter list in {path}: {key}")
            continue
        if ":" not in line:
            raise HarnessError(f"invalid frontmatter: {path}")
        key, value = line.split(":", 1)
        if key in result:
            raise HarnessError(f"duplicate frontmatter key in {path}: {key}")
        value = value.strip()
        if value in {">-", ">", "|", "|-"}:
            result[key] = ""
        elif value == "":
            result[key] = []
        else:
            result[key] = value.strip('"\'')
    return result, body


def load_hook(root):
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("pin_reviewer_models", root / HOOK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def instruction_files(root):
    names = {os.fsdecode(p) for p in git(root, "ls-files", "--cached", "--others", "--exclude-standard", "-z").split(b"\0") if p}
    return {Path(p) for p in names if Path(p).name in {"AGENTS.md", "CLAUDE.md"} and (root / p).is_file()}


def harness_documents(root):
    """Maintained instruction documents whose links and path references must resolve."""
    docs = [root / ".claude/README.md", root / ".claude/INVARIANTS.md"]
    docs += sorted((root / ".claude/agents").glob("*.md"))
    docs += sorted((root / ".claude/skills").glob("*/**/*.md"))
    docs += [root / p for p in sorted(instruction_files(root))]
    return docs


def path_references(root):
    """Backticked repository paths in harness documents, as (document, path) pairs."""
    refs = []
    for doc in harness_documents(root):
        if not doc.is_file():
            continue
        for token in re.findall(r"`([^`\s]+)`", doc.read_text()):
            token = re.sub(r":\d+(-\d+)?$", "", token)
            if token.startswith(PATH_PREFIXES) and not re.search(r"[<>*{}$]", token):
                refs.append((doc, token))
    return refs


def check_agents(root, skills, require, errors):
    agents_dir = root / ".claude/agents"
    require({p.stem for p in agents_dir.glob("*.md")} == set(AGENT_POLICY), "agent coverage must match AGENT_POLICY exactly")
    for role, (model, effort, preload) in sorted(AGENT_POLICY.items()):
        try:
            meta, body = frontmatter(agents_dir / f"{role}.md")
            require(meta.get("name") == role, f"agent name mismatch: {role}")
            require(bool(meta.get("description")), f"agent description missing: {role}")
            require(set(meta) == AGENT_KEYS, f"unexpected agent settings (memory/hooks/permissions): {role}")
            require({t.strip() for t in str(meta.get("tools", "")).split(",")} == AGENT_TOOLS, f"unsafe agent tools: {role}")
            require({t.strip() for t in str(meta.get("disallowedTools", "")).split(",")} == AGENT_DENY, f"agent tool exclusions missing: {role}")
            require(meta.get("model") == model and meta.get("effort") == effort, f"agent model policy mismatch: {role}")
            listed = meta.get("skills")
            require(isinstance(listed, list) and set(listed) == preload, f"agent preloaded skills mismatch: {role}")
            require(all(s in skills for s in listed or []), f"agent preloads a missing skill: {role}")
            require(".claude/INVARIANTS.md" in body, f"agent body must point to the charter: {role}")
        except (OSError, ValueError, HarnessError) as exc:
            errors.append(str(exc))


def check_settings(root, require, errors):
    try:
        settings = json.loads((root / ".claude/settings.json").read_text())
        require(settings.get("model") == "opus" and settings.get("effortLevel") == "high", "project model policy mismatch")
        denies = settings.get("permissions", {}).get("deny", [])
        require("Edit(/.claude/INVARIANTS.md)" in denies, "charter edit guard missing")
        require(not any(r.startswith(("Write(", "NotebookEdit(", "MultiEdit(")) for r in denies), "Write/NotebookEdit path rules are never consulted; use Edit(path)")
        missing = sorted(SECRET_READ_DENIES - set(denies))
        require(not missing, f"secret files must be read-denied (absolute forms): {missing}")
        hooks = settings.get("hooks", {}).get("PreToolUse", [])
        require(any("Agent" in h.get("matcher", "") and any(c.get("command") == HOOK_COMMAND for c in h.get("hooks", [])) for h in hooks), "model-pin hook not registered for Agent spawns with the fail-closed command")
        hook = load_hook(root)
        require(hook.OPUS_ONLY == {r for r, (m, _, _) in AGENT_POLICY.items() if m == "opus"}, "hook Opus pins disagree with AGENT_POLICY")
        require(hook.NO_HAIKU == {r for r, (m, _, _) in AGENT_POLICY.items() if m != "opus"}, "hook Haiku ban disagrees with AGENT_POLICY")
    except (OSError, ValueError, AttributeError, ImportError, SyntaxError) as exc:
        errors.append(f"settings/hook: {exc}")


def check_documents(root, require, errors):
    for doc in harness_documents(root):
        if not doc.is_file():
            errors.append(f"missing document: {doc.relative_to(root)}")
            continue
        for target in re.findall(r"\]\(([^)]+)\)", doc.read_text()):
            target = target.split(' "', 1)[0].strip("<>")
            if re.match(r"[a-zA-Z][\w+.-]*:", target) or target.startswith(("#", "/")):
                continue
            target = unquote(target.split("#", 1)[0])
            if target:
                require((doc.parent / target).exists(), f"broken link in {doc.relative_to(root)}: {target}")
    for doc, target in path_references(root):
        if (root / target).exists():
            continue
        excused = target.startswith(MACHINE_LOCAL) or target.rstrip("/") + "/" in MACHINE_LOCAL
        if excused:
            try:
                git(root, "check-ignore", "--no-index", "-q", target)
            except HarnessError:
                excused = False
        if not excused:
            errors.append(f"stale path in {doc.relative_to(root)}: {target}")


def routed_paths_present(root):
    """Lowercased tracked and untracked files plus their directories, as routing compares them."""
    files = {os.fsdecode(p).lower() for p in git(root, "ls-files", "--cached", "--others", "--exclude-standard", "-z").split(b"\0") if p}
    dirs = {str(parent).lower() + "/" for f in files for parent in Path(f).parents if str(parent) != "."}
    return files | dirs


def check(root):
    errors = []

    def require(condition, message):
        if not condition:
            errors.append(message)

    skills_dir = root / ".claude/skills"
    skills = {p.name for p in skills_dir.iterdir() if not p.name.startswith(".")} if skills_dir.is_dir() else set()
    check_agents(root, skills, require, errors)
    check_settings(root, require, errors)
    require(REQUIRED_SKILLS <= skills, f"required skills missing: {sorted(REQUIRED_SKILLS - skills)}")
    for name in sorted(skills):
        folder = skills_dir / name
        require(folder.is_dir() and not folder.is_symlink(), f"skill must be a real directory: {name}")
        try:
            meta, _ = frontmatter(folder / "SKILL.md")
            description = meta.get("description") or ""
            require(meta.get("name") == name and bool(description), f"invalid skill identity: {name}")
            require(len(description) <= MAX_SKILL_DESCRIPTION, f"skill description over {MAX_SKILL_DESCRIPTION} chars (always in context): {name}")
        except (OSError, HarnessError) as exc:
            errors.append(str(exc))
    # Instruction budgets: the root file loads into every session and every subagent.
    instructions = instruction_files(root)
    root_size = (root / "CLAUDE.md").stat().st_size if (root / "CLAUDE.md").is_file() else 0
    require(0 < root_size < ROOT_INSTRUCTIONS_LIMIT, f"root CLAUDE.md must be below {ROOT_INSTRUCTIONS_LIMIT} bytes: {root_size}")
    max_chain = 0
    for path in instructions:
        size = sum((root / q).stat().st_size for q in instructions if q.parent == path.parent or q.parent in path.parent.parents)
        max_chain = max(max_chain, size)
        require(size < CHAIN_LIMIT, f"instruction chain exceeds {CHAIN_LIMIT} bytes at {path}: {size}")
    try:
        git(root, "check-ignore", "--no-index", f"{TASK_DIR}probe.md")
    except HarnessError:
        errors.append("task state ignore rule missing")
    try:
        git(root, "check-ignore", "--no-index", "-q", ".claude/agents/probe.md")
        errors.append(".claude/ is ignored: the shared harness would never be committed")
    except HarnessError:
        pass
    # Routing compares lowercase paths; a renamed module must not silently drop out of its route.
    present = routed_paths_present(root)
    for routed in (*ANNOTATION_PATHS, *ANNOTATION_EXEMPT, *STORE_PATHS, *CONTRACT_PATHS, DIVERGENCES, ".github/workflows/release.yml"):
        require(routed in present, f"routed path no longer exists: {routed}")
    check_documents(root, require, errors)
    return {
        "ok": not errors,
        "errors": errors,
        "root_instruction_bytes": root_size,
        "max_instruction_chain_bytes": max_chain,
        "limits": "Static structure only; live agent loading, effective tools and review need separate verification.",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("check")
    snap = sub.add_parser("snapshot")
    snap.add_argument("--base", help="comparison whose start the --diff-out diff begins at (default HEAD)")
    snap.add_argument("--diff-out", type=Path, help="write the exact reviewed binary diff here")
    snap.add_argument("--allow", action="append", default=[], metavar="PATH", help="a confirmed-synthetic file that only looks sensitive (repeatable)")
    routing = sub.add_parser("route")
    routing.add_argument("--base", required=True, help="explicit revision, left..right, or left...right")
    routing.add_argument("--verbose", action="store_true", help="include per-file routing detail")
    args = parser.parse_args(argv)
    try:
        root = Path(os.fsdecode(git(args.root, "rev-parse", "--show-toplevel")).strip())
        if args.command == "check":
            result = check(root)
        elif args.command == "snapshot":
            result = snapshot(root, args.base, args.diff_out, args.allow)
        else:
            result = route(root, args.base, args.verbose)
        print(json.dumps(result, indent=2, ensure_ascii=True))
        return 0 if result.get("ok", True) else 1
    except (HarnessError, OSError, ValueError) as exc:
        print(f"harness error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
