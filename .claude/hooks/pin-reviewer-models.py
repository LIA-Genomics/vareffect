#!/usr/bin/env python3
"""PreToolUse hook: block spawning a vareffect reviewer on a model below its policy floor.

A per-invocation `model` overrides an agent's frontmatter, and environment settings can
override or remap both, so the pins in .claude/agents/ are only enforceable here. Exit 2
blocks the call and shows the reason to Claude. Anything that is not a vareffect reviewer spawn
passes through; an unreadable spawn request fails closed.
"""

import json
import os
import re
import sys

OPUS_ONLY = {"vareffect-bioinformatics"}
NO_HAIKU = {"vareffect-reviewer", "vareffect-researcher"}
# Exact aliases, or a full Opus/Fable model id including Bedrock (`us.anthropic.claude-opus-…-v1:0`,
# inference-profile ARNs) and Vertex (`claude-opus-…@…`) forms. No prefix match: `opusplan` runs Sonnet.
OPUS_ALIASES = {"opus", "fable", "opus[1m]", "fable[1m]"}
OPUS_ID = re.compile(r"(^|[/.:])claude-(opus|fable)-[0-9][0-9a-z.:@-]*(\[1m\])?$")


def opus_class(model):
    model = model.strip().lower()
    return model in OPUS_ALIASES or bool(OPUS_ID.search(model))


def environment_downgrade(env):
    """Return why the environment would run a pinned specialist below Opus, else None."""
    forced = env.get("CLAUDE_CODE_SUBAGENT_MODEL_FORCE", "").strip().lower() not in {"", "0", "false"}
    subagent_model = env.get("CLAUDE_CODE_SUBAGENT_MODEL", "").strip()
    if forced and not (subagent_model and opus_class(subagent_model)):
        return "CLAUDE_CODE_SUBAGENT_MODEL_FORCE is set without an Opus CLAUDE_CODE_SUBAGENT_MODEL (set that to an Opus model id, or unset the force flag)"
    remap = env.get("ANTHROPIC_DEFAULT_OPUS_MODEL", "").strip()
    if remap and not opus_class(remap):
        return f"ANTHROPIC_DEFAULT_OPUS_MODEL maps the opus alias to {remap} (point it at an Opus model id for this provider)"
    return None


def verdict(event, env=None):
    """Return a block reason for a disallowed reviewer spawn, else None."""
    env = os.environ if env is None else env
    if not isinstance(event, dict):
        return "unreadable Agent call; vareffect model pins cannot be checked"
    if event.get("tool_name") not in {"Agent", "Task"}:
        return None
    tool_input = event.get("tool_input")
    if not isinstance(tool_input, dict):
        return "unreadable Agent call; vareffect model pins cannot be checked"
    role = tool_input.get("subagent_type")
    role = role.strip().lower() if isinstance(role, str) else ""
    if role not in OPUS_ONLY | NO_HAIKU:
        return None
    model = tool_input.get("model")
    if model is not None and not isinstance(model, str):
        return f"{role}: unreadable model parameter; vareffect model pins cannot be checked"
    model = (model or "").strip().lower()
    if role in OPUS_ONLY:
        if model and not opus_class(model):
            return f"{role} is pinned to Opus by vareffect policy; omit the model parameter (requested: {model})."
        downgrade = environment_downgrade(env)
        if downgrade:
            return f"{role} is pinned to Opus by vareffect policy, but {downgrade}."
    elif "haiku" in model:
        return f"{role} may not run on Haiku by vareffect policy; omit the model or choose sonnet/opus."
    return None


def main():
    try:
        event = json.load(sys.stdin)
    except (ValueError, OSError):
        print("vareffect model-pin hook could not read the Agent call; blocking to keep reviewer pins enforced.", file=sys.stderr)
        return 2
    try:
        reason = verdict(event)
    except Exception as exc:  # fail closed: a crashed guard must not let a pinned spawn through
        reason = f"vareffect model-pin hook failed ({type(exc).__name__}); blocking."
    if reason:
        print(reason, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
