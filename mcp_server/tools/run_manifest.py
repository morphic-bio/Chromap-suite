"""Run manifest helpers for recipe-driven Chromap executions."""

from __future__ import annotations

import json
import os
import platform
import shlex
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from string import Formatter
from typing import Any
from uuid import uuid4

from ..config import get_config
from ..schemas.recipe import RecipeEntry
from .preflight import preflight_recipe
from .recipes import get_recipe, load_recipe_registry
from .workflows import render_workflow_command


_ENV_PREFIXES = ("CHROMAP", "MCP", "OMP", "SLURM", "STAR")


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _utc_stamp(ts: datetime | None = None) -> str:
    return (ts or _utc_now()).strftime("%Y%m%dT%H%M%SZ")


def _iso(ts: datetime | None) -> str | None:
    return ts.isoformat().replace("+00:00", "Z") if ts else None


def _repo_root() -> Path:
    return Path(get_config().paths.repo_root).expanduser().resolve()


def _run_git(args: list[str]) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=str(_repo_root()),
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def _git_state() -> dict[str, Any]:
    status = _run_git(["status", "--porcelain"])
    return {
        "commit": _run_git(["rev-parse", "HEAD"]),
        "dirty": bool(status),
    }


def _safe_format(template: str, params: dict[str, Any]) -> str:
    values: dict[str, Any] = {}
    for _, field_name, _, _ in Formatter().parse(template):
        if field_name and field_name not in values:
            values[field_name] = params.get(field_name, "{" + field_name + "}")
    try:
        return template.format(**values)
    except Exception:
        return template


def _render_template_list(template: list[str], params: dict[str, Any]) -> list[str]:
    return [_safe_format(str(token), params) for token in template]


def _recipe_for_workflow(workflow_id: str) -> RecipeEntry | None:
    for recipe in load_recipe_registry().recipes:
        if recipe.workflow_id == workflow_id:
            return recipe
    return None


def recipe_id_for_workflow(workflow_id: str) -> str | None:
    recipe = _recipe_for_workflow(workflow_id)
    return recipe.id if recipe else None


def _output_paths(recipe: RecipeEntry, params: dict[str, Any]) -> list[str]:
    paths: list[str] = []
    for output in recipe.outputs:
        rendered = _safe_format(output.path_template, params).strip()
        if rendered:
            paths.append(rendered)
    return paths


def _input_paths(recipe: RecipeEntry, params: dict[str, Any]) -> dict[str, str]:
    out: dict[str, str] = {}
    for input_def in recipe.inputs:
        if input_def.name not in params:
            continue
        value = params[input_def.name]
        if value in (None, ""):
            continue
        if input_def.type in ("file", "directory"):
            out[input_def.name] = str(value)
    return out


def _reference_paths(params: dict[str, Any]) -> dict[str, str]:
    keys = ("reference", "ref", "genome", "index")
    return {key: str(params[key]) for key in keys if params.get(key)}


def _relevant_environment(env_overrides: dict[str, str] | None = None) -> dict[str, str]:
    env: dict[str, str] = {}
    for key, value in os.environ.items():
        if key in ("PATH", "LD_LIBRARY_PATH") or key.startswith(_ENV_PREFIXES):
            env[key] = value
    for key, value in (env_overrides or {}).items():
        env[str(key)] = str(value)
    return env


def _binary_version(binary: str) -> str | None:
    candidates = ([binary, "--version"], [binary, "-v"])
    for cmd in candidates:
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=5,
            )
        except Exception:
            continue
        text = (result.stdout or result.stderr).strip()
        if text:
            return text.splitlines()[0]
    return None


def _binary_info(argv: list[str]) -> list[dict[str, Any]]:
    if not argv:
        return []
    binary = argv[0]
    return [
        {
            "name": Path(binary).name,
            "path": binary,
            "version": _binary_version(binary),
        }
    ]


def _next_manifest_dir(recipe_id: str, started_at: datetime) -> Path:
    root = Path(get_config().paths.artifact_log_root).expanduser().resolve() / "mcp_runs"
    base = root / f"{_utc_stamp(started_at)}_{recipe_id}"
    candidate = base
    while candidate.exists():
        candidate = root / f"{_utc_stamp(started_at)}_{recipe_id}_{uuid4().hex[:8]}"
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def _render_recipe_command(
    recipe: RecipeEntry,
    params: dict[str, Any],
    argv: list[str] | None,
    shell_preview: str | None,
    env_overrides: dict[str, str] | None,
) -> tuple[list[str], str, dict[str, str]]:
    if argv is not None:
        final_argv = list(argv)
        final_shell = shell_preview or " ".join(shlex.quote(a) for a in final_argv)
        return final_argv, final_shell, dict(env_overrides or {})

    if recipe.workflow_id:
        rendered = render_workflow_command(recipe.workflow_id, params)
        return (
            list(rendered.argv or []),
            rendered.shell_preview,
            dict(rendered.env_overrides or {}),
        )

    final_argv = _render_template_list(recipe.command_template, params)
    return final_argv, " ".join(shlex.quote(a) for a in final_argv), dict(env_overrides or {})


def build_run_manifest(
    recipe_id: str,
    params: dict[str, Any],
    *,
    execution_status: str = "dry_run",
    argv: list[str] | None = None,
    shell_preview: str | None = None,
    env_overrides: dict[str, str] | None = None,
    artifact_dir: Path | None = None,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
    pid: int | None = None,
    exit_code: int | None = None,
    error: str | None = None,
    started_at: datetime | None = None,
    ended_at: datetime | None = None,
) -> dict[str, Any]:
    """Build a complete run manifest dictionary without writing it."""
    recipe = get_recipe(recipe_id)
    if recipe is None:
        raise ValueError(f"Unknown recipe: {recipe_id}")

    started_at = started_at or _utc_now()
    registry = load_recipe_registry()
    final_argv, final_shell, final_env = _render_recipe_command(
        recipe, params, argv, shell_preview, env_overrides
    )
    preflight = preflight_recipe(recipe_id, params)
    artifact_dir_str = str(artifact_dir) if artifact_dir else None
    log_paths = {
        "stdout": str(stdout_path) if stdout_path else None,
        "stderr": str(stderr_path) if stderr_path else None,
    }
    manifest = {
        "manifest_version": 1,
        "recipe_id": recipe.id,
        "registry_version": registry.registry_version,
        "workflow_id": recipe.workflow_id,
        "runtime_class": recipe.runtime_class,
        "benchmark_policy": recipe.benchmark_policy,
        "benchmark_execution_policy": (
            "serial_only" if recipe.runtime_class == "benchmark" else None
        ),
        "execution_status": execution_status,
        "argv": final_argv,
        "shell_preview": final_shell,
        "git": _git_state(),
        "binaries": _binary_info(final_argv),
        "input_paths": _input_paths(recipe, params),
        "reference_paths": _reference_paths(params),
        "output_paths": _output_paths(recipe, params),
        "artifact_dir": artifact_dir_str,
        "log_paths": log_paths,
        "started_at": _iso(started_at),
        "ended_at": _iso(ended_at),
        "exit_code": exit_code,
        "pid": pid,
        "host": socket.gethostname(),
        "user": os.environ.get("USER") or os.environ.get("LOGNAME"),
        "platform": platform.platform(),
        "environment": _relevant_environment(final_env),
        "preflight": preflight.model_dump(),
        "classification": {
            "runtime_class": recipe.runtime_class,
            "benchmark_policy": recipe.benchmark_policy,
            "smoke_coverage": recipe.smoke_coverage,
        },
        "params": dict(params),
    }
    if error:
        manifest["error"] = error
    return manifest


def write_run_manifest(
    recipe_id: str,
    params: dict[str, Any],
    *,
    execution_status: str = "dry_run",
    argv: list[str] | None = None,
    shell_preview: str | None = None,
    env_overrides: dict[str, str] | None = None,
    artifact_dir: Path | None = None,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
    pid: int | None = None,
    exit_code: int | None = None,
    error: str | None = None,
    started_at: datetime | None = None,
    ended_at: datetime | None = None,
) -> dict[str, Any]:
    """Write a recipe run manifest and return path plus manifest payload."""
    started_at = started_at or _utc_now()
    artifact_dir = artifact_dir or _next_manifest_dir(recipe_id, started_at)
    stdout_path = stdout_path or artifact_dir / "stdout.log"
    stderr_path = stderr_path or artifact_dir / "stderr.log"
    stdout_path.touch(exist_ok=True)
    stderr_path.touch(exist_ok=True)
    manifest = build_run_manifest(
        recipe_id,
        params,
        execution_status=execution_status,
        argv=argv,
        shell_preview=shell_preview,
        env_overrides=env_overrides,
        artifact_dir=artifact_dir,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        pid=pid,
        exit_code=exit_code,
        error=error,
        started_at=started_at,
        ended_at=ended_at,
    )
    manifest_path = artifact_dir / "run.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {
        "manifest_path": str(manifest_path),
        "artifact_dir": str(artifact_dir),
        "manifest": manifest,
    }


def update_run_manifest(
    manifest_path: str | Path,
    **updates: Any,
) -> dict[str, Any]:
    """Update fields in an existing run manifest."""
    path = Path(manifest_path)
    manifest = json.loads(path.read_text())
    manifest.update(updates)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


__all__ = [
    "build_run_manifest",
    "recipe_id_for_workflow",
    "update_run_manifest",
    "write_run_manifest",
]
