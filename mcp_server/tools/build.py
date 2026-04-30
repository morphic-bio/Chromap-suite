"""Build helpers for the Chromap-suite MCP server."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Any

BUILD_STATE_FILE = ".chromap_build_state.json"

TARGETS: dict[str, dict[str, Any]] = {
    "chromap": {
        "make_target": None,
        "outputs": ["chromap", "chromap_callpeaks"],
        "source_dirs": ["src", "third_party/libMACS3/include"],
    },
    "libchromap": {
        "make_target": "libchromap.a",
        "outputs": ["libchromap.a"],
        "source_dirs": ["src", "third_party/libMACS3/include"],
    },
    "runner": {
        "make_target": "chromap_lib_runner",
        "outputs": ["chromap_lib_runner"],
        "source_dirs": ["src", "third_party/libMACS3/include"],
    },
}


def get_source_hash(source_dir: Path, extensions: list[str] | None = None) -> str:
    """Compute a lightweight source tree hash from paths, mtimes, and sizes."""
    if extensions is None:
        extensions = [".cc", ".cpp", ".c", ".h", ".hh", ".hpp"]

    hasher = hashlib.sha256()
    for ext in extensions:
        for src_file in sorted(source_dir.rglob(f"*{ext}")):
            try:
                stat = src_file.stat()
            except OSError:
                continue
            hasher.update(f"{src_file}:{stat.st_mtime}:{stat.st_size}".encode())
    return hasher.hexdigest()[:16]


def load_build_state(repo_root: Path) -> dict[str, Any]:
    """Load the persisted build state."""
    state_file = repo_root / BUILD_STATE_FILE
    if state_file.exists():
        try:
            return json.loads(state_file.read_text())
        except (OSError, json.JSONDecodeError):
            pass
    return {}


def save_build_state(repo_root: Path, state: dict[str, Any]) -> None:
    """Persist build state under the repository root."""
    (repo_root / BUILD_STATE_FILE).write_text(json.dumps(state, indent=2) + "\n")


def _target_config(target: str) -> dict[str, Any] | None:
    return TARGETS.get(target)


def _target_source_hash(repo_root: Path, target: str) -> str:
    cfg = _target_config(target)
    if not cfg:
        return ""
    pieces: list[str] = []
    for rel in cfg["source_dirs"]:
        src_dir = repo_root / rel
        if src_dir.exists():
            pieces.append(get_source_hash(src_dir))
    makefile = repo_root / "Makefile"
    if makefile.exists():
        stat = makefile.stat()
        pieces.append(f"Makefile:{stat.st_mtime}:{stat.st_size}")
    return hashlib.sha256("|".join(pieces).encode()).hexdigest()[:16]


def needs_rebuild(repo_root: Path, target: str) -> tuple[bool, str]:
    """Return whether *target* needs a rebuild and why."""
    cfg = _target_config(target)
    if not cfg:
        return True, f"Unknown target '{target}'. Valid targets: {', '.join(TARGETS)}"

    missing = [p for p in cfg["outputs"] if not (repo_root / p).exists()]
    if missing:
        return True, f"Build output missing: {', '.join(missing)}"

    state = load_build_state(repo_root)
    current_hash = _target_source_hash(repo_root, target)
    saved_hash = state.get(target, {}).get("source_hash", "")
    if saved_hash != current_hash:
        return True, "Source files changed since last build"

    return False, "Build is up to date"


def build_target(
    repo_root: Path,
    target: str,
    clean: bool = False,
    force: bool = False,
    parallel: int | None = None,
) -> dict[str, Any]:
    """Build a Chromap-suite target with stale-build tracking."""
    result: dict[str, Any] = {
        "target": target,
        "success": False,
        "clean": clean,
        "forced": force,
        "start_time": datetime.now().isoformat(),
        "duration_seconds": 0,
        "output": "",
        "error": "",
    }
    start = time.time()

    cfg = _target_config(target)
    if not cfg:
        result["error"] = f"Unknown target: {target}. Valid targets: {', '.join(TARGETS)}"
        result["duration_seconds"] = time.time() - start
        return result

    if not force and not clean:
        needs_build, reason = needs_rebuild(repo_root, target)
        if not needs_build:
            result.update(
                {
                    "success": True,
                    "skipped": True,
                    "reason": reason,
                    "duration_seconds": time.time() - start,
                }
            )
            return result

    if parallel is None:
        parallel = os.cpu_count() or 4

    try:
        if clean:
            proc = subprocess.run(
                ["make", "clean"],
                cwd=repo_root,
                capture_output=True,
                text=True,
                timeout=300,
            )
            result["clean_output"] = proc.stdout + proc.stderr
            if proc.returncode != 0:
                result["error"] = f"make clean failed: {proc.stderr}"
                result["returncode"] = proc.returncode
                result["duration_seconds"] = time.time() - start
                return result

        cmd = ["make", f"-j{parallel}"]
        if cfg["make_target"]:
            cmd.append(str(cfg["make_target"]))

        proc = subprocess.run(
            cmd,
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=1800,
        )
        result["output"] = proc.stdout
        result["error"] = proc.stderr
        result["returncode"] = proc.returncode

        if proc.returncode == 0:
            result["success"] = True
            state = load_build_state(repo_root)
            state[target] = {
                "source_hash": _target_source_hash(repo_root, target),
                "build_time": datetime.now().isoformat(),
                "clean_build": clean,
                "make_target": cfg["make_target"] or "all",
            }
            save_build_state(repo_root, state)
    except subprocess.TimeoutExpired:
        result["error"] = "Build timed out after 30 minutes"
    except Exception as exc:  # pragma: no cover - defensive wrapper for MCP callers.
        result["error"] = str(exc)

    result["duration_seconds"] = time.time() - start
    return result


def ensure_clean_build(repo_root: Path, target: str = "chromap") -> dict[str, Any]:
    """Force a clean build for *target*."""
    return build_target(repo_root, target, clean=True, force=True)


def smart_build(repo_root: Path, target: str = "chromap") -> dict[str, Any]:
    """Build *target* only when source hashes or outputs indicate staleness."""
    return build_target(repo_root, target, clean=False, force=False)
