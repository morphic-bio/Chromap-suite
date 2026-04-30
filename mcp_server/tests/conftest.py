"""Pytest fixtures for the Chromap MCP server."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Generator

import pytest

import mcp_server.config as config_module
from mcp_server.config import load_config


@pytest.fixture(autouse=True)
def reset_config_state() -> Generator[None, None, None]:
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    with tempfile.TemporaryDirectory() as tmp:
        yield Path(tmp)


@pytest.fixture
def loaded_default_config():
    return load_config(Path(__file__).resolve().parents[1] / "config.yaml")
