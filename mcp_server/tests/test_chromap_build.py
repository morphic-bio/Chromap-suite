from pathlib import Path

from mcp_server.tools.build import needs_rebuild


def test_unknown_build_target_is_reported(temp_dir):
    needed, reason = needs_rebuild(temp_dir, "not-a-target")
    assert needed is True
    assert "Unknown target" in reason


def test_missing_chromap_output_needs_rebuild(temp_dir):
    (temp_dir / "src").mkdir()
    (temp_dir / "src" / "chromap.cc").write_text("int main() { return 0; }\n")
    (temp_dir / "Makefile").write_text("all:\n\t@true\n")

    needed, reason = needs_rebuild(Path(temp_dir), "chromap")
    assert needed is True
    assert "Build output missing" in reason
