import json
from pathlib import Path

from starlette.testclient import TestClient

from mcp_server.app import build_http_app, write_recipe_run_manifest as manifest_tool
from mcp_server.tools.run_manifest import write_run_manifest


def _touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x\n")
    return str(path)


def _params(temp_dir: Path) -> dict:
    return {
        "reference": _touch(temp_dir / "ref.fa"),
        "index": _touch(temp_dir / "ref.index"),
        "read1": str(temp_dir / "r1.fastq.gz"),
        "read2": str(temp_dir / "r2.fastq.gz"),
        "output": str(temp_dir / "out.bed"),
    }


def test_dry_run_manifest_records_required_fields(loaded_default_config, temp_dir):
    loaded_default_config.paths.artifact_log_root = temp_dir / "artifacts"

    result = write_run_manifest("chromap_atac_bed", _params(temp_dir))
    manifest_path = Path(result["manifest_path"])
    manifest = json.loads(manifest_path.read_text())

    assert manifest_path.name == "run.json"
    assert manifest["manifest_version"] == 1
    assert manifest["recipe_id"] == "chromap_atac_bed"
    assert manifest["registry_version"] == 1
    assert manifest["execution_status"] == "dry_run"
    assert manifest["argv"][0].endswith("chromap")
    assert "--preset" in manifest["argv"]
    assert manifest["git"]["commit"]
    assert isinstance(manifest["git"]["dirty"], bool)
    assert manifest["preflight"]["valid"] is True
    assert manifest["output_paths"] == [str(temp_dir / "out.bed")]
    assert Path(manifest["log_paths"]["stdout"]).exists()
    assert Path(manifest["log_paths"]["stderr"]).exists()


def test_manifest_tool_writes_dry_run_manifest(loaded_default_config, temp_dir):
    loaded_default_config.paths.artifact_log_root = temp_dir / "artifacts"

    result = manifest_tool.fn("chromap_atac_bed", _params(temp_dir))

    assert "manifest_path" in result
    assert result["manifest"]["execution_status"] == "dry_run"
    assert result["manifest"]["recipe_id"] == "chromap_atac_bed"


def test_run_manifest_records_serial_benchmark_policy(loaded_default_config, temp_dir):
    loaded_default_config.paths.artifact_log_root = temp_dir / "artifacts"

    result = write_run_manifest(
        "chromap_macs3_frag_peaks",
        {
            "macs3_frag_peaks_output": str(temp_dir / "peaks.narrowPeak"),
            "macs3_frag_summits_output": str(temp_dir / "summits.bed"),
            "macs3_pvalue": "0.01",
        },
        execution_status="dry_run",
    )
    manifest = result["manifest"]

    assert manifest["benchmark_policy"] == "serial_required"
    assert manifest["classification"]["benchmark_policy"] == "serial_required"


def test_launchpad_launch_writes_manifest_and_logs(
    loaded_default_config,
    temp_dir,
    monkeypatch,
):
    loaded_default_config.paths.artifact_log_root = temp_dir / "artifacts"
    params = _params(temp_dir)

    class FakePopen:
        def __init__(self, argv, **kwargs):
            self.argv = argv
            self.kwargs = kwargs
            self.pid = 4242
            assert kwargs["stdout"] is not None
            assert kwargs["stderr"] is not None

    monkeypatch.setattr("mcp_server.launchpad.api.subprocess.Popen", FakePopen)

    with TestClient(build_http_app()) as client:
        response = client.post(
            "/launchpad/api/workflows/chromap_atac_bed/launch",
            json={"params": params},
        )

    assert response.status_code == 200
    payload = response.json()
    assert payload["pid"] == 4242
    manifest_path = Path(payload["manifest_path"])
    manifest = json.loads(manifest_path.read_text())
    assert manifest["execution_status"] == "running"
    assert manifest["pid"] == 4242
    assert Path(payload["stdout_log"]).exists()
    assert Path(payload["stderr_log"]).exists()
