from pathlib import Path

from starlette.testclient import TestClient

from mcp_server.app import build_http_app


def test_launchpad_defaults_to_chromap_recipe_filter():
    app_js = (Path(__file__).resolve().parents[1] / "launchpad" / "static" / "app.js").read_text()
    assert "sortRecipesForLaunchpad" in app_js
    assert "/recipes" in app_js
    assert ("star" + "_genome_generate") not in app_js


def test_launchpad_page_uses_chromap_labels():
    index = (Path(__file__).resolve().parents[1] / "launchpad" / "static" / "index.html").read_text()
    assert "Chromap Launchpad" in index
    assert "basic Chromap CLI" in index
    assert "Expected outputs" in index
    assert "Write dry-run manifest" in index


def _touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x\n")
    return str(path)


def _recipe_params(temp_dir: Path) -> dict:
    return {
        "reference": _touch(temp_dir / "ref.fa"),
        "index": _touch(temp_dir / "ref.index"),
        "read1": str(temp_dir / "r1.fastq.gz"),
        "read2": str(temp_dir / "r2.fastq.gz"),
        "output": str(temp_dir / "out.bed"),
    }


def test_launchpad_recipe_api_lists_enabled_recipes_only(loaded_default_config):
    with TestClient(build_http_app()) as client:
        response = client.get("/launchpad/api/recipes")

    assert response.status_code == 200
    payload = response.json()
    ids = {recipe["id"] for recipe in payload["recipes"]}
    assert "chromap_atac_bed" in ids
    assert "chromap_macs3_frag_peaks" not in ids


def test_launchpad_recipe_api_can_list_metadata_only_but_not_execute(loaded_default_config):
    with TestClient(build_http_app()) as client:
        listed = client.get("/launchpad/api/recipes?include_disabled=1")
        schema = client.get("/launchpad/api/recipes/chromap_macs3_frag_peaks/schema")

    assert listed.status_code == 200
    ids = {recipe["id"] for recipe in listed.json()["recipes"]}
    assert "chromap_macs3_frag_peaks" in ids
    assert schema.status_code == 400
    assert schema.json()["code"] == "NOT_EXECUTABLE"


def test_launchpad_recipe_schema_is_registry_driven(loaded_default_config):
    with TestClient(build_http_app()) as client:
        response = client.get("/launchpad/api/recipes/chromap_atac_bed/schema")

    assert response.status_code == 200
    payload = response.json()
    names = [p["name"] for p in payload["parameters"]]
    assert payload["recipe_id"] == "chromap_atac_bed"
    assert payload["workflow_id"] == "chromap_atac_bed"
    assert names == ["reference", "index", "read1", "read2", "output"]
    assert payload["outputs"][0]["artifact_type"] == "bed"


def test_launchpad_recipe_preflight_render_and_manifest(loaded_default_config, temp_dir):
    loaded_default_config.paths.artifact_log_root = temp_dir / "artifacts"
    params = _recipe_params(temp_dir)

    with TestClient(build_http_app()) as client:
        preflight = client.post(
            "/launchpad/api/recipes/chromap_atac_bed/preflight",
            json={"params": params},
        )
        render = client.post(
            "/launchpad/api/recipes/chromap_atac_bed/render",
            json={"params": params},
        )
        manifest = client.post(
            "/launchpad/api/recipes/chromap_atac_bed/manifest",
            json={"params": params},
        )

    assert preflight.status_code == 200
    assert preflight.json()["valid"] is True
    assert render.status_code == 200
    assert render.json()["recipe_id"] == "chromap_atac_bed"
    assert "--preset" in render.json()["argv"]
    assert manifest.status_code == 200
    assert manifest.json()["manifest"]["execution_status"] == "dry_run"
