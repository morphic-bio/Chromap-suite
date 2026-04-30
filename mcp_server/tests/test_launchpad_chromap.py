from pathlib import Path


def test_launchpad_defaults_to_chromap_recipe_filter():
    app_js = (Path(__file__).resolve().parents[1] / "launchpad" / "static" / "app.js").read_text()
    assert "CHROMAP_WORKFLOW_ORDER" in app_js
    assert 'startsWith("chromap_")' in app_js
    assert ("star" + "_genome_generate") not in app_js


def test_launchpad_page_uses_chromap_labels():
    index = (Path(__file__).resolve().parents[1] / "launchpad" / "static" / "index.html").read_text()
    assert "Chromap Launchpad" in index
    assert "basic Chromap CLI" in index
