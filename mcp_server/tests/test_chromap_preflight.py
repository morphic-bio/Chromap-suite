from pathlib import Path

from mcp_server.app import preflight_recipe as preflight_recipe_tool
from mcp_server.tools.preflight import preflight_recipe


def _touch(path: Path, content: str = "x\n") -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return str(path)


def _base_params(temp_dir: Path) -> dict:
    return {
        "reference": _touch(temp_dir / "ref.fa", ">chr1\nACGT\n"),
        "index": _touch(temp_dir / "ref.index"),
        "read1": _touch(temp_dir / "r1.fastq.gz"),
        "read2": _touch(temp_dir / "r2.fastq.gz"),
        "output": str(temp_dir / "out.bed"),
    }


def test_recipe_preflight_passes_for_atac_bed(loaded_default_config, temp_dir):
    result = preflight_recipe("chromap_atac_bed", _base_params(temp_dir))

    assert result.valid
    assert not result.errors
    statuses = {check.rule_id: check.status for check in result.checks}
    assert statuses["reference_fasta_exists"] == "pass"
    assert statuses["chromap_index_exists"] == "pass"
    assert statuses["input_paths_trusted"] == "pass"
    assert statuses["read_pair_lists_match"] == "pass"
    assert statuses["output_parent_trusted_or_creatable"] == "pass"


def test_recipe_preflight_reports_missing_index(loaded_default_config, temp_dir):
    params = _base_params(temp_dir)
    params["index"] = str(temp_dir / "missing.index")

    result = preflight_recipe("chromap_atac_bed", params)

    assert not result.valid
    assert any(check.rule_id == "chromap_index_exists" for check in result.checks)
    assert any("Path does not exist" in error for error in result.errors)


def test_recipe_preflight_reports_mismatched_read_lists(loaded_default_config, temp_dir):
    params = _base_params(temp_dir)
    params["read1"] = f"{temp_dir / 'r1a.fq.gz'},{temp_dir / 'r1b.fq.gz'}"
    params["read2"] = str(temp_dir / "r2a.fq.gz")

    result = preflight_recipe("chromap_atac_bed", params)

    assert not result.valid
    mismatch = [c for c in result.checks if c.rule_id == "read_pair_lists_match"][0]
    assert mismatch.status == "fail"
    assert "lane counts differ" in mismatch.message


def test_recipe_preflight_reports_output_collision(loaded_default_config, temp_dir):
    params = {
        **_base_params(temp_dir),
        "barcode": _touch(temp_dir / "bc.fastq.gz"),
        "output": str(temp_dir / "out.bam"),
        "atac_fragments": str(temp_dir / "out.bam"),
    }

    result = preflight_recipe("chromap_atac_bam_fragments", params)

    assert not result.valid
    collision = [c for c in result.checks if c.rule_id == "secondary_output_not_primary"][0]
    assert collision.status == "fail"


def test_recipe_preflight_reports_barcode_lane_mismatch(loaded_default_config, temp_dir):
    params = {
        **_base_params(temp_dir),
        "barcode": f"{temp_dir / 'bc1.fq.gz'},{temp_dir / 'bc2.fq.gz'}",
        "output": str(temp_dir / "out.bam"),
        "atac_fragments": str(temp_dir / "fragments.tsv.gz"),
    }

    result = preflight_recipe("chromap_atac_bam_fragments", params)

    assert not result.valid
    barcode = [c for c in result.checks if c.rule_id == "barcode_list_matches_read1"][0]
    assert barcode.status == "fail"


def test_recipe_preflight_rejects_untrusted_output(loaded_default_config, temp_dir):
    params = _base_params(temp_dir)
    params["output"] = "/root/chromap/out.bed"

    result = preflight_recipe("chromap_atac_bed", params)

    assert not result.valid
    output = [c for c in result.checks if c.rule_id == "output_parent_trusted_or_creatable"][0]
    assert output.status == "fail"
    assert "trusted roots" in output.message


def test_recipe_preflight_rejects_untrusted_input(loaded_default_config, temp_dir):
    params = _base_params(temp_dir)
    params["read1"] = "/etc/passwd"

    result = preflight_recipe("chromap_atac_bed", params)

    assert not result.valid
    input_paths = [c for c in result.checks if c.rule_id == "input_paths_trusted"][0]
    assert input_paths.status == "fail"
    assert "outside trusted roots" in input_paths.message


def test_recipe_preflight_rejects_shell_metacharacters(loaded_default_config, temp_dir):
    params = _base_params(temp_dir)
    params["output"] = str(temp_dir / "out.bed;touch BAD")

    result = preflight_recipe("chromap_atac_bed", params)

    assert not result.valid
    output = [c for c in result.checks if c.rule_id == "output_parent_trusted_or_creatable"][0]
    assert output.status == "fail"
    assert "shell metacharacters" in output.message


def test_recipe_preflight_enforces_hic_pairs_suffix(loaded_default_config, temp_dir):
    params = _base_params(temp_dir)
    params["output"] = str(temp_dir / "hic.cool")

    result = preflight_recipe("chromap_hic_pairs", params)

    assert not result.valid
    hic = [c for c in result.checks if c.rule_id == "hic_stops_at_pairs"][0]
    assert hic.status == "fail"
    assert ".pairs" in hic.message


def test_recipe_preflight_mcp_tool(loaded_default_config, temp_dir):
    result = preflight_recipe_tool.fn("chromap_atac_bed", _base_params(temp_dir))

    assert result["valid"] is True
    assert result["recipe_id"] == "chromap_atac_bed"
    assert any(check["rule_id"] == "chromap_index_exists" for check in result["checks"])
