from pathlib import Path

from mcp_server.config import get_workflow_schemas
from mcp_server.tools.workflows import (
    get_workflow_parameter_schema,
    list_workflows,
    render_workflow_command,
    validate_workflow_parameters,
)


EXPECTED_WORKFLOWS = {
    "chromap_index",
    "chromap_atac_bed",
    "chromap_atac_bam_fragments",
    "chromap_hic_pairs",
}


def test_default_config_loads_chromap_workflows(loaded_default_config):
    schemas = get_workflow_schemas()
    assert set(schemas) == EXPECTED_WORKFLOWS
    listed = {w.id for w in list_workflows(authenticated=False).workflows}
    assert listed == EXPECTED_WORKFLOWS


def test_chromap_index_renders_build_command(loaded_default_config):
    result = render_workflow_command(
        "chromap_index",
        {
            "reference": "ref.fa",
            "index": "ref.chromap.idx",
            "num_threads": 16,
        },
    )
    assert Path(result.argv[0]).name == "chromap"
    assert result.argv == [
        result.argv[0],
        "--build-index",
        "--ref",
        "ref.fa",
        "--output",
        "ref.chromap.idx",
        "--num-threads",
        "16",
    ]


def test_atac_bam_fragments_renders_peak_options(loaded_default_config):
    result = render_workflow_command(
        "chromap_atac_bam_fragments",
        {
            "reference": "ref.fa",
            "index": "ref.idx",
            "read1": "r1.fastq.gz",
            "read2": "r2.fastq.gz",
            "barcode": "bc.fastq.gz",
            "output": "out.bam",
            "atac_fragments": "fragments.tsv.gz",
            "call_macs3_frag_peaks": True,
            "macs3_frag_peaks_output": "peaks.narrowPeak",
            "macs3_frag_summits_output": "summits.bed",
            "macs3_frag_peaks_source": "memory",
        },
    )
    assert "--BAM" in result.argv
    assert "--atac-fragments" in result.argv
    assert "--call-macs3-frag-peaks" in result.argv
    assert result.argv[result.argv.index("--macs3-frag-peaks-source") + 1] == "memory"


def test_hic_pairs_renders_pairs_command(loaded_default_config):
    result = render_workflow_command(
        "chromap_hic_pairs",
        {
            "reference": "ref.fa",
            "index": "ref.idx",
            "read1": "hic_R1.fastq.gz",
            "read2": "hic_R2.fastq.gz",
            "output": "out.pairs.gz",
        },
    )
    assert "--preset" in result.argv
    assert result.argv[result.argv.index("--preset") + 1] == "hic"
    assert "--pairs" in result.argv


def test_schema_validation_reports_required_params(loaded_default_config):
    validation = validate_workflow_parameters("chromap_atac_bed", {}, check_paths=False)
    assert not validation.valid
    assert any("reference" in error for error in validation.errors)


def test_parameter_schema_has_groups(loaded_default_config):
    schema = get_workflow_parameter_schema("chromap_atac_bed")
    groups = {g.name for g in schema.parameter_groups}
    assert {"inputs", "output", "mapping"} <= groups
