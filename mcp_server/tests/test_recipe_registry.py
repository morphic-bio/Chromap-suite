from mcp_server.config import get_workflow_schemas
from mcp_server.app import describe_recipe, list_recipes as list_recipes_tool
from mcp_server.tools.recipes import get_recipe, list_recipes, load_recipe_registry


PUBLIC_WORKFLOWS = {
    "chromap_index",
    "chromap_atac_bed",
    "chromap_atac_bam_fragments",
    "chromap_hic_pairs",
}

PLANNED_RECIPES = {
    "chromap_chip_tagalign",
    "chromap_sorted_bam",
    "chromap_y_noy_split",
    "chromap_macs3_frag_peaks",
    "chromap_atac_evidence_from_peaks",
    "chromap_lib_runner_parity",
}


def test_recipe_registry_loads():
    registry = load_recipe_registry()
    assert registry.registry_version == 1
    assert {recipe.id for recipe in registry.recipes} >= PUBLIC_WORKFLOWS
    assert {recipe.id for recipe in registry.recipes} >= PLANNED_RECIPES


def test_enabled_recipes_cover_public_workflows(loaded_default_config):
    workflow_ids = set(get_workflow_schemas())
    enabled = {recipe.workflow_id for recipe in list_recipes(enabled_only=True)}

    assert workflow_ids == PUBLIC_WORKFLOWS
    assert enabled == workflow_ids


def test_enabled_recipe_metadata_is_complete(loaded_default_config):
    for recipe in list_recipes(enabled_only=True):
        assert recipe.workflow_id in get_workflow_schemas()
        assert recipe.command_template
        assert recipe.outputs
        assert recipe.preflight
        assert recipe.smoke_coverage
        assert recipe.docs
        assert recipe.runtime_class in {"smoke", "interactive", "long", "benchmark"}
        assert recipe.benchmark_policy in {
            "serial_required",
            "parallel_safe",
            "not_benchmark",
        }


def test_current_public_recipes_have_handoff_or_output_metadata():
    for recipe_id in PUBLIC_WORKFLOWS:
        recipe = get_recipe(recipe_id)
        assert recipe is not None
        assert recipe.enabled
        assert recipe.outputs
        assert recipe.handoff_artifacts


def test_planned_recipes_are_metadata_only():
    for recipe_id in PLANNED_RECIPES:
        recipe = get_recipe(recipe_id)
        assert recipe is not None
        assert not recipe.enabled
        assert recipe.outputs
        assert recipe.smoke_coverage


def test_recipe_mcp_tools_expose_registry(loaded_default_config):
    listed = list_recipes_tool.fn(enabled_only=True)
    recipe_ids = {recipe["id"] for recipe in listed["recipes"]}
    assert recipe_ids == PUBLIC_WORKFLOWS

    described = describe_recipe.fn("chromap_hic_pairs")
    assert described["id"] == "chromap_hic_pairs"
    assert described["outputs"][0]["artifact_type"] == "pairs"
    assert "hic_stops_at_pairs" in described["preflight"]
