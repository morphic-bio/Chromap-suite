"""Preflight validation for MCP server."""

import re
from pathlib import Path
from typing import Any, Optional

from ..config import get_config
from ..schemas.recipe import RecipeEntry
from ..schemas.config import ScriptConfig
from ..schemas.responses import (
    PreflightCheck,
    PreflightResponse,
    RecipePreflightCheck,
    RecipePreflightResponse,
)
from ..schemas.run_config import RunConfig
from .recipes import get_recipe
from .utils import (
    find_binary,
    get_disk_space_gb,
    is_path_allowed,
    validate_path,
)


def _csv_items(value: Any) -> list[str]:
    """Return comma-separated path-ish values as a list."""
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if str(item)]
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


_SHELL_METACHARS = set(";&|`$<>")


def _check(
    rule_id: str,
    status: str,
    message: str,
    path: Optional[str] = None,
    suggested_fix: Optional[str] = None,
    **details: Any,
) -> RecipePreflightCheck:
    return RecipePreflightCheck(
        rule_id=rule_id,
        status=status,
        message=message,
        path=path,
        suggested_fix=suggested_fix,
        details={k: v for k, v in details.items() if v is not None},
    )


def _validate_existing_file(rule_id: str, params: dict[str, Any], key: str) -> RecipePreflightCheck:
    value = params.get(key)
    if not value:
        return _check(
            rule_id,
            "fail",
            f"Missing required path parameter: {key}",
            suggested_fix=f"Provide {key}.",
        )
    valid, error = validate_path(value, must_exist=True, must_be_file=True, must_be_readable=True)
    if not valid:
        return _check(
            rule_id,
            "fail",
            error or f"Invalid file path for {key}: {value}",
            path=str(value),
            suggested_fix=f"Provide an existing readable file for {key} under a trusted root.",
        )
    return _check(rule_id, "pass", f"{key} exists and is readable", path=str(value))


def _render_template(template: str, params: dict[str, Any]) -> Optional[str]:
    names = set(re.findall(r"{([a-zA-Z_][a-zA-Z0-9_]*)}", template))
    rendered = template
    for name in names:
        value = params.get(name)
        if value in (None, ""):
            return None
        rendered = rendered.replace("{" + name + "}", str(value))
    return rendered


def _output_paths(recipe: RecipeEntry, params: dict[str, Any]) -> list[str]:
    paths: list[str] = []
    for output in recipe.outputs:
        rendered = _render_template(output.path_template, params)
        if rendered:
            paths.append(rendered)
    return paths


def _check_output_parent(recipe: RecipeEntry, params: dict[str, Any]) -> RecipePreflightCheck:
    paths = _output_paths(recipe, params)
    if not paths:
        return _check(
            "output_parent_trusted_or_creatable",
            "fail",
            "No concrete output paths could be derived from recipe parameters",
            suggested_fix="Provide the recipe output path parameters.",
        )

    checked: list[str] = []
    for output_path in paths:
        if any(ch in output_path for ch in _SHELL_METACHARS):
            return _check(
                "output_parent_trusted_or_creatable",
                "fail",
                "Output path contains shell metacharacters",
                path=str(output_path),
                suggested_fix="Use literal output paths only; do not pass shell syntax through recipe parameters.",
                outputs=paths,
            )
        parent = Path(output_path).parent
        valid, error = validate_path(parent, must_be_writable=True)
        if not valid:
            return _check(
                "output_parent_trusted_or_creatable",
                "fail",
                error or f"Output parent is not writable: {parent}",
                path=str(parent),
                suggested_fix="Use an output path under a trusted writable artifact directory.",
                outputs=paths,
            )
        checked.append(str(parent))

    return _check(
        "output_parent_trusted_or_creatable",
        "pass",
        "Output parent directories are trusted and writable or creatable",
        details_checked=checked,
        outputs=paths,
    )


def _check_read_pair_lists_match(params: dict[str, Any]) -> RecipePreflightCheck:
    read1 = _csv_items(params.get("read1"))
    read2 = _csv_items(params.get("read2"))
    if not read1 or not read2:
        return _check(
            "read_pair_lists_match",
            "fail",
            "read1 and read2 are both required for paired-end recipes",
            suggested_fix="Provide read1 and read2 FASTQ path lists.",
        )
    if len(read1) != len(read2):
        return _check(
            "read_pair_lists_match",
            "fail",
            f"read1/read2 lane counts differ: {len(read1)} vs {len(read2)}",
            suggested_fix="Provide the same number of comma-separated read1 and read2 FASTQs.",
            read1=read1,
            read2=read2,
        )
    return _check(
        "read_pair_lists_match",
        "pass",
        f"read1/read2 lane counts match ({len(read1)})",
        read1_count=len(read1),
        read2_count=len(read2),
    )


def _check_barcode_list_matches_read1(params: dict[str, Any]) -> RecipePreflightCheck:
    read1 = _csv_items(params.get("read1"))
    barcode = _csv_items(params.get("barcode"))
    if not barcode:
        return _check(
            "barcode_list_matches_read1",
            "fail",
            "barcode FASTQ list is required for this recipe",
            suggested_fix="Provide barcode FASTQ path(s).",
        )
    if len(read1) != len(barcode):
        return _check(
            "barcode_list_matches_read1",
            "fail",
            f"read1/barcode lane counts differ: {len(read1)} vs {len(barcode)}",
            suggested_fix="Provide one barcode FASTQ for each read1 FASTQ.",
            read1_count=len(read1),
            barcode_count=len(barcode),
        )
    return _check(
        "barcode_list_matches_read1",
        "pass",
        f"read1/barcode lane counts match ({len(read1)})",
        read1_count=len(read1),
        barcode_count=len(barcode),
    )


def _check_secondary_output_not_primary(params: dict[str, Any]) -> RecipePreflightCheck:
    primary = params.get("output")
    secondary_candidates = [
        params.get("atac_fragments"),
        params.get("noY_output"),
        params.get("Y_output"),
        params.get("y_read_names_output"),
    ]
    for secondary in secondary_candidates:
        if primary and secondary and Path(str(primary)) == Path(str(secondary)):
            return _check(
                "secondary_output_not_primary",
                "fail",
                "Primary and secondary output paths must differ",
                path=str(primary),
                suggested_fix="Choose distinct output paths.",
            )
    return _check("secondary_output_not_primary", "pass", "Primary and secondary outputs differ")


def _check_atac_fragments_requires_barcoded_paired_bam(params: dict[str, Any]) -> RecipePreflightCheck:
    missing = [
        key
        for key in ("read1", "read2", "barcode", "output", "atac_fragments")
        if not params.get(key)
    ]
    if missing:
        return _check(
            "atac_fragments_requires_barcoded_paired_bam",
            "fail",
            f"Missing ATAC fragments requirements: {', '.join(missing)}",
            suggested_fix="Provide paired reads, barcode reads, BAM/CRAM output, and fragments path.",
            missing=missing,
        )
    output = str(params["output"]).lower()
    if not (output.endswith(".bam") or output.endswith(".cram")):
        return _check(
            "atac_fragments_requires_barcoded_paired_bam",
            "fail",
            "--atac-fragments requires BAM or CRAM primary output",
            path=str(params["output"]),
            suggested_fix="Use an output path ending in .bam or .cram.",
        )
    return _check(
        "atac_fragments_requires_barcoded_paired_bam",
        "pass",
        "ATAC fragments recipe has paired reads, barcode reads, and BAM/CRAM output",
    )


def _check_hic_stops_at_pairs(params: dict[str, Any]) -> RecipePreflightCheck:
    output = str(params.get("output", ""))
    if output and not (output.endswith(".pairs") or output.endswith(".pairs.gz")):
        return _check(
            "hic_stops_at_pairs",
            "fail",
            "Hi-C Chromap recipe should emit .pairs or .pairs.gz only",
            path=output,
            suggested_fix="Use an output path ending in .pairs or .pairs.gz.",
        )
    return _check("hic_stops_at_pairs", "pass", "Hi-C recipe stops at pairs output")


def _check_write_index_requires_sort_bam(params: dict[str, Any]) -> RecipePreflightCheck:
    if _truthy(params.get("write_index")) and not _truthy(params.get("sort_bam")):
        return _check(
            "write_index_requires_sort_bam",
            "fail",
            "--write-index requires --sort-bam",
            suggested_fix="Enable sort_bam or disable write_index.",
        )
    return _check("write_index_requires_sort_bam", "pass", "write_index/sort_bam settings are compatible")


def _check_macs3_required_path(rule_id: str, params: dict[str, Any], key: str) -> RecipePreflightCheck:
    if not params.get(key):
        return _check(
            rule_id,
            "fail",
            f"Missing required MACS3 output parameter: {key}",
            suggested_fix=f"Provide {key}.",
        )
    return _check(rule_id, "pass", f"{key} is set", path=str(params[key]))


def _check_macs3_pvalue_positive(params: dict[str, Any]) -> RecipePreflightCheck:
    raw = params.get("macs3_frag_pvalue", params.get("pvalue", 0.01))
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return _check(
            "macs3_pvalue_positive",
            "fail",
            f"MACS3 p-value is not numeric: {raw}",
            suggested_fix="Use a positive numeric p-value.",
        )
    if value <= 0:
        return _check(
            "macs3_pvalue_positive",
            "fail",
            "MACS3 p-value must be positive",
            suggested_fix="Use a p-value greater than 0.",
            value=value,
        )
    return _check("macs3_pvalue_positive", "pass", "MACS3 p-value is positive", value=value)


def _check_y_noy_requires_sam_bam_cram(params: dict[str, Any]) -> RecipePreflightCheck:
    output = str(params.get("output", "")).lower()
    if output and not (
        output.endswith(".sam")
        or output.endswith(".sam.gz")
        or output.endswith(".bam")
        or output.endswith(".cram")
    ):
        return _check(
            "y_noy_requires_sam_bam_cram",
            "fail",
            "Y/noY split requires SAM, BAM, or CRAM output",
            path=output,
            suggested_fix="Use .sam, .sam.gz, .bam, or .cram output.",
        )
    return _check("y_noy_requires_sam_bam_cram", "pass", "Y/noY output format is compatible")


def _check_binary_exists(rule_id: str, binary_name: str) -> RecipePreflightCheck:
    path = find_binary(binary_name)
    if path is None:
        return _check(
            rule_id,
            "fail",
            f"Required binary not found under trusted roots: {binary_name}",
            suggested_fix=f"Build or install {binary_name} under a trusted root.",
        )
    return _check(rule_id, "pass", f"{binary_name} found", path=str(path))


def _check_input_paths_trusted(recipe: RecipeEntry, params: dict[str, Any]) -> RecipePreflightCheck:
    path_values: list[str] = []
    for input_def in recipe.inputs:
        if input_def.name not in params:
            continue
        value = params.get(input_def.name)
        if value in (None, ""):
            continue
        if input_def.type in ("file", "directory", "string_list"):
            path_values.extend(_csv_items(value))
        elif input_def.type == "string" and input_def.name in {"read1", "read2", "barcode"}:
            path_values.extend(_csv_items(value))

    for raw_path in path_values:
        if any(ch in raw_path for ch in _SHELL_METACHARS):
            return _check(
                "input_paths_trusted",
                "fail",
                "Input path contains shell metacharacters",
                path=raw_path,
                suggested_fix="Use literal file paths only; do not pass shell syntax through recipe parameters.",
            )
        if not is_path_allowed(raw_path):
            return _check(
                "input_paths_trusted",
                "fail",
                "Input path is outside trusted roots",
                path=raw_path,
                suggested_fix="Place inputs under a configured trusted root.",
            )
    return _check(
        "input_paths_trusted",
        "pass",
        "Recipe input paths are under trusted roots and contain no shell metacharacters",
        path_count=len(path_values),
    )


def _run_recipe_rule(recipe: RecipeEntry, params: dict[str, Any], rule_id: str) -> RecipePreflightCheck:
    if rule_id == "reference_fasta_exists":
        return _validate_existing_file(rule_id, params, "reference")
    if rule_id == "chromap_index_exists":
        return _validate_existing_file(rule_id, params, "index")
    if rule_id == "output_parent_trusted_or_creatable":
        return _check_output_parent(recipe, params)
    if rule_id == "read_pair_lists_match":
        return _check_read_pair_lists_match(params)
    if rule_id == "barcode_list_matches_read1":
        return _check_barcode_list_matches_read1(params)
    if rule_id == "secondary_output_not_primary":
        return _check_secondary_output_not_primary(params)
    if rule_id == "atac_fragments_requires_barcoded_paired_bam":
        return _check_atac_fragments_requires_barcoded_paired_bam(params)
    if rule_id == "hic_stops_at_pairs":
        return _check_hic_stops_at_pairs(params)
    if rule_id == "write_index_requires_sort_bam":
        return _check_write_index_requires_sort_bam(params)
    if rule_id == "macs3_peaks_output_required":
        return _check_macs3_required_path(rule_id, params, "macs3_frag_peaks_output")
    if rule_id == "macs3_summits_output_required":
        return _check_macs3_required_path(rule_id, params, "macs3_frag_summits_output")
    if rule_id == "macs3_pvalue_positive":
        return _check_macs3_pvalue_positive(params)
    if rule_id == "y_noy_requires_sam_bam_cram":
        return _check_y_noy_requires_sam_bam_cram(params)
    if rule_id == "chromap_binary_exists":
        return _check_binary_exists(rule_id, "chromap")
    if rule_id == "chromap_lib_runner_exists":
        return _check_binary_exists(rule_id, "chromap_lib_runner")

    return _check(
        rule_id,
        "fail",
        f"Unknown recipe preflight rule: {rule_id}",
        suggested_fix="Implement this rule before enabling recipe execution.",
    )


def preflight_recipe(recipe_id: str, params: dict[str, Any]) -> RecipePreflightResponse:
    """Run recipe-driven preflight checks without executing Chromap."""
    recipe = get_recipe(recipe_id)
    if recipe is None:
        raise ValueError(f"Unknown recipe: {recipe_id}")

    checks: list[RecipePreflightCheck] = []
    supplied = dict(params or {})

    for input_def in recipe.inputs:
        if input_def.required and supplied.get(input_def.name) in (None, ""):
            checks.append(
                _check(
                    f"required_input:{input_def.name}",
                    "fail",
                    f"Missing required input: {input_def.name}",
                    suggested_fix=f"Provide {input_def.name}.",
                )
            )

    checks.append(_check_input_paths_trusted(recipe, supplied))

    for rule_id in recipe.preflight:
        checks.append(_run_recipe_rule(recipe, supplied, rule_id))

    warnings = [check.message for check in checks if check.status == "warn"]
    errors = [check.message for check in checks if check.status == "fail"]
    return RecipePreflightResponse(
        recipe_id=recipe.id,
        valid=not errors,
        checks=checks,
        warnings=warnings,
        errors=errors,
    )


def check_script_allowed(run_config: RunConfig) -> PreflightCheck:
    """Check if the requested script is in the allowlist.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None:
        allowed_scripts = [s.name for s in config.scripts]
        return PreflightCheck(
            name="script_allowed",
            passed=False,
            message=f"Script '{run_config.script}' is not in the allowlist",
            details={
                "requested": run_config.script,
                "allowed": allowed_scripts,
            },
        )

    # Check if script is marked as runnable
    if not script.runnable:
        return PreflightCheck(
            name="script_allowed",
            passed=False,
            message=f"Script '{run_config.script}' is not currently runnable",
            details={"reason": "Script marked as not runnable in config"},
        )

    return PreflightCheck(
        name="script_allowed",
        passed=True,
        details={"script": script.name, "module": script.module},
    )


def check_script_path(run_config: RunConfig) -> PreflightCheck:
    """Check if the script file exists and is within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None:
        return PreflightCheck(
            name="script_path",
            passed=False,
            message="Script not found in config",
        )

    # Resolve script path
    script_path = Path(script.path)
    if not script_path.is_absolute():
        script_path = config.paths.repo_root / script_path

    # For build commands (like 'make'), check if the command exists
    if script.path in ("make", "cmake", "ninja"):
        # find_binary enforces trusted roots by default
        binary = find_binary(script.path, enforce_trusted_roots=True)
        if binary:
            return PreflightCheck(
                name="script_path",
                passed=True,
                details={"path": str(binary), "type": "build_command"},
            )
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Build command '{script.path}' not found in trusted roots",
        )

    # Check script path is within trusted roots
    if not is_path_allowed(script_path):
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Script path outside trusted roots: {script_path}",
            details={"path": str(script_path)},
        )

    # Check script file exists
    if not script_path.exists():
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Script file not found: {script_path}",
            details={"path": str(script_path)},
        )

    # Check script is a file, not a directory
    if not script_path.is_file():
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Script path is not a file: {script_path}",
            details={"path": str(script_path)},
        )

    return PreflightCheck(
        name="script_path",
        passed=True,
        details={"path": str(script_path)},
    )


def check_working_dir(run_config: RunConfig) -> PreflightCheck:
    """Check if the script's working_dir is within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None:
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message="Script not found in config",
        )

    # If no working_dir specified, defaults to repo_root which is trusted
    if not script.working_dir:
        return PreflightCheck(
            name="working_dir",
            passed=True,
            message="No working_dir specified (defaults to repo_root)",
            details={"path": str(config.paths.repo_root)},
        )

    # Resolve working_dir path
    working_dir = Path(script.working_dir)
    if not working_dir.is_absolute():
        working_dir = config.paths.repo_root / working_dir

    # Check working_dir is within trusted roots
    if not is_path_allowed(working_dir):
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message=f"Script working_dir outside trusted roots: {working_dir}",
            details={"path": str(working_dir), "script": script.name},
        )

    # Check working_dir exists
    if not working_dir.exists():
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message=f"Script working_dir does not exist: {working_dir}",
            details={"path": str(working_dir), "script": script.name},
        )

    # Check working_dir is a directory
    if not working_dir.is_dir():
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message=f"Script working_dir is not a directory: {working_dir}",
            details={"path": str(working_dir), "script": script.name},
        )

    return PreflightCheck(
        name="working_dir",
        passed=True,
        details={"path": str(working_dir)},
    )


def check_dataset(run_config: RunConfig) -> PreflightCheck:
    """Check if the requested dataset exists and is within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    if not run_config.dataset_id:
        return PreflightCheck(
            name="dataset",
            passed=True,
            message="No dataset specified",
        )

    config = get_config()
    dataset = config.get_dataset(run_config.dataset_id)

    if dataset is None:
        available = [d.id for d in config.datasets]
        return PreflightCheck(
            name="dataset",
            passed=False,
            message=f"Dataset '{run_config.dataset_id}' not found",
            details={"requested": run_config.dataset_id, "available": available},
        )

    dataset_path = Path(dataset.path)

    # Check dataset path is within trusted roots
    if not is_path_allowed(dataset_path):
        return PreflightCheck(
            name="dataset",
            passed=False,
            message=f"Dataset path outside trusted roots: {dataset.path}",
            details={"dataset_id": dataset.id, "path": str(dataset.path)},
        )

    # Check dataset path exists
    if not dataset_path.exists():
        return PreflightCheck(
            name="dataset",
            passed=False,
            message=f"Dataset path does not exist: {dataset.path}",
            details={"dataset_id": dataset.id, "path": str(dataset.path)},
        )

    return PreflightCheck(
        name="dataset",
        passed=True,
        details={"dataset_id": dataset.id, "path": str(dataset.path)},
    )


def check_output_directory(run_config: RunConfig) -> PreflightCheck:
    """Check if the output directory is writable.

    This check validates that:
    1. The path is within trusted roots
    2. If the directory exists, it's writable
    3. If it doesn't exist, a writable parent exists in the chain

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()

    # Determine output directory
    if run_config.out_dir:
        out_dir = Path(run_config.out_dir)
    else:
        out_dir = config.paths.temp_root / "mcp_runs"

    # Check path is within trusted roots
    if not is_path_allowed(out_dir):
        return PreflightCheck(
            name="output_directory",
            passed=False,
            message=f"Output directory outside trusted roots: {out_dir}",
            details={"path": str(out_dir)},
        )

    # Check if directory exists or can be created
    valid, error = validate_path(out_dir, must_be_writable=True)
    if not valid:
        return PreflightCheck(
            name="output_directory",
            passed=False,
            message=error or f"Output directory not writable: {out_dir}",
            details={"path": str(out_dir)},
        )

    # Additional check: verify we can actually create missing parent chain
    # validate_path finds an existing ancestor, but we should report how deep that is
    if not out_dir.exists():
        ancestor = out_dir.parent
        depth = 1
        while not ancestor.exists() and ancestor.parent != ancestor:
            ancestor = ancestor.parent
            depth += 1

        if depth > 3:
            # Warn if the directory chain is deep - may indicate a typo
            return PreflightCheck(
                name="output_directory",
                passed=True,
                message=f"Output directory will be created ({depth} levels deep)",
                details={
                    "path": str(out_dir),
                    "existing_ancestor": str(ancestor),
                    "depth": depth,
                    "warning": "Deep directory chain - verify path is correct",
                },
            )

    return PreflightCheck(
        name="output_directory",
        passed=True,
        details={"path": str(out_dir), "exists": out_dir.exists()},
    )


def check_binaries(run_config: RunConfig) -> PreflightCheck:
    """Check if required binaries are present.

    Skips this check for scripts in the "build" module since those scripts
    are meant to build the required binaries.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()

    # Get the script to check its module
    script = config.get_script(run_config.script)

    # Skip binary checks for build scripts - they create the binaries
    if script and script.module == "build":
        return PreflightCheck(
            name="binaries_present",
            passed=True,
            message="Skipped for build scripts",
            details={"reason": "Build scripts create binaries, not require them"},
        )

    missing: list[str] = []
    found: dict[str, str] = {}

    for binary_config in config.required_binaries:
        binary_path = find_binary(binary_config.name)
        if binary_path:
            found[binary_config.name] = str(binary_path)
        else:
            missing.append(binary_config.name)

    if missing:
        return PreflightCheck(
            name="binaries_present",
            passed=False,
            message=f"Required binaries not found: {', '.join(missing)}",
            details={"missing": missing, "found": found},
        )

    return PreflightCheck(
        name="binaries_present",
        passed=True,
        details=found,
    )


def check_disk_space(run_config: RunConfig) -> PreflightCheck:
    """Check if there is sufficient disk space.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()

    # Determine output directory for space check
    if run_config.out_dir:
        check_path = Path(run_config.out_dir)
    else:
        check_path = config.paths.temp_root

    # Get available space
    available_gb = get_disk_space_gb(check_path)

    # Determine required space based on script's test suite
    script = config.get_script(run_config.script)
    required_gb = config.disk_space.default_min_gb

    if script:
        # Find the test suite this script belongs to
        for suite in config.test_suites:
            if script.name in suite.scripts:
                if suite.disk_space_min_gb is not None:
                    required_gb = suite.disk_space_min_gb
                break

    if available_gb < required_gb:
        return PreflightCheck(
            name="disk_space",
            passed=False,
            message=f"Insufficient disk space: {available_gb:.1f}GB available, {required_gb}GB required",
            details={
                "available_gb": round(available_gb, 1),
                "required_gb": required_gb,
                "check_path": str(check_path),
            },
        )

    return PreflightCheck(
        name="disk_space",
        passed=True,
        details={
            "available_gb": round(available_gb, 1),
            "required_gb": required_gb,
        },
    )


def check_fixtures(run_config: RunConfig) -> PreflightCheck:
    """Check if required fixtures are present and within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None or not script.fixtures:
        return PreflightCheck(
            name="fixtures_present",
            passed=True,
            message="No fixtures required",
        )

    missing: list[str] = []
    found: list[str] = []
    outside_roots: list[str] = []

    for fixture_path in script.fixtures:
        path = Path(fixture_path)

        # Check fixture is within trusted roots first
        if not is_path_allowed(path):
            outside_roots.append(fixture_path)
        elif path.exists():
            found.append(fixture_path)
        else:
            missing.append(fixture_path)

    # Report paths outside trusted roots as a security error
    if outside_roots:
        return PreflightCheck(
            name="fixtures_present",
            passed=False,
            message=f"Fixture paths outside trusted roots: {len(outside_roots)}",
            details={"outside_roots": outside_roots, "found": found, "missing": missing},
        )

    if missing:
        return PreflightCheck(
            name="fixtures_present",
            passed=False,
            message=f"Required fixtures not found: {len(missing)} missing",
            details={"missing": missing, "found": found},
        )

    return PreflightCheck(
        name="fixtures_present",
        passed=True,
        details={"fixtures": found},
    )


def check_paths_valid(run_config: RunConfig) -> PreflightCheck:
    """Check that all specified paths are within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    invalid_paths: list[str] = []

    # Check out_dir if specified
    if run_config.out_dir:
        if not is_path_allowed(run_config.out_dir):
            invalid_paths.append(run_config.out_dir)

    if invalid_paths:
        return PreflightCheck(
            name="paths_valid",
            passed=False,
            message="Paths outside trusted roots",
            details={"invalid_paths": invalid_paths},
        )

    return PreflightCheck(
        name="paths_valid",
        passed=True,
    )


def run_preflight(run_config: RunConfig) -> PreflightResponse:
    """Run all preflight checks.

    Args:
        run_config: Run configuration to validate.

    Returns:
        PreflightResponse with all check results.
    """
    checks: list[PreflightCheck] = []
    warnings: list[str] = []
    errors: list[str] = []

    # Run all checks
    check_functions = [
        check_script_allowed,
        check_script_path,
        check_working_dir,  # Validate script's working_dir is within trusted roots
        check_dataset,
        check_output_directory,
        check_binaries,
        check_disk_space,
        check_fixtures,
        check_paths_valid,
    ]

    for check_fn in check_functions:
        try:
            result = check_fn(run_config)
            checks.append(result)
            if not result.passed:
                errors.append(result.message or f"{result.name} check failed")
        except Exception as e:
            checks.append(
                PreflightCheck(
                    name=check_fn.__name__.replace("check_", ""),
                    passed=False,
                    message=f"Check error: {str(e)}",
                )
            )
            errors.append(f"{check_fn.__name__}: {str(e)}")

    # Determine overall validity
    valid = all(check.passed for check in checks)

    return PreflightResponse(
        valid=valid,
        checks=checks,
        warnings=warnings,
        errors=errors,
    )
