# Workflow Schema Authoring Guide

How to add a new workflow to the Chromap-suite MCP server.

## Overview

A workflow is a structured parameter contract for a shell script. Adding one
requires three tracked metadata updates:

1. a workflow YAML schema,
2. a `mcp_server/config.yaml` workflow entry,
3. a `mcp_server/recipes/registry.yaml` recipe entry.

No code changes are needed for a basic command-rendering recipe.

## Step 1: Create the schema YAML

Create a file in `mcp_server/workflows/<workflow_id>.yaml`.

### Minimal example

```yaml
id: "my_workflow"
title: "My Workflow"
summary: "One-line description of what this workflow does."
entry_script: "scripts/my_workflow.sh"

parameters:
  - name: "input_dir"
    cli_flag: "--input-dir"
    type: "directory"
    required: true
    description: "Input data directory."
    path_must_exist: true
    category: "inputs"

  - name: "threads"
    cli_flag: "--threads"
    type: "int"
    default: 8
    description: "Thread count."
    category: "performance"

  - name: "dry_run"
    cli_flag: "--dry-run"
    type: "bool"
    default: false
    description: "Dry run mode."
    category: "mode"

constraints:
  - kind: "positive"
    params: ["threads"]
    message: "Must be a positive integer."

rendering:
  flag_order: ["input_dir", "threads", "dry_run"]
```

### Full field reference

#### Top-level fields

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `id` | yes | | Unique workflow identifier (snake_case) |
| `title` | yes | | Human-readable display name |
| `summary` | yes | | One-line description |
| `kind` | no | `"shell_workflow"` | Workflow kind |
| `entry_script` | yes | | Path to entry script (relative to repo root) |
| `supported_modes` | no | `["local", "dry-run-capable"]` | Execution modes |
| `caveats` | no | `[]` | Known caveats (shown to agents) |
| `default_output_layout` | no | `""` | Description of output directory structure |
| `stages` | no | `[]` | Semantic stages within the workflow |
| `parameter_groups` | no | `[]` | Ordered groups for display |
| `parameters` | no | `[]` | Parameter definitions |
| `constraints` | no | `[]` | Cross-parameter constraints |
| `rendering` | no | defaults | Rendering configuration |

#### Parameter fields

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `name` | yes | | Identifier (snake_case, must be unique within workflow) |
| `cli_flag` | yes | | CLI flag string (e.g. `--my-flag`) |
| `type` | yes | | One of: `string`, `int`, `float`, `bool`, `enum`, `file`, `directory`, `string_list` |
| `required` | no | `false` | Whether the parameter must be supplied |
| `default` | no | `null` | Default value; `null` means no default |
| `description` | no | `""` | Human-readable description |
| `choices` | no | `null` | Allowed values (only for `enum` type) |
| `repeatable` | no | `false` | Whether the flag can appear multiple times |
| `path_must_exist` | no | `false` | For `file`/`directory` types: validate existence at validation time |
| `category` | no | `"general"` | Logical group (matches `parameter_groups.name`) |
| `stage` | no | `"top_level"` | Which workflow stage owns this parameter |
| `source` | no | `"workflow_wrapper"` | Which script layer defines this parameter |
| `env_var` | no | `null` | If set, export the value as this env var in the rendered command |
| `skip_when_default` | no | `false` | If true, omit the CLI flag when the value equals the schema default |
| `is_output_root` | no | `false` | If true, report this parameter's value as the workflow `output_root` |

#### Parameter types

| Type | CLI rendering | Validation |
|------|--------------|------------|
| `string` | `--flag value` | Non-empty string |
| `int` | `--flag 42` | Must parse as integer |
| `float` | `--flag 0.5` | Must parse as float |
| `bool` | `--flag` (when true, omitted when false) | Boolean-like value |
| `enum` | `--flag value` | Must be in `choices` list |
| `file` | `--flag /path/to/file` | Optionally checked for existence |
| `directory` | `--flag /path/to/dir` | Optionally checked for existence |
| `string_list` | `--flag "a,b,c"` | CSV or list |

#### Rendering hints

These three fields on parameters control how the renderer handles special cases
without hardcoding workflow-specific logic:

**`env_var`**: When set, the parameter's value is added to `env_overrides` in the
render response under the given name. Use this for parameters that the entry
script exports as environment variables rather than passing as CLI flags.

```yaml
  - name: "my_seed"
    cli_flag: "--seed"
    type: "int"
    default: 1
    env_var: "MY_SEED"  # Will appear in env_overrides as MY_SEED=<value>
```

**`skip_when_default`**: When true, the CLI flag is omitted from the rendered
argv when the value equals the schema default. Use this for parameters where
the default value means "disabled" and the script already handles the default
internally.

```yaml
  - name: "downsample_reads"
    cli_flag: "--downsample-reads"
    type: "int"
    default: 0              # 0 means no downsampling
    skip_when_default: true  # Don't emit --downsample-reads 0
```

**`is_output_root`**: Marks which parameter holds the workflow's output
directory. Its value is returned as `output_root` in the render response so
agents can find where results will land.

```yaml
  - name: "out_dir"
    cli_flag: "--out-dir"
    type: "directory"
    is_output_root: true
```

#### Constraint kinds

| Kind | Behavior |
|------|----------|
| `mutual_exclusion` | Error if more than one of `params` is set |
| `group_required` | Warn if some but not all of `params` are set |
| `dependency` | Warn if `params[0]` is set but `params[1]` makes it irrelevant |
| `positive` | Error if any of `params` is <= 0 |
| `non_negative` | Error if any of `params` is < 0 |

#### Rendering configuration

| Field | Default | Description |
|-------|---------|-------------|
| `bool_style` | `"flag_only"` | How booleans render (currently only `flag_only` is supported) |
| `flag_order` | `[]` | Deterministic emission order (by param name); if empty, uses declaration order |
| `omit_absent_optionals` | `true` | Skip flags for optional params with no value |
| `csv_style` | `"quoted"` | How `string_list` values are joined |

#### Stages

Stages are informational metadata describing the workflow's semantic phases.
They don't affect validation or rendering but help agents understand the
workflow structure.

```yaml
stages:
  - name: "alignment"
    title: "Chromap mapping"
    script: "scripts/align.sh"
    description: "Run Chromap mapping."
  - name: "qc"
    title: "Quality Control"
    script: "scripts/qc.sh"
    gated_by: "skip_qc"  # This stage is skipped when skip_qc=true
```

#### Parameter groups

Groups control display ordering. Each group has a name, title, and a list of
parameter names. Parameters can only belong to one group.

```yaml
parameter_groups:
  - name: "inputs"
    title: "Input Paths"
    parameters: ["input_dir", "reference"]
  - name: "performance"
    title: "Performance"
    parameters: ["threads", "memory"]
```

## Step 2: Register in config.yaml

Add an entry under the `workflows:` section:

```yaml
workflows:
  - id: "my_workflow"
    title: "My Workflow"
    summary: "One-line description."
    entry_script: "scripts/my_workflow.sh"
    kind: "shell_workflow"
    schema_file: "mcp_server/workflows/my_workflow.yaml"
```

## Step 3: Register in the recipe registry

Add a matching entry to `mcp_server/recipes/registry.yaml`.

The registry is the shared metadata layer intended for MCP tools, Launchpad,
docs, tests, and future preflight/run-manifest work. A workflow YAML says how
to render a command; a recipe entry says why the command exists, what artifacts
it produces, how it is tested, and what safety checks apply.

Minimal executable recipe example:

```yaml
  - id: my_workflow
    title: My Workflow
    purpose: One-line description of the operator intent.
    enabled: true
    workflow_id: my_workflow
    runtime_class: interactive
    benchmark_policy: not_benchmark
    command_template:
      - scripts/my_workflow.sh
      - --input-dir
      - "{input_dir}"
    inputs:
      - name: input_dir
        type: directory
        required: true
        description: Input data directory.
    outputs:
      - name: result_dir
        artifact_type: directory
        path_template: "{out_dir}"
        description: Workflow output directory.
    preflight:
      - input_dir_exists
      - output_parent_trusted_or_creatable
    smoke_coverage:
      S0: tests/run_my_workflow_smoke.sh
    docs:
      - path: mcp_server/workflows/my_workflow.yaml
        description: Workflow schema.
```

Registry field policy:

| Field | Required for enabled recipes | Description |
|-------|------------------------------|-------------|
| `id` | yes | Stable recipe id; usually matches workflow id. |
| `title` | yes | Human-readable label. |
| `purpose` | yes | Operator intent, not just command syntax. |
| `enabled` | yes | `true` exposes the recipe as executable metadata. |
| `workflow_id` | yes | Workflow schema id for executable recipes. |
| `runtime_class` | yes | `smoke`, `interactive`, `long`, or `benchmark`. |
| `benchmark_policy` | yes | `not_benchmark`, `serial_required`, or `parallel_safe`. |
| `command_template` | yes | Argv template, not a shell string. |
| `inputs` | yes | Typed input metadata. |
| `outputs` | yes | Expected artifacts and path templates. |
| `preflight` | yes | Preflight rule ids. Rules may be implemented later. |
| `smoke_coverage` | yes | S0/S1/S2 test coverage reference. |
| `docs` | yes | Durable documentation references. |
| `handoff_artifacts` | recommended | STAR Suite or downstream handoff artifacts. |

Metadata-only planned recipes may set `enabled: false`. They still need output
and smoke-coverage metadata so the roadmap stays testable.

The MCP server will automatically pick up workflow YAML changes on the next
config load (or call `reload_config`). The recipe registry is loaded separately
by `mcp_server.tools.recipes`.

Execution safety rules for recipes:

- `command_template` is an argv template, not a shell command line.
- Enabled Launchpad execution requires a matching `workflow_id`.
- Recipe preflight must cover required inputs, output parents, and any
  mode-specific invariants before a workflow is marked executable.
- Path-like inputs and outputs must remain under trusted roots.
- Long or benchmark recipes must use `runtime_class: long` or
  `runtime_class: benchmark`; Launchpad requires explicit opt-in before running
  them.
- Handoff points to STAR Suite or downstream tools belong in
  `handoff_artifacts`, not in hidden UI text.

## Step 4: Verify

```python
# List workflows -- your new one should appear
client.call_tool("list_workflows", {})

# Get the parameter schema
client.call_tool("get_workflow_parameter_schema", {
    "workflow_id": "my_workflow"
})

# Validate params
client.call_tool("validate_workflow_parameters", {
    "workflow_id": "my_workflow",
    "params": {"input_dir": "/path/to/data", "threads": 4},
    "check_paths": True,
})

# Render command
client.call_tool("render_workflow_command", {
    "workflow_id": "my_workflow",
    "params": {"input_dir": "/path/to/data", "threads": 4, "dry_run": True},
})
```

Run the registry tests as well:

```bash
python3 -m pytest mcp_server/tests/test_recipe_registry.py -q
python3 -m pytest mcp_server/tests -q
```

## Step 4: Add tests (recommended)

Create `mcp_server/tests/test_<workflow_id>_e2e.py` following the pattern in
`test_ucsf_workflow_e2e.py`:

1. Use the real schema YAML with temporary directories and stub executables
2. Test the full contract: schema load -> validate -> render -> assert argv
3. Use `dry_run=true` -- never execute the real workflow in tests

## Checklist

- [ ] Schema YAML created in `mcp_server/workflows/`
- [ ] Config entry added to `config.yaml` under `workflows:`
- [ ] `list_workflows` returns the new workflow
- [ ] `get_workflow_parameter_schema` returns correct params, groups, constraints
- [ ] `validate_workflow_parameters` catches missing required params, bad types, constraint violations
- [ ] `render_workflow_command` produces correct argv with deterministic flag order
- [ ] E2E test passes with temp fixtures and dry-run mode
