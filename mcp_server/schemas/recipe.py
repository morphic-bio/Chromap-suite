"""Pydantic models for Chromap recipe registry metadata."""

from typing import Any, Optional

from pydantic import BaseModel, Field, model_validator


class RecipeInput(BaseModel):
    """A typed input exposed by a recipe."""

    name: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    type: str = Field(
        pattern=r"^(string|int|float|bool|enum|file|directory|string_list)$"
    )
    required: bool = False
    description: str = ""
    default: Any = None
    choices: Optional[list[str]] = None


class RecipeOutput(BaseModel):
    """An output artifact produced by a recipe."""

    name: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    artifact_type: str = Field(
        pattern=r"^(index|bed|tagalign|bam|cram|fragments|pairs|summary|peaks|names|directory|other)$"
    )
    path_template: str
    required: bool = True
    description: str = ""


class RecipeDocRef(BaseModel):
    """Reference to durable documentation for a recipe."""

    path: str
    description: str = ""


class RecipeEntry(BaseModel):
    """A single recipe registry entry.

    Recipes can be executable (`enabled: true`) or metadata-only placeholders
    for planned work (`enabled: false`).
    """

    id: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    title: str
    purpose: str
    enabled: bool = True
    workflow_id: Optional[str] = Field(
        default=None,
        description="MCP workflow id when this recipe is executable through workflows.",
    )
    runtime_class: str = Field(pattern=r"^(smoke|interactive|long|benchmark)$")
    benchmark_policy: str = Field(
        pattern=r"^(serial_required|parallel_safe|not_benchmark)$"
    )
    command_template: list[str] = Field(
        default_factory=list,
        description="Command argv template; placeholders use {parameter_name}.",
    )
    inputs: list[RecipeInput] = Field(default_factory=list)
    outputs: list[RecipeOutput] = Field(default_factory=list)
    preflight: list[str] = Field(default_factory=list)
    smoke_coverage: dict[str, str] = Field(default_factory=dict)
    docs: list[RecipeDocRef] = Field(default_factory=list)
    handoff_artifacts: list[str] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_enabled_recipe(self) -> "RecipeEntry":
        if self.enabled and not self.workflow_id:
            raise ValueError(f"enabled recipe {self.id} must declare workflow_id")
        if self.enabled and not self.command_template:
            raise ValueError(f"enabled recipe {self.id} must declare command_template")
        if self.enabled and not self.outputs:
            raise ValueError(f"enabled recipe {self.id} must declare outputs")
        if self.enabled and not self.preflight:
            raise ValueError(f"enabled recipe {self.id} must declare preflight rules")
        if self.enabled and not self.smoke_coverage:
            raise ValueError(f"enabled recipe {self.id} must declare smoke coverage")
        if self.enabled and not self.docs:
            raise ValueError(f"enabled recipe {self.id} must declare docs")
        return self


class RecipeRegistry(BaseModel):
    """Top-level recipe registry file."""

    registry_version: int = 1
    recipes: list[RecipeEntry]

    @model_validator(mode="after")
    def _validate_unique_recipe_ids(self) -> "RecipeRegistry":
        ids = [recipe.id for recipe in self.recipes]
        duplicates = sorted({recipe_id for recipe_id in ids if ids.count(recipe_id) > 1})
        if duplicates:
            raise ValueError(f"duplicate recipe ids: {', '.join(duplicates)}")
        return self
