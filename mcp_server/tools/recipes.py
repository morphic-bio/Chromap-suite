"""Recipe registry loading helpers."""

from pathlib import Path
from typing import Optional

import yaml

from ..schemas.recipe import RecipeEntry, RecipeRegistry


DEFAULT_RECIPE_REGISTRY = Path(__file__).resolve().parents[1] / "recipes" / "registry.yaml"


def load_recipe_registry(path: Optional[Path] = None) -> RecipeRegistry:
    """Load and validate the Chromap recipe registry."""
    registry_path = path or DEFAULT_RECIPE_REGISTRY
    if not registry_path.exists():
        raise FileNotFoundError(f"Recipe registry not found: {registry_path}")
    with open(registry_path) as handle:
        raw = yaml.safe_load(handle)
    return RecipeRegistry(**raw)


def list_recipes(enabled_only: bool = False, path: Optional[Path] = None) -> list[RecipeEntry]:
    """Return recipe registry entries."""
    registry = load_recipe_registry(path)
    if enabled_only:
        return [recipe for recipe in registry.recipes if recipe.enabled]
    return registry.recipes


def get_recipe(recipe_id: str, path: Optional[Path] = None) -> Optional[RecipeEntry]:
    """Return one recipe by id."""
    for recipe in list_recipes(enabled_only=False, path=path):
        if recipe.id == recipe_id:
            return recipe
    return None
