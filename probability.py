"""Exact theoretical distributions and comparison utilities."""

from collections.abc import Mapping

import pandas as pd

from genetics import InheritanceMode, genotype_probabilities, phenotype_probabilities


def theoretical_distributions(parent_one: str, parent_two: str, mode: InheritanceMode | str) -> tuple[dict[str, float], dict[str, float]]:
    """Return exact genotype and phenotype probability distributions."""
    genotypes = genotype_probabilities(parent_one, parent_two)
    return genotypes, phenotype_probabilities(genotypes, mode)


def comparison_frame(theoretical: Mapping[str, float], observed: Mapping[str, float], category: str) -> pd.DataFrame:
    """Build a tidy comparison table, including categories absent from observations."""
    categories = list(dict.fromkeys([*theoretical.keys(), *observed.keys()]))
    rows = [{
        "category": category_name,
        "theoretical_probability": float(theoretical.get(category_name, 0.0)),
        "observed_frequency": float(observed.get(category_name, 0.0)),
        "absolute_error": abs(float(observed.get(category_name, 0.0)) - float(theoretical.get(category_name, 0.0))),
    } for category_name in categories]
    return pd.DataFrame(rows)
