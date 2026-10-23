"""Vectorized Monte Carlo simulation for single-gene autosomal crosses."""

from collections import Counter

import numpy as np

from genetics import InheritanceMode, gamete_probabilities, validate_genotype

ALLOWED_SAMPLE_SIZES = (100, 1_000, 10_000, 100_000)


def _validate_parameters(sample_size: int, seed: int | None) -> None:
    if isinstance(sample_size, bool) or not isinstance(sample_size, (int, np.integer)) or sample_size not in ALLOWED_SAMPLE_SIZES:
        raise ValueError(f"sample_size must be one of {ALLOWED_SAMPLE_SIZES}.")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0 or seed > np.iinfo(np.uint32).max):
        raise ValueError("seed must be None or an integer between 0 and 4,294,967,295.")


def simulate_offspring(parent_one: str, parent_two: str, sample_size: int = 1_000,
                       seed: int | None = None, mode: InheritanceMode | str = InheritanceMode.COMPLETE_DOMINANCE) -> dict:
    """Sample offspring efficiently and aggregate genotype and phenotype frequencies."""
    _validate_parameters(sample_size, seed)
    parent_one, parent_two = validate_genotype(parent_one), validate_genotype(parent_two)
    mode = InheritanceMode(mode)
    rng = np.random.default_rng(seed)

    def sample_alleles(genotype: str) -> np.ndarray:
        alleles = np.array(list(gamete_probabilities(genotype).keys()))
        probabilities = np.array(list(gamete_probabilities(genotype).values()), dtype=float)
        return rng.choice(alleles, size=sample_size, p=probabilities)

    first, second = sample_alleles(parent_one), sample_alleles(parent_two)
    genotypes = np.where(first == second, first + second, "Aa")
    genotype_counts = Counter(genotypes.tolist())
    genotype_frequencies = {key: genotype_counts.get(key, 0) / sample_size for key in ("AA", "Aa", "aa") if genotype_counts.get(key, 0)}
    if mode is InheritanceMode.COMPLETE_DOMINANCE:
        phenotype_values = np.where(genotypes == "aa", "Recessive phenotype", "Dominant phenotype")
    elif mode is InheritanceMode.INCOMPLETE_DOMINANCE:
        phenotype_values = np.select(
            [genotypes == "AA", genotypes == "aa"],
            ["A-allele phenotype", "a-allele phenotype"],
            default="Intermediate phenotype",
        )
    else:
        phenotype_values = np.select(
            [genotypes == "AA", genotypes == "aa"],
            ["A-allele phenotype", "a-allele phenotype"],
            default="Both A- and a-allele phenotypes expressed",
        )
    phenotype_counts = Counter(phenotype_values.tolist())
    phenotype_frequencies = {key: count / sample_size for key, count in phenotype_counts.items()}
    return {
        "sample_size": int(sample_size), "seed": None if seed is None else int(seed),
        "genotypes": genotype_frequencies, "phenotypes": phenotype_frequencies,
        "genotype_counts": dict(genotype_counts), "phenotype_counts": dict(phenotype_counts),
    }
