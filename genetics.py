"""Core genetics types and allele-transmission rules for a single autosomal locus."""

from collections import Counter
from enum import Enum
from itertools import product
from typing import Mapping


class InheritanceMode(str, Enum):
    """Supported educational phenotype mappings."""

    COMPLETE_DOMINANCE = "Complete dominance"
    INCOMPLETE_DOMINANCE = "Incomplete dominance"
    CODOMINANCE = "Codominance"


GENOTYPES = ("AA", "Aa", "aa")


def validate_genotype(genotype: str) -> str:
    """Validate and canonicalize a diploid genotype at the modeled locus."""
    if not isinstance(genotype, str) or len(genotype) != 2:
        raise ValueError("Genotype must be a two-character string: AA, Aa, or aa.")
    if any(allele not in "Aa" for allele in genotype):
        raise ValueError("Alleles must be the hypothetical dominant A or recessive a.")
    canonical = "".join(sorted(genotype, key=lambda allele: allele.islower()))
    if canonical not in GENOTYPES:
        raise ValueError("Genotype must be AA, Aa, or aa.")
    return canonical


def gamete_probabilities(genotype: str) -> dict[str, float]:
    """Return the probability of transmitting each allele from a parent."""
    genotype = validate_genotype(genotype)
    counts = Counter(genotype)
    return {allele: count / 2 for allele, count in sorted(counts.items(), key=lambda pair: pair[0].islower())}


def combine_gametes(allele_one: str, allele_two: str) -> str:
    """Combine one allele from each parent into canonical genotype notation."""
    if allele_one not in "Aa" or len(allele_one) != 1 or allele_two not in "Aa" or len(allele_two) != 1:
        raise ValueError("Each gamete must contain exactly one A or a allele.")
    return validate_genotype(allele_one + allele_two)


def genotype_probabilities(parent_one: str, parent_two: str) -> dict[str, float]:
    """Calculate the exact offspring genotype distribution by enumerating gametes."""
    first, second = gamete_probabilities(parent_one), gamete_probabilities(parent_two)
    result = {genotype: 0.0 for genotype in GENOTYPES}
    for (allele_one, probability_one), (allele_two, probability_two) in product(first.items(), second.items()):
        result[combine_gametes(allele_one, allele_two)] += probability_one * probability_two
    return {genotype: probability for genotype, probability in result.items() if probability > 0}


def phenotype_for(genotype: str, mode: InheritanceMode | str) -> str:
    """Map genotype to an explicitly named abstract phenotype for the selected model."""
    genotype = validate_genotype(genotype)
    try:
        mode = InheritanceMode(mode)
    except ValueError as exc:
        raise ValueError(f"Unsupported inheritance mode: {mode}") from exc

    if mode is InheritanceMode.COMPLETE_DOMINANCE:
        return "Dominant phenotype" if "A" in genotype else "Recessive phenotype"
    if mode is InheritanceMode.INCOMPLETE_DOMINANCE:
        return {
            "AA": "A-allele phenotype",
            "Aa": "Intermediate phenotype",
            "aa": "a-allele phenotype",
        }[genotype]
    return {
        "AA": "A-allele phenotype",
        "Aa": "Both A- and a-allele phenotypes expressed",
        "aa": "a-allele phenotype",
    }[genotype]


def phenotype_probabilities(genotype_distribution: Mapping[str, float], mode: InheritanceMode | str) -> dict[str, float]:
    """Aggregate genotype probabilities into phenotype probabilities."""
    result: dict[str, float] = {}
    for genotype, probability in genotype_distribution.items():
        phenotype = phenotype_for(genotype, mode)
        result[phenotype] = result.get(phenotype, 0.0) + float(probability)
    return result


def punnett_square(parent_one: str, parent_two: str) -> list[dict[str, str]]:
    """Return one row per unique gamete pairing, with probability mass."""
    first, second = gamete_probabilities(parent_one), gamete_probabilities(parent_two)
    return [
        {"parent_one_allele": a, "parent_two_allele": b,
         "genotype": combine_gametes(a, b), "probability": pa * pb}
        for (a, pa), (b, pb) in product(first.items(), second.items())
    ]
