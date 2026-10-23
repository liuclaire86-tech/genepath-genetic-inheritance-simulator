import pytest

from genetics import (
    InheritanceMode, combine_gametes, gamete_probabilities,
    genotype_probabilities, phenotype_for, phenotype_probabilities,
    validate_genotype,
)


@pytest.mark.parametrize(("raw", "expected"), [("AA", "AA"), ("Aa", "Aa"), ("aA", "Aa"), ("aa", "aa")])
def test_validate_and_canonicalize_genotype(raw, expected):
    assert validate_genotype(raw) == expected


@pytest.mark.parametrize("invalid", ["", "A", "AAA", "AB", "bb", None, 12])
def test_invalid_genotypes_raise(invalid):
    with pytest.raises(ValueError):
        validate_genotype(invalid)


def test_gamete_probabilities():
    assert gamete_probabilities("AA") == {"A": 1.0}
    assert gamete_probabilities("Aa") == {"A": 0.5, "a": 0.5}
    assert gamete_probabilities("aa") == {"a": 1.0}


def test_combine_gametes_and_invalid_alleles():
    assert combine_gametes("a", "A") == "Aa"
    with pytest.raises(ValueError):
        combine_gametes("B", "a")


@pytest.mark.parametrize(("parents", "expected"), [
    (("Aa", "Aa"), {"AA": .25, "Aa": .5, "aa": .25}),
    (("AA", "aa"), {"Aa": 1.0}),
    (("Aa", "aa"), {"Aa": .5, "aa": .5}),
    (("AA", "Aa"), {"AA": .5, "Aa": .5}),
    (("aa", "AA"), {"Aa": 1.0}),
])
def test_exact_genotype_probabilities(parents, expected):
    actual = genotype_probabilities(*parents)
    assert actual == pytest.approx(expected)
    assert sum(actual.values()) == pytest.approx(1.0)


@pytest.mark.parametrize(("mode", "heterozygote"), [
    (InheritanceMode.COMPLETE_DOMINANCE, "Dominant phenotype"),
    (InheritanceMode.INCOMPLETE_DOMINANCE, "Intermediate phenotype"),
    (InheritanceMode.CODOMINANCE, "Both A- and a-allele phenotypes expressed"),
])
def test_heterozygote_phenotype_mapping(mode, heterozygote):
    assert phenotype_for("Aa", mode) == heterozygote


def test_complete_dominance_distribution_for_aaxaa():
    genotypes = genotype_probabilities("Aa", "Aa")
    assert phenotype_probabilities(genotypes, InheritanceMode.COMPLETE_DOMINANCE) == pytest.approx(
        {"Dominant phenotype": .75, "Recessive phenotype": .25}
    )


@pytest.mark.parametrize("mode", list(InheritanceMode))
def test_phenotype_probabilities_sum_to_one(mode):
    distribution = phenotype_probabilities(genotype_probabilities("Aa", "Aa"), mode)
    assert sum(distribution.values()) == pytest.approx(1.0)


def test_invalid_inheritance_mode_raises():
    with pytest.raises(ValueError):
        phenotype_for("Aa", "polygenic")
