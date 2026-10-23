import pytest

from genetics import InheritanceMode, genotype_probabilities, phenotype_probabilities
from simulation import ALLOWED_SAMPLE_SIZES, simulate_offspring


def test_seed_reproduces_same_simulation():
    first = simulate_offspring("Aa", "Aa", 1_000, 1729, InheritanceMode.CODOMINANCE)
    second = simulate_offspring("Aa", "Aa", 1_000, 1729, InheritanceMode.CODOMINANCE)
    assert first == second


@pytest.mark.parametrize("size", [0, 99, 101, 1_000_000, True, 1.5])
def test_invalid_sample_sizes_raise(size):
    with pytest.raises(ValueError):
        simulate_offspring("Aa", "Aa", size)


@pytest.mark.parametrize("seed", [-1, 4_294_967_296, 1.2, True])
def test_invalid_seeds_raise(seed):
    with pytest.raises(ValueError):
        simulate_offspring("Aa", "Aa", 100, seed)


def test_invalid_parent_genotype_raises():
    with pytest.raises(ValueError):
        simulate_offspring("AB", "Aa", 100)


def test_simulated_distribution_approaches_theory_at_large_sample():
    result = simulate_offspring("Aa", "Aa", 100_000, 2026, InheritanceMode.COMPLETE_DOMINANCE)
    expected = genotype_probabilities("Aa", "Aa")
    for genotype, probability in expected.items():
        assert result["genotypes"][genotype] == pytest.approx(probability, abs=0.01)
    expected_phenotypes = phenotype_probabilities(expected, InheritanceMode.COMPLETE_DOMINANCE)
    for phenotype, probability in expected_phenotypes.items():
        assert result["phenotypes"][phenotype] == pytest.approx(probability, abs=0.01)


def test_sample_size_options_are_supported():
    assert ALLOWED_SAMPLE_SIZES == (100, 1_000, 10_000, 100_000)
