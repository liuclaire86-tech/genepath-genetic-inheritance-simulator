import pytest

from genetics import InheritanceMode
from probability import comparison_frame, theoretical_distributions


def test_theoretical_distributions_for_classic_cross():
    genotypes, phenotypes = theoretical_distributions("Aa", "Aa", InheritanceMode.COMPLETE_DOMINANCE)
    assert genotypes == pytest.approx({"AA": .25, "Aa": .5, "aa": .25})
    assert phenotypes == pytest.approx({"Dominant phenotype": .75, "Recessive phenotype": .25})


def test_comparison_frame_includes_missing_observed_categories_and_absolute_error():
    frame = comparison_frame({"AA": .5, "Aa": .5}, {"AA": .6}, "Genotype").set_index("category")
    assert frame.loc["AA", "absolute_error"] == pytest.approx(.1)
    assert frame.loc["Aa", "observed_frequency"] == 0
