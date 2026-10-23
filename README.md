# GenePath — Genetic Inheritance Simulator

GenePath is an interactive educational application for exploring how allele transmission creates genotype and phenotype probability distributions in a simplified single-gene autosomal model. Configure hypothetical parents, inspect exact predictions, simulate offspring with NumPy, and compare the sample with the theoretical expectations.

> **Educational model only.** GenePath does not predict real human traits or provide clinical or diagnostic guidance.

## Problem and objectives

Inheritance lessons often move between Punnett squares, probability, and observed outcomes without showing how these ideas connect computationally. GenePath makes the assumptions explicit and lets learners compare exact probability with a reproducible Monte Carlo sample. The project demonstrates modular Python design, scientific modeling, numerical simulation, data analysis, visualization, testing, and documentation.

## Features

- Validate and canonicalize `AA`, `Aa`, and `aa` diploid genotypes.
- Calculate gamete probabilities and enumerate exact offspring genotype probabilities for any supported cross.
- Model complete dominance, incomplete dominance, and codominance with distinct, documented abstract phenotype mappings.
- Render a Punnett square and exact genotype and phenotype distributions.
- Simulate 100, 1,000, 10,000, or 100,000 hypothetical offspring using vectorized NumPy sampling and a configurable seed.
- Compare theoretical and observed frequencies with absolute errors in interactive Plotly charts and tables.
- Export genotype and phenotype comparison tables as CSV with simulation parameters.

## Technology stack

Python 3.11+, Streamlit, NumPy, Pandas, Plotly, and Pytest.

## Architecture

```text
app.py             Streamlit interface, controls, results, and downloads
genetics.py        Genotype validation, gametes, crosses, phenotype rules
probability.py     Exact distributions and tidy comparison data
simulation.py      Vectorized sampling, validation, and frequency aggregation
visualization.py   Plotly comparison charts and Punnett square
tests/             Unit, edge-case, reproducibility, and statistical tests
```

The interface calls independent scientific modules. The exact probability algorithm enumerates parental gamete combinations and multiplies their probabilities; it does not special-case any cross. The simulation samples one allele from each parent for every offspring, canonicalizes the resulting genotype, maps phenotypes under the selected mode, and aggregates frequencies. A fixed seed makes a run reproducible.

## Install and run on macOS

Requires Python 3.11 or newer.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Streamlit prints a local URL, normally `http://localhost:8501`.

## Run tests

```bash
python -m pytest -v
```

## Sample scientific experiment

Set both parents to `Aa`, choose complete dominance, and calculate the theoretical predictions. The algorithm gives `AA = 25%`, `Aa = 50%`, and `aa = 25%`; under the complete-dominance mapping that becomes dominant phenotype `75%` and recessive phenotype `25%`. Run a sample of 10,000 with a chosen seed. The simulated frequencies should usually be near the exact values, but will vary because each offspring is a random draw. Reuse the same seed and settings to reproduce the same result.

## Model assumptions and scientific limitations

- The initial model represents one autosomal gene with two hypothetical alleles, `A` and `a`.
- Each parent contributes one allele with probabilities determined by its genotype; transmission between parents is modeled independently.
- Genotype is the allele pair. Phenotype is a separate, simplified mapping from genotype.
- Phenotype mappings assume complete penetrance within the selected educational model.
- Complete dominance: `AA` and `Aa` map to “Dominant phenotype”; `aa` maps to “Recessive phenotype.”
- Incomplete dominance: `AA` and `aa` map to their respective allele-associated abstract phenotypes; `Aa` maps to “Intermediate phenotype.”
- Codominance: homozygotes map to their respective allele-associated abstract phenotypes; `Aa` maps to “Both A- and a-allele phenotypes expressed.”
- Phenotype descriptions are hypothetical model labels, not claims about established human traits.
- The model does not represent environmental influences, most polygenic traits, gene interactions, sex-linked inheritance, linkage, mutation, variable penetrance, or the full complexity of human genetics.
- A Monte Carlo sample is probabilistic. A simulated outcome is not a guarantee about real offspring, and this application is not a clinical tool.

## Screenshots

*Screenshot placeholder: add a screenshot of the theoretical predictions dashboard after running the app.*

*Screenshot placeholder: add a screenshot of the simulation comparison and CSV export section after running the app.*

## Future improvements

- Add confidence intervals and goodness-of-fit summaries for simulated versus theoretical frequencies.
- Add optional downloadable experiment reports and support for multiple independently run experiments.
- Add further inheritance models only with clearly defined assumptions and dedicated tests.
- Improve accessibility, mobile layout, and localization.
- Add continuous integration for linting and test execution.

## License

No license has been added to this project. Until a license is selected and included, reuse permissions are not explicitly granted.
