"""GenePath interactive Streamlit application."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from genetics import GENOTYPES, InheritanceMode, genotype_probabilities, phenotype_probabilities
from probability import comparison_frame
from simulation import ALLOWED_SAMPLE_SIZES, simulate_offspring
from visualization import comparison_chart, punnett_figure

st.set_page_config(page_title="GenePath | Inheritance Simulator", page_icon="🧬", layout="wide")
st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Manrope:wght@600;700;800&display=swap');
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;color:#203246}.stApp{background:#F5F8FA}.block-container{max-width:1220px;padding-top:2rem;padding-bottom:3rem}h1,h2,h3{font-family:'Manrope',sans-serif;color:#17324D}.hero{background:linear-gradient(120deg,#17324D,#176C72);padding:30px 34px;border-radius:20px;color:white;margin-bottom:20px}.hero h1{color:white;margin:0;font-size:2.35rem}.hero p{color:#DAE9EC;font-size:1.05rem;margin:8px 0 0}.note{background:#E9F3F3;border-left:4px solid #168C8C;border-radius:6px;padding:12px 16px;color:#304958}[data-testid="stMetric"]{background:white;padding:16px;border:1px solid #E4EAEE;border-radius:14px;box-shadow:0 4px 14px #17324d08}
</style>""", unsafe_allow_html=True)
st.markdown('<div class="hero"><h1>GenePath</h1><p>Genetic Inheritance Simulator · explore a simplified single-gene model</p></div>', unsafe_allow_html=True)

with st.sidebar:
    st.markdown("## Model setup")
    parent_one = st.selectbox("Parent 1 genotype", GENOTYPES, index=1, help="Diploid genotype at one hypothetical autosomal locus.")
    parent_two = st.selectbox("Parent 2 genotype", GENOTYPES, index=1, help="Each parent transmits one of their two alleles.")
    mode = st.selectbox("Inheritance mode", list(InheritanceMode), format_func=lambda item: item.value,
                        help="Controls the abstract genotype-to-phenotype mapping.")
    st.divider()
    sample_size = st.selectbox("Hypothetical offspring", ALLOWED_SAMPLE_SIZES, index=1, format_func=lambda n: f"{n:,}")
    seed_value = st.number_input("Random seed", min_value=0, max_value=4_294_967_295, value=42, step=1,
                                 help="Reusing the same seed and settings reproduces the same simulation.")

genotype_theory = genotype_probabilities(parent_one, parent_two)
phenotype_theory = phenotype_probabilities(genotype_theory, mode)
overview, prediction, lab, export, notes = st.tabs(["Overview", "Theoretical predictions", "Simulation lab", "Results & export", "Scientific notes"])

with overview:
    st.subheader("Explore inheritance with transparent assumptions")
    st.write("GenePath calculates exact allele-transmission probabilities, then uses Monte Carlo sampling to show how a finite set of hypothetical offspring can vary around those expectations.")
    a, b, c = st.columns(3)
    a.metric("Locus model", "Single-gene", help="One autosomal locus with two alleles.")
    b.metric("Cross", f"{parent_one} × {parent_two}")
    c.metric("Inheritance rule", mode.value)
    st.markdown('<div class="note">Genotype describes the allele pair (for example, Aa). Phenotype is the abstract modeled outcome mapped from that genotype.</div>', unsafe_allow_html=True)
    st.markdown("#### Workflow")
    st.write("Choose parental genotypes and an inheritance mode in the sidebar. Review exact expectations, then run a reproducible sample to compare observed frequencies.")

with prediction:
    st.subheader("Theoretical predictions")
    st.write(f"Each parent contributes one allele. For **{parent_one} × {parent_two}**, the algorithm enumerates possible gamete pairs and multiplies their probabilities.")
    left, right = st.columns([1, 1.1])
    with left:
        st.markdown("#### Punnett square")
        st.plotly_chart(punnett_figure(parent_one, parent_two), width="stretch", config={"displayModeBar": False})
        st.caption("Rows and columns show allele copies; repeated letters are retained for homozygous parents.")
    with right:
        st.markdown("#### Exact genotype probabilities")
        gt = pd.DataFrame([{"Genotype": key, "Probability": value, "Percent": f"{value:.1%}"} for key, value in genotype_theory.items()])
        st.dataframe(gt, hide_index=True, width="stretch")
        st.markdown("#### Phenotype probabilities")
        pt = pd.DataFrame([{"Abstract phenotype": key, "Probability": value, "Percent": f"{value:.1%}"} for key, value in phenotype_theory.items()])
        st.dataframe(pt, hide_index=True, width="stretch")
    st.caption("Phenotype labels are educational abstractions, not claims about any human trait.")

with lab:
    st.subheader("Monte Carlo simulation")
    st.write(f"Simulate **{sample_size:,}** hypothetical offspring using NumPy. Each parent transmits an allele sampled according to their genotype.")
    if st.button("Run simulation", type="primary"):
        st.session_state["simulation_result"] = simulate_offspring(parent_one, parent_two, int(sample_size), int(seed_value), mode)
        st.session_state["simulation_signature"] = (parent_one, parent_two, mode.value, int(sample_size), int(seed_value))
    signature = (parent_one, parent_two, mode.value, int(sample_size), int(seed_value))
    result = st.session_state.get("simulation_result") if st.session_state.get("simulation_signature") == signature else None
    if result:
        st.caption(f"Seed: {result['seed']} · sample size: {result['sample_size']:,}. A finite random sample will generally differ from exact expectations.")
        g1, g2 = st.columns(2)
        with g1:
            st.markdown("#### Genotypes")
            st.plotly_chart(comparison_chart(genotype_theory, result["genotypes"], "Genotype"), width="stretch")
            frame = comparison_frame(genotype_theory, result["genotypes"], "Genotype").rename(columns={"category":"Genotype", "theoretical_probability":"Theory", "observed_frequency":"Simulated", "absolute_error":"Absolute error"})
            st.dataframe(frame, hide_index=True, width="stretch")
        with g2:
            st.markdown("#### Phenotypes")
            st.plotly_chart(comparison_chart(phenotype_theory, result["phenotypes"], "Phenotype"), width="stretch")
            frame = comparison_frame(phenotype_theory, result["phenotypes"], "Phenotype").rename(columns={"category":"Phenotype", "theoretical_probability":"Theory", "observed_frequency":"Simulated", "absolute_error":"Absolute error"})
            st.dataframe(frame, hide_index=True, width="stretch")
    else:
        st.info("Choose settings and select **Run simulation** to generate results. Changing a setting requires a new run.")

with export:
    st.subheader("Results and export")
    signature = (parent_one, parent_two, mode.value, int(sample_size), int(seed_value))
    result = st.session_state.get("simulation_result") if st.session_state.get("simulation_signature") == signature else None
    if not result:
        st.info("Run a simulation with the current settings to enable CSV downloads.")
    else:
        st.write("Each CSV includes theoretical probability, simulated frequency, absolute error, and the simulation parameters.")
        params = {"parent_one": parent_one, "parent_two": parent_two, "inheritance_mode": mode.value, "sample_size": result["sample_size"], "seed": result["seed"]}
        genotype_df = comparison_frame(genotype_theory, result["genotypes"], "category").rename(columns={"category":"genotype"})
        phenotype_df = comparison_frame(phenotype_theory, result["phenotypes"], "category").rename(columns={"category":"phenotype"})
        for frame in (genotype_df, phenotype_df):
            for key, value in params.items():
                frame[key] = value
        st.download_button("Download genotype distribution CSV", genotype_df.to_csv(index=False).encode("utf-8"), "genepath_genotype_distribution.csv", "text/csv")
        st.download_button("Download phenotype distribution CSV", phenotype_df.to_csv(index=False).encode("utf-8"), "genepath_phenotype_distribution.csv", "text/csv")
        st.markdown("#### Current genotype table")
        st.dataframe(genotype_df, hide_index=True, width="stretch")
        st.markdown("#### Current phenotype table")
        st.dataframe(phenotype_df, hide_index=True, width="stretch")

with notes:
    st.subheader("Scientific notes and model boundaries")
    st.markdown("""
- This is a **simplified single-gene autosomal model** with two hypothetical alleles, `A` and `a`.
- Allele transmission follows the probabilities determined by each parent's genotype; the parents are treated as independent contributors.
- Phenotype mappings assume **complete penetrance** within each selected educational model.
- Complete dominance maps `AA` and `Aa` to a dominant phenotype, and `aa` to a recessive phenotype.
- Incomplete dominance maps `Aa` to an abstract intermediate phenotype. Codominance maps `Aa` to an abstract outcome where both allele-associated phenotypes are expressed.
- These labels do not describe established human traits. The application does not model environmental influences, most polygenic traits, gene interactions, sex-linked inheritance, or the full complexity of human genetics.
- Monte Carlo outcomes are random samples. Larger samples tend to be closer to theoretical probabilities, but simulated values are not guarantees about real offspring.
""")
    st.caption("For education only. GenePath is not a clinical or diagnostic tool.")
