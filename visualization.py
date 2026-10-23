"""Plotly chart builders used by the Streamlit dashboard."""

import plotly.graph_objects as go

from probability import comparison_frame

NAVY, TEAL, GRAY = "#17324D", "#168C8C", "#9AA8B5"


def comparison_chart(theoretical: dict[str, float], observed: dict[str, float], label: str) -> go.Figure:
    frame = comparison_frame(theoretical, observed, label)
    figure = go.Figure()
    figure.add_bar(name="Theoretical", x=frame["category"], y=frame["theoretical_probability"], marker_color=NAVY)
    figure.add_bar(name="Simulated", x=frame["category"], y=frame["observed_frequency"], marker_color=TEAL)
    figure.update_layout(barmode="group", yaxis_title="Probability / frequency", yaxis_tickformat=".0%",
                         legend_title_text="Distribution", margin=dict(l=10, r=10, t=25, b=10),
                         plot_bgcolor="white", paper_bgcolor="white")
    return figure


def punnett_figure(parent_one: str, parent_two: str) -> go.Figure:
    """Render a compact 2x2 Punnett square including repeated homozygous gametes."""
    first, second = list(parent_one), list(parent_two)
    rows = [["" for _ in range(3)]]
    for allele_one in first:
        row = [allele_one]
        for allele_two in second:
            pair = allele_one + allele_two
            row.append("AA" if pair == "AA" else "aa" if pair == "aa" else "Aa")
        rows.append(row)
    columns = list(map(list, zip(*rows)))
    fig = go.Figure(data=[go.Table(
        header=dict(values=["Parent 1 ↓ / Parent 2 →", *second], fill_color=NAVY,
                    font=dict(color="white", size=14), align="center"),
        cells=dict(values=columns, fill_color="#EAF0F5", font=dict(color=NAVY, size=17), align="center", height=38),
    )])
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0), height=190)
    return fig
