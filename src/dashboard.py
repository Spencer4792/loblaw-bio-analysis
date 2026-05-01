"""Streamlit dashboard. Run: streamlit run src/dashboard.py"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from analysis import (
    CELL_POPULATIONS,
    baseline_breakdown,
    connect,
    load_wide,
    melanoma_male_responder_b_cell_mean,
    responder_frequencies,
    responder_stats,
    summary_table,
)

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "loblaw_bio.db"


st.set_page_config(
    page_title="Loblaw Bio Cell Count Dashboard",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("Loblaw Bio Cell Count Analytics")
st.caption("Interactive view of the miraclib clinical-trial cell-count data.")

if not DB_PATH.exists():
    csv_path = ROOT / "cell-count.csv"
    if csv_path.exists():
        import sys
        sys.path.insert(0, str(ROOT))
        from load_data import load as _load_db
        _load_db(csv_path, DB_PATH)
    else:
        st.error(
            f"Database not found at `{DB_PATH.name}` and cell-count.csv missing. "
            "Place cell-count.csv at the repo root and reload."
        )
        st.stop()


@st.cache_data(show_spinner=False)
def _load_all() -> dict[str, pd.DataFrame]:
    with connect(DB_PATH) as conn:
        return {
            "wide": load_wide(conn),
            "summary": summary_table(conn),
            "responder_long": responder_frequencies(conn),
            "responder_stats": responder_stats(conn),
            "baseline": baseline_breakdown(conn),
            "b_cell_mean": melanoma_male_responder_b_cell_mean(conn),
        }


data = _load_all()
wide = data["wide"]


st.sidebar.header("Filters")

def _multi(label, options, default=None):
    options = sorted([o for o in options if pd.notna(o)])
    return st.sidebar.multiselect(label, options, default=default or options)

projects = _multi("Project", wide["project_id"].dropna().unique().tolist())
conditions = _multi("Condition", wide["condition"].dropna().unique().tolist())
treatments = _multi("Treatment", wide["treatment"].dropna().unique().tolist())
sample_typ = _multi("Sample type", wide["sample_type"].dropna().unique().tolist())
responses = _multi("Response", wide["response"].dropna().unique().tolist())
sexes = _multi("Sex", wide["sex"].dropna().unique().tolist())

filtered = wide[
    wide["project_id"].isin(projects)
    & wide["condition"].isin(conditions)
    & wide["treatment"].isin(treatments)
    & wide["sample_type"].isin(sample_typ)
    & wide["response"].isin(responses + [None])  # keep nulls if user keeps them
    & wide["sex"].isin(sexes + [None])
]

st.sidebar.markdown(
    f"**{len(filtered):,}** samples / **{filtered['subject_id'].nunique():,}** subjects "
    "after filters."
)


tab_overview, tab_summary, tab_response, tab_baseline = st.tabs(
    ["Overview", "Part 2 · Summary table",
     "Part 3 · Responder analysis", "Part 4 · Baseline subset"]
)


with tab_overview:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples",  f"{len(filtered):,}")
    c2.metric("Subjects", f"{filtered['subject_id'].nunique():,}")
    c3.metric("Projects", f"{filtered['project_id'].nunique():,}")
    c4.metric("Conditions", f"{filtered['condition'].nunique():,}")

    st.subheader("Sample table")
    show_cols = [
        "sample_id", "subject_id", "project_id", "condition", "sex",
        "treatment", "response", "sample_type", "time_from_treatment_start",
        *CELL_POPULATIONS, "total_count",
    ]
    st.dataframe(
        filtered[show_cols].sort_values(["project_id", "subject_id", "sample_id"]),
        use_container_width=True,
        hide_index=True,
    )


with tab_summary:
    st.subheader("Relative frequency of each cell population (Part 2)")

    summary = data["summary"]
    summary_filtered = summary[summary["sample"].isin(filtered["sample_id"])]
    st.dataframe(
        summary_filtered,
        use_container_width=True,
        hide_index=True,
        column_config={
            "count": st.column_config.NumberColumn(format="%d"),
            "total_count": st.column_config.NumberColumn(format="%d"),
            "percentage": st.column_config.NumberColumn(format="%.2f%%"),
        },
    )

    st.download_button(
        "Download summary as CSV",
        summary_filtered.to_csv(index=False).encode(),
        file_name="summary_table.csv",
        mime="text/csv",
    )

    st.subheader("Distribution by population")
    fig = px.box(
        summary_filtered,
        x="population", y="percentage",
        points="all", color="population",
        labels={"percentage": "Relative frequency (%)"},
    )
    fig.update_layout(showlegend=False, height=450)
    st.plotly_chart(fig, use_container_width=True)


with tab_response:
    st.subheader("Responders vs non-responders, melanoma PBMC samples on miraclib (Part 3)")

    long = data["responder_long"]
    if long.empty:
        st.info("No melanoma / miraclib / PBMC samples with response yes/no in this database.")
    else:
        c1, c2 = st.columns([2, 1])
        with c1:
            fig = px.box(
                long,
                x="population", y="percentage", color="response",
                points="all",
                category_orders={"response": ["yes", "no"]},
                color_discrete_map={"yes": "steelblue", "no": "crimson"},
                labels={"percentage": "Relative frequency (%)",
                        "response": "Response"},
            )
            fig.update_layout(height=480, boxmode="group")
            st.plotly_chart(fig, use_container_width=True)
        with c2:
            n_yes = long.loc[long["response"] == "yes", "sample"].nunique()
            n_no  = long.loc[long["response"] == "no",  "sample"].nunique()
            st.metric("Responder samples", n_yes)
            st.metric("Non-responder samples", n_no)

        st.subheader("Statistical comparison")
        st.caption("Mann-Whitney U, two-sided. `q_value_bh` is the BH-FDR adjusted p-value across the five populations.")
        stats_df = data["responder_stats"]
        st.dataframe(
            stats_df.style.format({
                "median_responders": "{:.2f}",
                "median_non_responders": "{:.2f}",
                "mean_responders": "{:.2f}",
                "mean_non_responders": "{:.2f}",
                "mannwhitney_u": "{:.1f}",
                "p_value": "{:.4f}",
                "q_value_bh": "{:.4f}",
                "rank_biserial_effect_size": "{:.3f}",
            }),
            use_container_width=True,
            hide_index=True,
        )

        sig = stats_df[stats_df["p_value"] < 0.05]["population"].tolist()
        if sig:
            st.success(
                "Populations with a significant difference (p < 0.05): "
                + ", ".join(sig)
            )
        else:
            st.warning("No populations reached p < 0.05 in this dataset.")


with tab_baseline:
    st.subheader("Melanoma PBMC samples at baseline on miraclib (Part 4)")
    bd = data["baseline"]
    sub = bd["subset"]
    st.markdown(f"**{len(sub)}** samples · **{sub['subject_id'].nunique()}** subjects")

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**By project**")
        st.dataframe(bd["by_project"], use_container_width=True, hide_index=True)
    with c2:
        st.markdown("**By response (subjects)**")
        st.dataframe(bd["by_response"], use_container_width=True, hide_index=True)
    with c3:
        st.markdown("**By sex (subjects)**")
        st.dataframe(bd["by_sex"], use_container_width=True, hide_index=True)

    st.subheader("Average B cells for melanoma males, responders, t=0")
    avg = data["b_cell_mean"]
    if avg is None:
        st.info("No qualifying samples in this database.")
    else:
        st.metric("Mean B-cell count", f"{avg:.2f}")

    st.markdown("**Underlying samples**")
    st.dataframe(sub, use_container_width=True, hide_index=True)
