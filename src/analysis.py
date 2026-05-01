"""Analysis functions for parts 2-4 plus the pipeline runner."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "loblaw_bio.db"

CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def connect(db_path: Path | str = DEFAULT_DB) -> sqlite3.Connection:
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(
            f"Database not found at {db_path}. Run `python load_data.py` first."
        )
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def load_wide(conn: sqlite3.Connection) -> pd.DataFrame:
    """One row per sample with cell counts and metadata."""
    counts = pd.read_sql_query("SELECT * FROM cell_counts", conn)
    wide = counts.pivot(index="sample_id", columns="population", values="count")
    for pop in CELL_POPULATIONS:
        if pop not in wide.columns:
            wide[pop] = 0
    wide = wide[CELL_POPULATIONS].fillna(0).astype(int).reset_index()

    samples = pd.read_sql_query("SELECT * FROM samples", conn)
    subjects = pd.read_sql_query("SELECT * FROM subjects", conn)

    df = (
        samples
        .merge(subjects, on="subject_id", how="left")
        .merge(wide, on="sample_id", how="left")
    )
    df["total_count"] = df[CELL_POPULATIONS].sum(axis=1)
    return df


def summary_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 2: per-sample, per-population relative frequencies."""
    wide = load_wide(conn)[["sample_id", "total_count", *CELL_POPULATIONS]]
    long = wide.melt(
        id_vars=["sample_id", "total_count"],
        value_vars=CELL_POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long["percentage"] = np.where(
        long["total_count"] > 0,
        100.0 * long["count"] / long["total_count"],
        0.0,
    )
    long = long.rename(columns={"sample_id": "sample"})
    return long[["sample", "total_count", "population", "count", "percentage"]]


def melanoma_miraclib_pbmc(conn: sqlite3.Connection) -> pd.DataFrame:
    wide = load_wide(conn)
    keep = (
        wide["condition"].str.lower().eq("melanoma")
        & wide["treatment"].str.lower().eq("miraclib")
        & wide["sample_type"].str.upper().eq("PBMC")
    )
    return wide.loc[keep].copy()


def responder_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 3: long-format frequencies for melanoma/miraclib/PBMC samples with yes/no response."""
    df = melanoma_miraclib_pbmc(conn)
    df = df[df["response"].isin(["yes", "no"])].copy()
    if df.empty:
        return pd.DataFrame(columns=[
            "sample", "subject_id", "response", "population", "percentage",
        ])

    for pop in CELL_POPULATIONS:
        df[f"{pop}_pct"] = np.where(
            df["total_count"] > 0,
            100.0 * df[pop] / df["total_count"],
            0.0,
        )
    long = df.melt(
        id_vars=["sample_id", "subject_id", "response"],
        value_vars=[f"{p}_pct" for p in CELL_POPULATIONS],
        var_name="population",
        value_name="percentage",
    )
    long["population"] = long["population"].str.replace("_pct$", "", regex=True)
    return long.rename(columns={"sample_id": "sample"})


def _bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(list(pvals), dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def responder_stats(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 3: Mann-Whitney U per population, with BH-FDR q-values."""
    long = responder_frequencies(conn)
    if long.empty:
        return pd.DataFrame()

    rows = []
    for pop, grp in long.groupby("population"):
        yes = grp.loc[grp["response"] == "yes", "percentage"].to_numpy()
        no  = grp.loc[grp["response"] == "no",  "percentage"].to_numpy()
        if len(yes) == 0 or len(no) == 0:
            rows.append({
                "population": pop, "n_responders": len(yes),
                "n_non_responders": len(no),
                "median_responders": np.nan, "median_non_responders": np.nan,
                "mean_responders": np.nan, "mean_non_responders": np.nan,
                "mannwhitney_u": np.nan, "p_value": np.nan,
                "rank_biserial_effect_size": np.nan,
            })
            continue
        u, p = stats.mannwhitneyu(yes, no, alternative="two-sided")
        rb = 1 - (2 * u) / (len(yes) * len(no))  # rank-biserial effect size
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responders": float(np.median(yes)),
            "median_non_responders": float(np.median(no)),
            "mean_responders": float(np.mean(yes)),
            "mean_non_responders": float(np.mean(no)),
            "mannwhitney_u": float(u),
            "p_value": float(p),
            "rank_biserial_effect_size": float(rb),
        })

    out = pd.DataFrame(rows).sort_values("p_value", na_position="last")
    out["q_value_bh"] = _bh_fdr(out["p_value"].fillna(1.0).to_numpy())
    out["significant_p05"] = out["p_value"].lt(0.05)
    return out.reset_index(drop=True)


def baseline_subset(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 4: melanoma PBMC samples at baseline (time=0) on miraclib."""
    sql = """
        SELECT s.sample_id, s.subject_id, s.sample_type,
               s.time_from_treatment_start,
               sub.project_id, sub.condition, sub.sex,
               sub.treatment, sub.response
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE LOWER(sub.condition)  = 'melanoma'
          AND LOWER(sub.treatment)  = 'miraclib'
          AND UPPER(s.sample_type)  = 'PBMC'
          AND s.time_from_treatment_start = 0
    """
    return pd.read_sql_query(sql, conn)


def baseline_breakdown(conn: sqlite3.Connection) -> dict[str, pd.DataFrame]:
    sub = baseline_subset(conn)

    by_project = (
        sub.groupby("project_id")["sample_id"].nunique()
           .rename("n_samples").reset_index()
           .sort_values("project_id").reset_index(drop=True)
    )

    # response/sex breakdowns are by subject, not sample
    subj = sub.drop_duplicates("subject_id")

    by_response = (
        subj["response"].fillna("unknown").value_counts()
            .rename_axis("response").reset_index(name="n_subjects")
    )
    by_sex = (
        subj["sex"].fillna("unknown").value_counts()
            .rename_axis("sex").reset_index(name="n_subjects")
    )

    return {
        "subset": sub,
        "by_project": by_project,
        "by_response": by_response,
        "by_sex": by_sex,
    }


def melanoma_male_responder_b_cell_mean(conn: sqlite3.Connection) -> float | None:
    """Mean B-cell count for melanoma male responders at time=0."""
    sql = """
        SELECT cc.count AS b_cell_count
        FROM samples s
        JOIN subjects sub  ON sub.subject_id = s.subject_id
        JOIN cell_counts cc ON cc.sample_id = s.sample_id
        WHERE LOWER(sub.condition) = 'melanoma'
          AND UPPER(sub.sex) = 'M'
          AND LOWER(sub.response) = 'yes'
          AND LOWER(sub.treatment) = 'miraclib'
          AND UPPER(s.sample_type) = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND cc.population = 'b_cell'
    """
    df = pd.read_sql_query(sql, conn)
    if df.empty:
        return None
    return round(float(df["b_cell_count"].mean()), 2)


def make_boxplot(long_df: pd.DataFrame, out_path: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")  # headless
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    pops = sorted(long_df["population"].unique())
    width = 0.35
    positions = np.arange(len(pops))

    for i, resp in enumerate(["yes", "no"]):
        data = [
            long_df.loc[
                (long_df["population"] == pop) & (long_df["response"] == resp),
                "percentage",
            ].to_numpy()
            for pop in pops
        ]
        offset = (i - 0.5) * width
        bp = ax.boxplot(
            data,
            positions=positions + offset,
            widths=width * 0.9,
            patch_artist=True,
            showfliers=True,
        )
        color = "tab:blue" if resp == "yes" else "tab:red"
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
        for med in bp["medians"]:
            med.set_color("black")

    ax.set_xticks(positions)
    ax.set_xticklabels(pops)
    ax.set_ylabel("Relative frequency (%)")
    ax.set_title(
        "Cell-population frequencies, melanoma PBMC samples on miraclib\n"
        "responders (blue) vs non-responders (red)"
    )
    legend = [
        plt.Rectangle((0, 0), 1, 1, color="tab:blue", alpha=0.55, label="responders"),
        plt.Rectangle((0, 0), 1, 1, color="tab:red", alpha=0.55, label="non-responders"),
    ]
    ax.legend(handles=legend, loc="upper right")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def run_pipeline(db_path: Path = DEFAULT_DB,
                 out_dir: Path = ROOT / "outputs") -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with connect(db_path) as conn:
        # Part 2
        summary = summary_table(conn)
        summary.to_csv(out_dir / "summary_table.csv", index=False)

        # Part 3
        long = responder_frequencies(conn)
        long.to_csv(out_dir / "responder_frequencies.csv", index=False)
        stats_df = responder_stats(conn)
        stats_df.to_csv(out_dir / "responder_stats.csv", index=False)
        if not long.empty:
            make_boxplot(long, out_dir / "responder_boxplot.png")

        # Part 4
        breakdown = baseline_breakdown(conn)
        breakdown["subset"].to_csv(out_dir / "baseline_subset.csv", index=False)
        breakdown["by_project"].to_csv(out_dir / "baseline_by_project.csv", index=False)
        breakdown["by_response"].to_csv(out_dir / "baseline_by_response.csv", index=False)
        breakdown["by_sex"].to_csv(out_dir / "baseline_by_sex.csv", index=False)

        avg_b = melanoma_male_responder_b_cell_mean(conn)
        msg = (f"{avg_b:.2f}" if avg_b is not None else "n/a (no qualifying samples)")
        (out_dir / "melanoma_male_responder_b_cell_mean.txt").write_text(msg + "\n")

    print("Pipeline complete. Outputs written to:", out_dir)
    print(f"  Average B cells (melanoma males, responders, t=0): {msg}")


if __name__ == "__main__":
    run_pipeline()
