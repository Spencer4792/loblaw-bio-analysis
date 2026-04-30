"""
Build the SQLite database for the cell-count analysis.

Reads cell-count.csv from the repository root and writes loblaw_bio.db
with four tables: projects, subjects, samples, cell_counts.

    python load_data.py
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_CSV = ROOT / "cell-count.csv"
DEFAULT_DB = ROOT / "loblaw_bio.db"

CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT,
    age        INTEGER,
    sex        TEXT,
    treatment  TEXT,
    response   TEXT
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT,
    time_from_treatment_start INTEGER  -- days from treatment start; 0 = baseline
);

CREATE TABLE cell_counts (
    sample_id   TEXT NOT NULL REFERENCES samples(sample_id),
    population  TEXT NOT NULL,
    count       INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_project   ON subjects(project_id);
CREATE INDEX idx_subjects_condition ON subjects(condition);
CREATE INDEX idx_subjects_treatment ON subjects(treatment);
CREATE INDEX idx_subjects_response  ON subjects(response);
CREATE INDEX idx_samples_subject    ON samples(subject_id);
CREATE INDEX idx_samples_type_time  ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_cell_counts_pop    ON cell_counts(population);
"""

# just in case
COLUMN_ALIASES = {
    "project": ["project", "project_id"],
    "subject": ["subject", "subject_id", "patient", "patient_id"],
    "condition": ["condition", "indication", "disease"],
    "age": ["age"],
    "sex": ["sex", "gender"],
    "treatment": ["treatment", "drug"],
    "response": ["response"],
    "sample": ["sample", "sample_id"],
    "sample_type": ["sample_type", "tissue"],
    "time_from_treatment_start": ["time_from_treatment_start", "time", "timepoint"],
}


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={c: c.strip().lower() for c in df.columns})
    rename = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in df.columns and alias != canonical:
                rename[alias] = canonical
                break
    return df.rename(columns=rename)


def _validate(df: pd.DataFrame) -> None:
    required = {
        "project", "subject", "condition", "sex", "treatment",
        "sample", "sample_type", "time_from_treatment_start",
    } | set(CELL_POPULATIONS)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "cell-count.csv is missing expected columns: "
            f"{sorted(missing)}. Found columns: {sorted(df.columns)}"
        )

def load(csv_path: Path = DEFAULT_CSV, db_path: Path = DEFAULT_DB) -> Path:
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Could not find {csv_path}. Place cell-count.csv in the "
            "repository root (or pass --csv)."
        )

    df = pd.read_csv(csv_path)
    df = _normalize_columns(df)
    _validate(df)

    for pop in CELL_POPULATIONS:
        df[pop] = pd.to_numeric(df[pop], errors="coerce").fillna(0).astype(int)
    df["time_from_treatment_start"] = pd.to_numeric(
        df["time_from_treatment_start"], errors="coerce"
    )
    if "age" in df.columns:
        df["age"] = pd.to_numeric(df["age"], errors="coerce")

    str_cols = ["project", "subject", "condition", "sex", "treatment",
                "response", "sample", "sample_type"]
    for c in str_cols:
        if c in df.columns:
            df[c] = df[c].astype("string").str.strip()
            df.loc[df[c] == "", c] = pd.NA

    projects = (
        df[["project"]].drop_duplicates().rename(columns={"project": "project_id"})
    )

    subject_cols = ["subject", "project", "condition", "sex",
                    "treatment", "response"]
    if "age" in df.columns:
        subject_cols.insert(3, "age")
    subjects = (
        df[subject_cols]
        .drop_duplicates(subset=["subject"])
        .rename(columns={"subject": "subject_id", "project": "project_id"})
    )

    samples = (
        df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
        .drop_duplicates(subset=["sample"])
        .rename(columns={"sample": "sample_id", "subject": "subject_id"})
    )

    long = df.melt(
        id_vars=["sample"],
        value_vars=CELL_POPULATIONS,
        var_name="population",
        value_name="count",
    ).rename(columns={"sample": "sample_id"})

    if db_path.exists():
        db_path.unlink()

    with sqlite3.connect(db_path) as conn:
        conn.executescript(SCHEMA)
        projects.to_sql("projects", conn, if_exists="append", index=False)
        subjects.to_sql("subjects", conn, if_exists="append", index=False)
        samples.to_sql("samples", conn, if_exists="append", index=False)
        long.to_sql("cell_counts", conn, if_exists="append", index=False)
        conn.commit()

        counts = {
            t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
            for t in ("projects", "subjects", "samples", "cell_counts")
        }

    print(f"Loaded {csv_path.name} into {db_path.name}")
    for table, n in counts.items():
        print(f"  {table:<12s} {n:>6d} rows")
    return db_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    args = parser.parse_args()
    try:
        load(args.csv, args.db)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
