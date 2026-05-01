# Cell Count Analysis for Teiko Technical Exam

Analytics tool and Streamlit dashboard for the miraclib clinical-trial cell data.

## Running it

cell-count.csv needs to be in the project root, then:

```
make setup
make pipeline
make dashboard
```

Dashboard at: http://localhost:8501

Codespaces should work the same, just follow the commands and the dashboard will be active.

## Schema

Four tables:

- `projects(project_id)`
- `subjects(subject_id, project_id, condition, age, sex, treatment, response)`
- `samples(sample_id, subject_id, sample_type, time_from_treatment_start)`
- `cell_counts(sample_id, population, count)`

Subject-level info (sex, treatment, response) shows once in `subjects` rather than being repeated on every sample row. `cell_counts` is in long-format instead of wide so adding a new population is a data-only change, not a total schema migration. Additionally, indexes are on the columns we actually filter on (so things like condition, treatment, response, sample_type + time).

At scale, SQLite handles 10k samples easy, and the schema is super portable to Postgres or DuckDB without any significant changes, you would just need to swap the connection. New metadata columns (things like batch, cohort, etc.) are just one `ALTER TABLE` away.

## Code structure

- `load_data.py` at the root because it's required there
- `src/analysis.py` holds functions that take a SQLite connection and return DataFrames
- `src/dashboard.py` imports from analysis.py so the dashboard and pipeline don't drift apart

## Dashboard

Runs at http://localhost:8501 after the `make dashboard` command.
