.PHONY: setup pipeline dashboard clean

PYTHON ?= python3
DB := loblaw_bio.db
CSV := cell-count.csv

# ---------------------------------------------------------------------------
# Install dependencies
# ---------------------------------------------------------------------------
setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

# ---------------------------------------------------------------------------
# Run the full pipeline: load CSV -> SQLite, then generate all outputs
# ---------------------------------------------------------------------------
pipeline:
	@if [ ! -f $(CSV) ]; then \
		echo "ERROR: $(CSV) not found in repository root."; \
		echo "Place cell-count.csv next to this Makefile and re-run."; \
		exit 1; \
	fi
	$(PYTHON) load_data.py
	$(PYTHON) -m src.analysis

# ---------------------------------------------------------------------------
# Launch the Streamlit dashboard
# ---------------------------------------------------------------------------
dashboard:
	@if [ ! -f $(DB) ]; then \
		echo "Database $(DB) not found — running pipeline first..."; \
		$(MAKE) pipeline; \
	fi
	streamlit run src/dashboard.py --server.port 8501 --server.address 0.0.0.0

# ---------------------------------------------------------------------------
clean:
	rm -f $(DB)
	rm -rf outputs/*.csv outputs/*.png outputs/*.txt
