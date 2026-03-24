PYTHON ?= python3
VENV_DIR ?= env
MAKEFLAGS += --no-print-directory

RUN_CONFIG ?= configs/run_params.yaml
REF_FRAMES ?= hs3 nnr sa
FORMULATIONS ?= 1 2
DP_REF_VALUES ?= 0 2.35e7
SUMMARY_SUITES ?= best-fit

.PHONY: venv install run-matrix run-matrix-maps matrix-summary run-matrix-with-summary run-matrix-maps-with-summary

# ---------------------------------------------------------------------------
# Setup (one-time)
# ---------------------------------------------------------------------------

# Create the Python virtual environment.
venv:
	$(PYTHON) -m venv $(VENV_DIR)

# Install dependencies into the virtual environment.
install:
	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install -r requirements.txt

# ---------------------------------------------------------------------------
# Building blocks (can be run independently)
# ---------------------------------------------------------------------------

# Run the parameter sweep for all formulations, reference frames, and DP values,
# skipping map generation. Fast — use this for iterating on parameters.
# Override with e.g. FORMULATIONS=1 or REF_FRAMES=nnr or DP_REF_VALUES=0 to run a subset.
run-matrix:
	@for form in $(FORMULATIONS); do \
		for ref in $(REF_FRAMES); do \
			for dp in $(DP_REF_VALUES); do \
				echo "[matrix] formulation=$$form ref=$$ref dp_ref=$$dp"; \
				$(PYTHON) compute_rates_misfit.py --config $(RUN_CONFIG) --formulation $$form --vt-ref $$ref --dp-ref $$dp --skip-map --out-prefix plots/$$ref || exit $$?; \
			done \
		done \
	done

# Run the parameter sweep for all formulations, reference frames, and DP values,
# including full global map generation. Slow — use this for final/publication outputs.
run-matrix-maps:
	@for form in $(FORMULATIONS); do \
		for ref in $(REF_FRAMES); do \
			for dp in $(DP_REF_VALUES); do \
				echo "[matrix-maps] formulation=$$form ref=$$ref dp_ref=$$dp"; \
				$(PYTHON) compute_rates_misfit.py --config $(RUN_CONFIG) --formulation $$form --vt-ref $$ref --dp-ref $$dp --out-prefix plots/$$ref || exit $$?; \
			done \
		done \
	done

# Aggregate sweep results across all frames and models into summary CSV and
# bar plots under plots/summary/. Can be re-run without re-running sweeps.
# SUMMARY_SUITES controls which output suites are included (default: param-sweep).
matrix-summary:
	$(PYTHON) matrix_summary.py --runs-dir plots --suites $(SUMMARY_SUITES) --output-dir plots/summary

# ---------------------------------------------------------------------------
# Canonical workflows
# ---------------------------------------------------------------------------

# Sweep (no maps) + summary. Normal daily workflow.
run-matrix-with-summary:
	@$(MAKE) run-matrix
	@$(MAKE) matrix-summary SUMMARY_SUITES=best-fit

# Sweep with full maps + summary. Use for final outputs.
run-matrix-maps-with-summary:
	@$(MAKE) run-matrix-maps
	@$(MAKE) matrix-summary SUMMARY_SUITES=best-fit
