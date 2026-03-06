PYTHON ?= python3
VENV_DIR ?= .venv

MATRIX_CONFIG ?= configs/matrix.yaml
MATRIX_MAP_CONFIG ?= configs/matrix_maps.yaml
REF_FRAMES ?= hs3 nnr sa
MATRIX_RUNS_DIR ?= results/param-sweep
MATRIX_SUMMARY_DIR ?= results/param-sweep/summary

.PHONY: venv install run-matrix run-matrix-maps matrix-summary run-matrix-with-summary run-matrix-maps-with-summary
venv:
	$(PYTHON) -m venv $(VENV_DIR)

install:
	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install -r requirements.txt

run-matrix:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_CONFIG) --vt-ref $$ref --out-prefix $(MATRIX_RUNS_DIR)/$$ref || exit $$?; \
	done

run-matrix-maps:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix-maps] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_MAP_CONFIG) --vt-ref $$ref --out-prefix results/maps/$$ref || exit $$?; \
	done

matrix-summary:
	$(PYTHON) matrix_summary.py --runs-dir $(MATRIX_RUNS_DIR) --output-dir $(MATRIX_SUMMARY_DIR)

run-matrix-with-summary: run-matrix matrix-summary

run-matrix-maps-with-summary:
	@$(MAKE) run-matrix-maps
	@$(MAKE) matrix-summary MATRIX_RUNS_DIR=results/maps MATRIX_SUMMARY_DIR=results/maps/summary
