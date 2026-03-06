PYTHON ?= python3
VENV_DIR ?= env
MAKEFLAGS += --no-print-directory

MATRIX_CONFIG ?= configs/matrix.yaml
MATRIX_MAP_CONFIG ?= configs/matrix_maps.yaml
REF_FRAMES ?= hs3 nnr sa
SUMMARY_SUITES ?= param-sweep

.PHONY: venv install run-matrix run-matrix-maps matrix-summary run-matrix-with-summary run-matrix-maps-with-summary
venv:
	$(PYTHON) -m venv $(VENV_DIR)

install:
	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install -r requirements.txt

run-matrix:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_CONFIG) --vt-ref $$ref --out-prefix plots/$$ref/param-sweep || exit $$?; \
	done

run-matrix-maps:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix-maps] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_MAP_CONFIG) --vt-ref $$ref --out-prefix plots/$$ref/maps || exit $$?; \
	done

matrix-summary:
	$(PYTHON) matrix_summary.py --runs-dir plots --suites $(SUMMARY_SUITES) --output-dir plots/summary

run-matrix-with-summary:
	@$(MAKE) run-matrix
	@$(MAKE) matrix-summary SUMMARY_SUITES=param-sweep

run-matrix-maps-with-summary:
	@$(MAKE) run-matrix
	@$(MAKE) run-matrix-maps
	@$(MAKE) matrix-summary SUMMARY_SUITES=param-sweep,maps
