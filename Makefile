PYTHON ?= python3
VENV_DIR ?= .venv

SMOKE_REF ?= sa
SMOKE_FORMULATION ?= 1
SMOKE_INCLUDE_DP ?= 1
SMOKE_DP_REF ?= 23.5e6
SMOKE_INCLUDE_RP ?= 0

SINGLE_REF ?= sa
SINGLE_FORMULATION ?= 1
SINGLE_INCLUDE_DP ?= 1
SINGLE_DP_REF ?= 23.5e6
SINGLE_ASTHEN_VISC ?= 1e21
SINGLE_LITH_VISC ?= 1e22

QUICK_PREDICTED ?= archive/generated/predictions/new/linear/rms_samodel.DP2.35e+07MPa.l20.0_a10e22.0.viscous_bending.txt
QUICK_OBSERVED ?= data/vt/tnew.sa.dat
QUICK_OUTPUT ?= results/quick/quick_check.png
QUICK_TITLE ?= Manual quick check

MISFIT_CONFIG ?= configs/misfit_smoke.yaml
SINGLE_CONFIG ?= configs/single_smoke.yaml
MATRIX_CONFIG ?= configs/matrix.yaml
MATRIX_MAP_CONFIG ?= configs/matrix_maps.yaml
REF_FRAMES ?= hs3 nnr sa
MATRIX_RUNS_DIR ?= results/param-sweep
MATRIX_SUMMARY_DIR ?= results/param-sweep/summary

.PHONY: venv install smoke single-smoke smoke-config single-smoke-config run-matrix run-matrix-smoke run-matrix-maps matrix-summary run-matrix-with-summary run-matrix-smoke-with-summary run-matrix-maps-with-summary quick-plot
venv:
	$(PYTHON) -m venv $(VENV_DIR)

install:
	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install -r requirements.txt

smoke:
	$(PYTHON) compute_rates_misfit.py \
		$(SMOKE_REF) \
		$(SMOKE_FORMULATION) \
		$(SMOKE_INCLUDE_DP) \
		$(SMOKE_DP_REF) \
		$(SMOKE_INCLUDE_RP) \
		--smoke --skip-map

single-smoke:
	$(PYTHON) compute_rates_single.py \
		$(SINGLE_REF) \
		$(SINGLE_FORMULATION) \
		$(SINGLE_INCLUDE_DP) \
		$(SINGLE_DP_REF) \
		$(SINGLE_ASTHEN_VISC) \
		$(SINGLE_LITH_VISC) \
		--skip-map

smoke-config:
	$(PYTHON) compute_rates_misfit.py --config $(MISFIT_CONFIG)

single-smoke-config:
	$(PYTHON) compute_rates_single.py --config $(SINGLE_CONFIG)

run-matrix:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_CONFIG) --vt-ref $$ref --out-prefix $(MATRIX_RUNS_DIR)/$$ref || exit $$?; \
	done

run-matrix-smoke:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix-smoke] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_CONFIG) --vt-ref $$ref --out-prefix $(MATRIX_RUNS_DIR)/$$ref --smoke || exit $$?; \
	done

run-matrix-maps:
	@for ref in $(REF_FRAMES); do \
		echo "[matrix-maps] $$ref"; \
		$(PYTHON) compute_rates_misfit.py --config $(MATRIX_MAP_CONFIG) --vt-ref $$ref --out-prefix results/maps/$$ref || exit $$?; \
	done

matrix-summary:
	$(PYTHON) matrix_summary.py --runs-dir $(MATRIX_RUNS_DIR) --output-dir $(MATRIX_SUMMARY_DIR)

run-matrix-with-summary: run-matrix matrix-summary

run-matrix-smoke-with-summary: run-matrix-smoke matrix-summary

run-matrix-maps-with-summary:
	@$(MAKE) run-matrix-maps
	@$(MAKE) matrix-summary MATRIX_RUNS_DIR=results/maps MATRIX_SUMMARY_DIR=results/maps/summary

quick-plot:
	$(PYTHON) quick_plot.py \
		--predicted $(QUICK_PREDICTED) \
		--observed $(QUICK_OBSERVED) \
		--output $(QUICK_OUTPUT) \
		--title "$(QUICK_TITLE)"
