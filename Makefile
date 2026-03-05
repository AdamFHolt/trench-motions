PYTHON ?= python3
VENV_DIR ?= .venv

SMOKE_REF ?= sa
SMOKE_FORMULATION ?= 1
SMOKE_INCLUDE_DP ?= 1
SMOKE_DP_REF ?= 23.5e6
SMOKE_TRANS_STRAIN ?= 1e-13
SMOKE_PSP_FACTOR ?= 0.25
SMOKE_INCLUDE_RP ?= 0

SINGLE_REF ?= sa
SINGLE_FORMULATION ?= 1
SINGLE_INCLUDE_DP ?= 1
SINGLE_DP_REF ?= 23.5e6
SINGLE_TRANS_STRAIN ?= 1e-13
SINGLE_PSP_FACTOR ?= 0.25
SINGLE_ASTHEN_VISC ?= 1e21
SINGLE_LITH_VISC ?= 1e22

QUICK_PREDICTED ?= predictions/new/linear/rms_samodel.DP2.35e+07MPa.l20.0_a10e22.0.viscous_bending.txt
QUICK_OBSERVED ?= data/vt/tnew.sa.dat
QUICK_OUTPUT ?= plots/new/quick/quick_check.png
QUICK_TITLE ?= Manual quick check

MISFIT_CONFIG ?= configs/misfit_smoke.yaml
SINGLE_CONFIG ?= configs/single_smoke.yaml
MATRIX_CONFIGS ?= \
	configs/matrix/linear_hs3.yaml \
	configs/matrix/linear_nnr.yaml \
	configs/matrix/linear_sa.yaml \
	configs/matrix/powerlaw_hs3.yaml \
	configs/matrix/powerlaw_nnr.yaml \
	configs/matrix/powerlaw_sa.yaml
MATRIX_RUNS_DIR ?= runs/matrix
MATRIX_SUMMARY_DIR ?= runs/matrix/summary

.PHONY: venv install smoke single-smoke smoke-config single-smoke-config run-matrix run-matrix-smoke matrix-summary run-matrix-with-summary run-matrix-smoke-with-summary quick-plot test
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
		$(SMOKE_TRANS_STRAIN) \
		$(SMOKE_PSP_FACTOR) \
		$(SMOKE_INCLUDE_RP) \
		--smoke --skip-map

single-smoke:
	$(PYTHON) compute_rates_single.py \
		$(SINGLE_REF) \
		$(SINGLE_FORMULATION) \
		$(SINGLE_INCLUDE_DP) \
		$(SINGLE_DP_REF) \
		$(SINGLE_TRANS_STRAIN) \
		$(SINGLE_PSP_FACTOR) \
		$(SINGLE_ASTHEN_VISC) \
		$(SINGLE_LITH_VISC) \
		--skip-map

smoke-config:
	$(PYTHON) compute_rates_misfit.py --config $(MISFIT_CONFIG)

single-smoke-config:
	$(PYTHON) compute_rates_single.py --config $(SINGLE_CONFIG)

run-matrix:
	@for cfg in $(MATRIX_CONFIGS); do \
		echo "[matrix] $$cfg"; \
		$(PYTHON) compute_rates_misfit.py --config $$cfg || exit $$?; \
	done

run-matrix-smoke:
	@for cfg in $(MATRIX_CONFIGS); do \
		echo "[matrix-smoke] $$cfg"; \
		$(PYTHON) compute_rates_misfit.py --config $$cfg --smoke || exit $$?; \
	done

matrix-summary:
	$(PYTHON) matrix_summary.py --runs-dir $(MATRIX_RUNS_DIR) --output-dir $(MATRIX_SUMMARY_DIR)

run-matrix-with-summary: run-matrix matrix-summary

run-matrix-smoke-with-summary: run-matrix-smoke matrix-summary

quick-plot:
	$(PYTHON) quick_plot.py \
		--predicted $(QUICK_PREDICTED) \
		--observed $(QUICK_OBSERVED) \
		--output $(QUICK_OUTPUT) \
		--title "$(QUICK_TITLE)"

test:
	$(PYTHON) -m unittest discover -s tests -v
