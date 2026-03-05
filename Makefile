PYTHON ?= python3

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

.PHONY: smoke single-smoke quick-plot test
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

quick-plot:
	$(PYTHON) quick_plot.py \
		--predicted $(QUICK_PREDICTED) \
		--observed $(QUICK_OBSERVED) \
		--output $(QUICK_OUTPUT) \
		--title "$(QUICK_TITLE)"

test:
	$(PYTHON) -m unittest discover -s tests -v
