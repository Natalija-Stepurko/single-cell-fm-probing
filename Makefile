# Common tasks. Every target runs through uv; set UV_PROJECT_ENVIRONMENT first if the main
# environment should live outside .venv (e.g. export UV_PROJECT_ENVIRONMENT=/scratch/.venv-scfm).
#
#   make setup        main environment, all extras and the dev tools
#   make scgpt-env    the locked scGPT environment at envs/scgpt/.venv
#   make test lint    pytest and ruff
#   make smoke        the whole pipeline on 2,000 cells with short control loops (~10 min; network)
#   make all          the study, stages data -> report (hours on CPU; network)
#   make reproduce    rebuild everything downstream of the tracked signatures, no models (~1.5 h)
#   make verify       compare results/ with results/MANIFEST.sha256
#   make page         rebuild docs/index.html from results/

UV  ?= uv
RUN := $(UV) run --locked
# the atlas query and the embeddings need the data and embed extras
RUN_FULL := $(UV) run --locked --all-extras

.PHONY: setup scgpt-env test lint smoke all reproduce verify page

setup:
	$(UV) sync --locked --all-extras --group dev

scgpt-env:
	cd envs/scgpt && env -u UV_PROJECT_ENVIRONMENT $(UV) sync --locked

test:
	$(RUN) pytest -q -rs

lint:
	$(RUN) ruff check src tests

smoke:
	SCFM_SMOKE=1 $(RUN_FULL) scfm run all

all:
	$(RUN_FULL) scfm run all

# bulk cohort from the pinned Xena downloads, then translate -> ladder -> stratify -> report (and the PFI
# translate -> ladder) from the tracked results/states; the atlas and the embeddings are not needed.
# Regenerated outputs are removed first, so a file that is no longer written shows as missing in verify.
REGENERATED := results/translate results/ladder results/stratify results/report results/translate_pfi \
               results/ladder_pfi
reproduce:
	rm -rf $(REGENERATED)
	$(RUN) scfm run data -- --skip-atlas
	$(RUN) scfm run translate
	$(RUN) scfm run ladder
	$(RUN) scfm run stratify
	$(RUN) scfm run report
	$(RUN) scfm run translate -- --endpoint PFI
	$(RUN) scfm run ladder -- --endpoint PFI

verify:
	$(RUN) python -m scfm.verify

page:
	$(RUN) python docs/site/build.py
