# single-cell-fm-probing

**Do the cell states that single-cell foundation models learn stratify cancer patients better
than a linear baseline — and does the answer survive the controls?**

Single-cell foundation models (scGPT, Geneformer) are trained to reconstruct masked gene
expression across tens of millions of cells. Whether they learn biology beyond what a linear
method on ~2,000 highly-variable genes already captures is disputed: one side reports strong
zero-shot generalisation, the other that HVG + PCA matches or beats them. That argument has been
had on embedding geometry and cell-type clustering. This study has it on a **clinical endpoint**.

## The design in one paragraph

Take a public breast-cancer tumour atlas. Embed every cell three ways — scGPT, Geneformer, and
HVG-PCA as the baseline. Cluster each embedding into cell states and turn each state into a
marker-gene signature. Score every signature in ~1,100 TCGA breast tumours (bulk RNA-seq, a
disjoint patient cohort) and ask whether it predicts survival. Read every score against a
**ladder**: a permutation null, a floor of random gene sets matched for size and expression, the
HVG-PCA baseline, and the PAM50 subtype call as the reference the field already uses. The
question is answered by the margin of the foundation models over the baseline — not by any
number on its own.

Full design, predictions and conventions: [`docs/DESIGN.md`](docs/DESIGN.md).
Literature and novelty: `research/literature.md`.

## Status

Design complete. Pipeline written (`scripts/01`–`06`, orchestrator `run.py`). **Nothing has run.**
Both data sources are public and were verified reachable on 2026-09-28.

## Pipeline

| Stage | Script | Does |
|---|---|---|
| 01 | `01_data.py` | CELLxGENE Census atlas (QC, HVGs, stratified subsample) + TCGA-BRCA bulk, clinical and survival from UCSC Xena; gene harmonisation |
| 02 | `02_embed.py` | per-cell embeddings: `hvg_pca`, `scgpt`, `geneformer`, one schema |
| 03 | `03_states.py` | Leiden clusters per representation → marker-gene signatures, swept over resolution × top-k |
| 04 | `04_translate.py` | signature scores in bulk; age- and stage-adjusted Cox; C-index; family-wise permutation null; matched-random floor; PAM50 reference |
| 05 | `05_ladder.py` | the ladder, patient-bootstrap intervals, predictions P1–P4 |
| 06 | `06_report.py` | figures, candidate shortlist, validation-strategy template |
| — | `run.py` | stages as named tools with a JSON run log (`run.py list`, `run.py all --dry-run`) |

Every stage writes `params.json` beside its outputs: arguments, command, git commit, library
versions, timestamp.

## Setup

```bash
# uv only. Override the machine-wide venv pointer so this project gets its own environment.
export UV_PROJECT_ENVIRONMENT=/scratch/.venv-scfm
export HF_HOME=/scratch/.hf-cache          # model weights land here, not on the root disk
uv sync
uv run python scripts/run.py list
uv run python scripts/run.py all --dry-run  # every stage prints its plan; nothing downloads
```

scGPT is installed per its own instructions behind stage 02 (it carries version-sensitive
dependencies); Geneformer comes from its HuggingFace repository. Both run on CPU with the light
checkpoints named in `scripts/config.py`.

## Data

| | Source | Used |
|---|---|---|
| Single-cell atlas | CELLxGENE Census, `disease == "breast cancer"`, primary data | counts, cell type, donor |
| Bulk RNA-seq | TCGA-BRCA, UCSC Xena `HiSeqV2` | ~1,100 primary tumours |
| Clinical + survival | Xena clinical matrix + TCGA-CDR | age, stage, PAM50, OS, PFI |

The two cohorts share no patients. Data and results are git-ignored.

## Statistical conventions

- The **patient** is the unit of independence: all intervals are patient bootstraps.
- Every score is read **above its null and its floor**; the floor is matched on size and mean
  expression because prognostic signal correlates with both.
- "Best of k" is corrected by a **family-wise** permutation null — the null distribution of the
  *best* signature per representation.
- Age and stage are **adjusted**, not stratified, in the primary model; subtype-stratified models
  are a separate prediction (P3).

## Repository layout

```
docs/DESIGN.md        the design: question, ladder, predictions, stages, conventions
research/             literature and novelty notes
scripts/              01–06 stages, run.py orchestrator, config.py, qc_common.py, surv_common.py
data/  results/       git-ignored
```

## Licence

MIT.
