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

Project page: <https://natalija-stepurko.github.io/single-cell-fm-probing/>.
Full design, predictions and conventions: [`docs/DESIGN.md`](docs/DESIGN.md).
Literature and novelty: `research/literature.md`.

## Status

The full run on 50,002 cells and 1,094 patients completed on 2026-09-29. The result is negative:
neither foundation model's cell states predict survival better than the HVG-PCA baseline, and no
representation's best signature clears its family-wise permutation null. A post-hoc sensitivity
analysis that also credits protective signatures gives the same answer. Numbers and figures are
on the project page; the tables behind them are in `results/`.

## Pipeline

| Stage | Script | Does |
|---|---|---|
| 01 | `01_data.py` | CELLxGENE Census atlas (QC, HVGs, stratified subsample) + TCGA-BRCA bulk, clinical and survival from UCSC Xena; gene harmonisation |
| 02 | `02_embed.py` | per-cell embeddings: `hvg_pca`, `scgpt` (via `scgpt_worker.py`), `geneformer`, one schema |
| 03 | `03_states.py` | Leiden clusters per representation → marker-gene signatures, swept over resolution × top-k |
| 04 | `04_translate.py` | signature scores in bulk; age- and stage-adjusted Cox; C-index; family-wise permutation null; matched-random floor; PAM50 reference |
| 05 | `05_ladder.py` | the ladder, patient-bootstrap intervals, predictions P1–P4 |
| 06 | `06_report.py` | figures, candidate shortlist, validation-strategy template |
| — | `run.py` | stages as named tools with a JSON run log (`run.py list`, `run.py all --dry-run`) |

Every stage writes `params.json` beside its outputs: arguments, command, git commit, library
versions, timestamp.

## Setup

Two environments, both built with `uv`. scGPT pins `torchtext` and old `scvi-tools`, so it has its
own Python 3.11 environment; stage 02 calls it as a subprocess.

```bash
# main environment: every stage except the scGPT forward pass
export UV_PROJECT_ENVIRONMENT=/scratch/.venv-scfm
export HF_HOME=/scratch/.hf-cache          # model weights land here, not on the root disk
uv sync

# scGPT environment
uv venv --python 3.11 /scratch/.venv-scgpt
uv pip install --python /scratch/.venv-scgpt/bin/python \
    --index-url https://download.pytorch.org/whl/cpu torch==2.3.1
uv pip install --python /scratch/.venv-scgpt/bin/python torchtext==0.18.0 "numpy<2" \
    "scanpy<1.11" "anndata<0.11" pandas scikit-learn scikit-misc numba "datasets<3" ipython
uv pip install --python /scratch/.venv-scgpt/bin/python --no-deps scgpt==0.2.4
```

Checkpoints download on first use from Hugging Face: Geneformer `Geneformer-V1-10M` with its V1
dictionaries (`ctheodoris/Geneformer`), and the scGPT whole-human checkpoint as released by the
authors' lab (`wanglab/scGPT-human`). Set `SCGPT_PYTHON` if the scGPT environment lives elsewhere.

```bash
uv run python scripts/run.py list
SCFM_SMOKE=1 uv run python scripts/run.py all   # 2,000 cells, short control loops, writes smoke/; ~8 min
uv run python scripts/run.py all                # the study
```

### Runtime

Measured on 4 physical CPU cores (Xeon Platinum 8573C, no GPU):

| Stage | Throughput or time |
|---|---|
| 01 atlas and bulk download | ~1 min for the smoke set |
| 02 Geneformer, bfloat16 | ~66 cells/s |
| 02 scGPT, bfloat16 | ~10 cells/s |
| 03–05 on the smoke set | ~75 s |

Embedding runs in bfloat16 on CPUs with AMX. Against fp32 on 300 cells, per-cell cosine similarity
was at least 0.9999 for both models and 98% of 15-nearest neighbours were unchanged.
`02_embed.py --precision fp32 --limit 300` repeats the check. Cells are sorted by length and packed
into batches, and results are written in shards, so an interrupted stage 02 resumes.

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
docs/index.html       the project page (GitHub Pages); rebuilt by docs/site/build.py
results/              tracked: cell list, signatures, scores, nulls, ladder, report figures, params
                      and run log (5 MB); embeddings and shards stay local
research/             literature and novelty notes
scripts/              01–06 stages, run.py orchestrator, config.py, qc_common.py, surv_common.py
data/  smoke/        git-ignored
```

## Licence

MIT.
