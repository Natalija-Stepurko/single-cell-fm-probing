# single-cell-fm-probing

**Do single-cell foundation models learn biology beyond highly-expressed genes — and where in the
network?**

A layer-resolved geometry-and-probe study of single-cell foundation models (**scGPT**,
**Geneformer**, optionally **UCE**) against the linear baselines (**HVG-PCA**, **scVI**,
**logistic-on-HVG**) that sceptics argue still match them.

> **Status: design and literature complete; pipeline stages not yet implemented.**
> What is in this repository today is the shared metric toolkit (`scripts/qc_common.py`), the
> literature review that motivates the design, and the design itself. The numbered stages in
> [Pipeline](#pipeline) are specified but not yet written. This README says so plainly rather than
> describing a pipeline that does not exist.

---

## Why this question, on this substrate

For protein language models, "does the model encode things it was never trained on" is a **settled,
founding result** — so any study there contributes a protocol, not a discovery. Single-cell
foundation models are the opposite: **the answer is genuinely disputed**, and the two poles are
both credible and recent.

- **The sceptic pole.** Kedzierska et al. (*Genome Biology* 2025) find Geneformer and scGPT are
  frequently beaten on zero-shot cell-type clustering by plain HVG selection, scVI and Harmony.
  Boiarsky et al. find **logistic regression on the top-2000 HVGs matches** both models across five
  datasets. Souza & Mehta (arXiv:2602.16696) reach a similar verdict on a CZ CELLxGENE slice.
- **The optimist pole.** UCE (Rosen et al., *Nature* 2026), trained on ~36M cells across dozens of
  tissues and 8 species, maps *unseen* species and cell types into a coherent embedding zero-shot.

So the dispute is not "is there a phenomenon" but **where, and for which properties, the foundation
models actually beat a linear baseline**. That is a question a unified, layer-resolved protocol can
adjudicate — which is the entire reason this repository exists.

The full two-pole survey, with the interpretability literature and the novelty argument, is in
**[`research/literature.md`](research/literature.md)**.

---

## The design, and why each piece is there

The protocol is built so that a null result is as interpretable as a positive one. Each component
answers an objection that would otherwise sink the study:

| Component | The objection it answers |
|---|---|
| **Layer-resolved** extraction, not just the final embedding | "You used the wrong layer." Zero-shot evaluations typically read one layer; if biology lives mid-network, a final-layer test understates the model. |
| **Baselines embedded in the same protocol** (HVG-PCA, scVI, logistic-on-HVG) | "Your baseline was weak." The contested baselines are run through the identical pipeline, not quoted from other papers. |
| **Geometry *and* decodability, reported separately** | The sceptics' "the UMAPs don't cluster" is a claim about *geometry*. Biology may be present but not geometrically organised — so k-NN purity / LVR / UMAP and probe accuracy are measured independently and their disagreement is a result in itself. |
| **Expression- and HVG-stratified probes**, plus HVG ablation | This is the disputed claim stated precisely: does signal exist *beyond* highly-expressed genes, and at which layer does it appear? |
| **Convergence toward the baseline** (CKA/SVCCA/mutual-kNN of each model's layers against HVG-PCA) | Recasts "one PCA still rules them all" as something measurable. Strong convergence to PCA supports the sceptics; a shared model subspace that PCA does *not* sit inside is the cleanest evidence against them. |
| **Donor/dataset-grouped splits** | "Your probe memorised batch or donor." Grouped splits keep accuracy from riding on nuisance structure. |
| **Three similarity metrics, not one** | CKA, SVCCA and mutual k-NN answer different questions (whole representation, linear subspace, local neighbourhoods) and can disagree substantially on the same pair. Reporting one invites the wrong conclusion. |

Full specification — datasets, label set, per-stage detail and caveats — is in
**[`docs/DESIGN.md`](docs/DESIGN.md)**.

---

## The code that exists today

### `scripts/qc_common.py`

The shared representational-similarity toolkit. Every function takes plain `[n_samples, n_features]`
arrays, so the same code serves gene-token embeddings, cell embeddings and baseline features
without modification.

| Function | What it computes | Why it is here |
|---|---|---|
| `column_center(X)` | Feature-column centering | Required before `linear_cka`; separated out so the centering is explicit and auditable rather than hidden. |
| `linear_cka(X, Y)` | Feature-space linear CKA (Kornblith et al. 2019) | Whole-representation similarity. Sensitive to high-variance directions, which is why it is never reported alone. |
| `svcca(Xa, Xb, var=0.99)` | SVD-denoise to `var` energy, then mean CCA correlation (Raghu et al. 2017) | Shared **linear subspace**. Note it needs a permutation null: CCA finds correlated directions between unrelated high-dimensional data, and that null grows with representation width. |
| `mutual_knn(X, Y, k=10)` | Mean fraction of shared neighbours (Huh et al. 2024) | Whether the two representations agree on *local neighbourhoods* — a much stricter criterion than subspace overlap, and the metric of the Platonic Representation Hypothesis. |
| `knn_purity(X, labels, k=15)` | Chance-corrected k-NN label purity | Unsupervised geometry: is a label organised in the embedding without any supervision? Chance-correction matters because label priors are badly skewed in single-cell data. |
| `lvr(X, y, k=15)` | Local variance ratio for a continuous target | The continuous analogue of purity — how much of a continuous target's variance is explained locally. |

---

## Setup

CPU-only; no GPU required for the planned first pass.

```bash
# uv only. Override the global venv so this project gets its own environment.
export UV_PROJECT_ENVIRONMENT=/scratch/.venv-scfm
export UV_CACHE_DIR=/scratch/.uv-cache HF_HOME=/scratch/.hf-cache
export TORCH_HOME=/scratch/.torch-hub TORCHDYNAMO_DISABLE=1
cd /data/scfm
uv sync
```

The `UV_PROJECT_ENVIRONMENT` override matters: if a machine-wide default is set in `~/.bashrc`,
`uv sync` would otherwise install this project's dependencies into a different environment.

---

## Pipeline

**Specified, not yet implemented.** Stages are numbered and resume-safe by design: outputs guarded
by existence checks so a re-run fills gaps rather than recomputing.

| Stage | Script | Purpose |
|---|---|---|
| 01 | `01_fetch_cells.py` | Download scIB `h5ad`, QC, harmonise the gene vocabulary, build a cell manifest + gene metadata. |
| 02a | `02_extract_embeddings_scgpt.py` | scGPT per-layer per-gene and cell embeddings. |
| 02b | `02_extract_embeddings_geneformer.py` | Geneformer per-layer embeddings. |
| 02c | `02c_baselines.py` | HVG-PCA / scVI / logistic-on-HVG baseline features — the contested comparators. |
| 02d | `02_extract_embeddings_uce.py` | *(optional)* UCE 4-layer embeddings (the CPU-feasible checkpoint). |
| 03 | `03_analyze_embeddings.py` | Per-layer geometry: k-NN purity (cell type and batch), LVR, UMAP, depth law. |
| 04 | `04_convergence.py` | CKA/SVCCA/mutual-kNN layer × layer, cross-model **and** model ↔ HVG-PCA. |
| 05 | `05_property_prediction.py` | Linear + XGBoost probes against the baselines; expression/HVG-stratified probes. |
| 06 | `06_significance.py` | Resampled confidence intervals on convergence peaks and model-vs-baseline gaps. |

---

## Data and compute

CPU-only (8 cores, ~165 GB RAM). Starting on scIB **Pancreas** (16.4k cells) and **Immune**
(33.5k cells); optionally scaling to Tabula Sapiens or CZ CELLxGENE. Heavy outputs (`results/`,
`cells/`, `*.pt`, `*.h5ad`) are git-ignored and live on a large disk; the virtual environment and
model caches live on scratch storage.

---

## Repository layout

```
single-cell-fm-probing/
  README.md                 # this file
  docs/DESIGN.md            # full protocol specification and caveats
  research/literature.md    # the two-pole literature survey and novelty argument
  scripts/qc_common.py      # shared similarity/geometry toolkit (implemented)
  pyproject.toml            # dependencies, managed with uv
  results/                  # git-ignored; heavy outputs live off-repo
  cells/                    # git-ignored; downloaded h5ad datasets
```

---

## Licence

MIT © 2026 Natalija Stepurko — see [`LICENSE`](LICENSE).
