# Project plan — do single-cell foundation models learn biology beyond highly-expressed genes?

*Focused design for this repository. The full literature review and novelty analysis are the
private working note `../research/literature.md`; the positioning is summarised in §1. This
pipeline applies a layer-resolved geometry-and-probe protocol to a contested substrate.*

## 1. The claim under test

Three single-cell foundation models trained on **self-supervised gene-expression** signals,
with **different objectives and architectures**:

- **scGPT** (value-aware masked modelling; transformer, 12 layers, d=512) — never given cell
  labels.
- **Geneformer** (rank-value masked modelling; BERT-style, 6/12 layers) — no expression
  magnitudes, only rank.
- **UCE** (optional third arm; genes tokenized by **ESM-2 protein embeddings**, 4-layer on
  CPU) — the cross-modality optimist pole.

**The open question (disputed, not settled — this is the point).** Do these models encode
biology *beyond what a linear model on the ~2,000 highly-variable genes (HVGs) already
captures*, and *where in the network* does any such signal live? One pole (**UCE**, Rosen et
al., *Nature* 2026) reports strong zero-shot generalisation; the other (**Kedzierska** *Genome
Biology* 2025; **Boiarsky**; **Souza & Mehta** arXiv:2602.16696) reports that **HVG+PCA
matches or beats** the models. The novelty is the *unified, layer-resolved,
geometry-vs-decodability protocol* that adjudicates *where* and *for which properties* the
foundation models do or do not beat the linear baseline — not the phenomenon itself.

**Falsifiable predictions**
- **P1 — Depth × linearity law indexed by locality.** Per-gene (local) targets peak early and
  decode *linearly*; per-cell (global) targets peak late and *non-linearly* (XGB−linear gap);
  best layer tracks locality.
- **P2 — Geometry-vs-decodability divergence.** There exist regimes where unsupervised
  k-NN/LVR/UMAP organisation and supervised probe accuracy disagree — e.g. cell-type is
  decodable by a probe yet the UMAP "does not cluster."
- **P3 — Convergence toward the linear baseline.** FM layers converge (CKA/SVCCA) toward
  HVG-PCA above permutation baseline. High convergence favours the skeptics; a shared FM map
  that HVG-PCA does *not* sit inside favours the optimists. With UCE, also test whether the
  protein-embedding-tokenized model converges with the expression-tokenized ones.

## 2. Model set (CPU-feasible)

| Arm | Model | Signal | Architecture | Layer axis |
|---|---|---|---|---|
| A | **scGPT** (whole-human) | value-masked | Transformer (~50M) | embed → 12 blocks |
| B | **Geneformer** (V1 6/12-layer; V2-316M only if GPU) | rank-masked | BERT (~10–38M) | embed → 6/12 blocks |
| C (opt.) | **UCE** (4-layer; 33-layer GPU-only) | ESM-2 gene-embedding → cell | Transformer | embed → 4 blocks |
| Baseline | **HVG-PCA / scVI / logistic-on-HVG** | — | linear / VAE | single "layer" |

**Why on CPU.** All three FMs run batch-1 on CPU with the light checkpoints. UCE is added only
after the core scGPT/Geneformer pair works, because it needs precomputed ESM-2 gene-embedding
tables and is the heaviest to stand up. The baseline arm is treated as a pseudo-model with one
"layer" so stage 04 can measure FM→baseline convergence with the same CKA/SVCCA machinery.

## 3. Dataset

- **Primary:** scIB **Immune** (33,506 human immune cells, 10 batches) and **Pancreas**
  (16,382 cells, 9 batches) — Luecken 2021, *Nature Methods*; h5ad on Figshare/`theislab/scib`.
  Ship `cell_type` (global target) and `batch` (nuisance target).
- **Scale-up (optional):** a **Tabula Sapiens** organ subset (multi-tissue) and a
  **CZ CELLxGENE Census** slice — the exact benchmark Souza & Mehta used — for direct
  comparability with the published skeptic result.
- **Redundancy control:** **donor/dataset-grouped** train/test splits so probe accuracy is not
  memorised composition or donor identity.

## 4. Labels (split by locality: per-gene vs per-cell)

| Locality | Property | Source |
|---|---|---|
| **Per-gene (expect early / linear)** | HVG status, mean-expression bin, gene-program / GO membership, TF vs non-TF, marker-gene status | computed from data + GO/Reactome/TF lists |
| **Per-cell (expect late / non-linear)** | Cell type, tissue, disease/condition, cell-cycle phase | dataset annotations |
| **Nuisance (want LOW organisation)** | Batch / donor / technology | dataset metadata |

The "beyond highly-expressed genes" dispute is operationalised as **expression/HVG-stratified**
probing and **HVG-ablation** of the input.

## 5. Pipeline stages

| Stage | Script | Output |
|---|---|---|
| 01 | `01_fetch_cells.py` | scIB h5ad → QC, gene-vocab harmonisation, `cells.jsonl` manifest + per-cell `.npz` + gene-metadata table. |
| 02a | `02_extract_embeddings_scgpt.py` | scGPT all-layer per-gene + cell `.pt` (fp16). |
| 02b | `02_extract_embeddings_geneformer.py` | Geneformer all-layer per-gene + cell `.pt`. |
| 02c | `02c_baselines.py` | HVG-PCA / scVI / logistic-on-HVG embeddings in the same schema. |
| 02d | `02_extract_embeddings_uce.py` (opt.) | UCE 4-layer embeddings. |
| 03 | `03_analyze_embeddings.py` | Per-layer k-NN purity (cell type / batch), per-gene LVR, UMAP, depth law. |
| 04 | `04_convergence.py` | CKA/SVCCA/mutual-kNN, layer×layer: scGPT↔Geneformer(↔UCE) **and each FM↔HVG-PCA**, with permutation baselines. |
| 05 | `05_property_prediction.py` | Linear + XGBoost probes per layer × pooling × property vs the HVG baseline; learning curves; beyond-HVG stratified + ablation probes. |
| 06 | `06_significance.py` | Resampled 95% CIs on stage-04 peaks and FM-vs-baseline probe gaps. |

`scripts/qc_common.py` provides CKA/SVCCA/mutual-kNN/k-NN-purity/LVR; all operate on plain
`[n, d]` arrays, so they apply unchanged to gene-token and cell embeddings. Cross-model
alignment in stage 04 is by **cell barcode** (strictly simpler than the protein port's residue
alignment). Stage 04's FM↔HVG-PCA grid is the genuinely new scientific piece and the direct
test of "one PCA still rules them all."

## 6. Analyses (the paper spine)

1. **Do FM cell embeddings beat HVG-PCA/scVI/logistic-on-HVG?** — best-layer probe table +
   learning curves.
2. **Local→global depth × linearity law** (P1).
3. **Geometry-vs-decodability divergence** (P2) — the sharpest result on this substrate.
4. **Convergence-toward-a-linear-baseline** (P3) — the measurable form of "one PCA rules them
   all"; with UCE, the cross-modality convergence twist.
5. **Beyond-highly-expressed-genes** — expression/HVG-stratified probes + HVG ablation.
6. **Embedding health** — anisotropy / effective rank / collapse as intrinsic reliability
   signals (single-cell embeddings are known anisotropic; diagnostics are load-bearing).

## 7. First milestone (skeleton)

1. `01` on a 5k-cell scIB Pancreas subset (`--limit`); confirm manifest + gene-vocab
   harmonisation + gene-metadata join for both models.
2. `02a`/`02b` on the subset; confirm all-layer per-gene + cell `.pt` extraction and
   resume-skip.
3. `02c` baselines; `03` on both models + baseline — validate the depth split on **cell type
   (global) vs HVG-status (local)**.
4. `04` first scGPT↔Geneformer CKA peak **and** the FM↔HVG-PCA convergence number
   (permutation baseline ≈ 0).
5. `05` probes vs baselines + beyond-HVG; scale to full Immune + Pancreas; add UCE; `06`
   significance.

## 8. Risks & caveats

- **Not a new phenomenon, not an empty field** — the 2026 interpretability wave exists; frame
  the contribution as the unified matched-protocol adjudication *between* named poles.
- **CPU throughput / tokenization** — batch-1 CPU inference; start small; UCE is heaviest
  (needs ESM-2 gene embeddings), add last; V2-316M and UCE-33L are GPU-only.
- **Gene-vocabulary mismatch** — scGPT vs Geneformer vs UCE tokenize differently; per-model
  vocab harmonisation is the main new engineering; cells stay barcode-aligned for stage 04.
- **Attention-as-edges is weak here** — single-cell attention encodes co-expression, not
  regulation (arXiv:2602.17532); keep the edge arm diagnostic.
- **Batch effects are a feature** — report cell-type purity and batch purity side by side.
- **Label/coverage skew** — report per-property coverage explicitly, not a single headline number.
