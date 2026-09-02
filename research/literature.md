# Do single-cell foundation models learn biology beyond highly-expressed genes? — literature and novelty

*Compiled 2026-07-10 via web search (arXiv, bioRxiv, Nature, Semantic Scholar). Citation
counts, venues and IDs are as reported on the retrieval date and should be re-verified
before submission. The pipeline it motivates is described in [`docs/DESIGN.md`](../docs/DESIGN.md).*

---

## 0. The question, and why this substrate

The general question is:

> **Do foundation models — trained on one narrow self-supervised signal — spontaneously
> develop internal representations of properties they never saw, where in the network does
> each property live, and is that organisation shared across independently built models?**

In protein modelling (masked sequence → structure/function) the answer is a **settled,
founding result**, so work there is necessarily methodological. **Single-cell foundation
models are different: the answer is openly disputed as of 2026.** Whether scGPT / Geneformer / UCE learn meaningful biology *beyond what a linear
model on highly-variable genes (HVGs) already captures* is an active fight, with a skeptic
pole, an optimist pole, and a fast-growing interpretability wave. That is exactly why this
substrate was chosen (Topic B): the same protocol here **adjudicates an open question** rather
than re-confirming a closed one.

### What the units map onto

| Unit in sequence/structure domains | Single-cell analogue |
|---|---|
| atom / residue token | **gene token** within a cell |
| protein chain (the pooled unit) | **cell** (CLS or mean-pooled over gene tokens) |
| per-atom / per-residue embedding | per-gene-token embedding |
| pooled structure embedding | cell embedding |
| edge / attention (Rao 2021) | gene–gene attention (co-expression, **not** regulation — see §3) |
| layer axis: encoder → blocks | embedding → N transformer blocks |
| **local** (secondary structure, burial) | **per-gene / within-cell**: HVG status, expression bin, gene-program/GO membership, marker status |
| **global** (fold class, subcellular localisation) | **per-cell**: cell type, tissue, batch (nuisance), cell-cycle phase, disease |
| composition / from-scratch baseline | **HVG-PCA, scVI, logistic-on-HVG** — the contested baselines the FMs must beat |

The locality axis (`residue-local vs whole-protein` in protein models) becomes
**per-gene vs per-cell**, and the "beyond highly-expressed genes" dispute maps onto
**expression/HVG-stratified** probing and geometry.

---

## 1. The dispute — skeptic pole (FMs do not beat simple baselines)

- **Kedzierska, Crawford, Amini, Lu.** *Assessing the limits of zero-shot foundation models
  in single-cell biology.* bioRxiv 2023.10.16.561085 → ***Genome Biology* 2025**
  (doi:10.1186/s13059-025-03574-x). The first systematic **zero-shot** evaluation of
  Geneformer and scGPT. Both are frequently beaten on cell-type clustering (AvgBIO) by
  **HVG selection, scVI, and Harmony**. Two hypotheses offered for the underperformance:
  the masked-LM pretraining does not yield useful *cell* embeddings, or the models failed to
  learn the pretraining task. **The anchor skeptic paper — origin of the "one PCA still rules
  them all" framing.**
- **Boiarsky, Singh, Buendia, Getz, Uhler.** *A Deep Dive into Single-Cell RNA Sequencing
  Foundation Models.* bioRxiv 2023/2024. **Logistic regression on the top-2000 HVGs matches**
  Geneformer/scGPT on cell-type classification across five datasets — motivates our
  HVG / logistic-on-HVG baseline arm as the thing to beat.
- **Souza & Mehta.** *Parameter-free representations outperform single-cell foundation models
  on downstream benchmarks.* arXiv:2602.16696 (2026). **HVG + PCA beats scGPT, Geneformer,
  and UCE** on the **CZ CELLxGENE** downstream benchmarks (cell-type classification,
  perturbation response). The freshest skeptic result; its benchmark fixes the CZ CELLxGENE
  cross-check dataset in our plan so probe numbers are directly comparable to a published
  result, not merely qualitatively aligned.

---

## 2. The dispute — optimist pole (FMs do learn generalizable biology)

- **Rosen, Yanay, … Quake (Tabula Sapiens Consortium).** *Universal cell embedding provides a
  foundation model for cell biology.* ***Nature* 2026** (doi:10.1038/s41586-026-10689-z;
  preprint bioRxiv 2023.11.28.568918; code `github.com/snap-stanford/UCE`). **UCE** — trained
  on ~36M cells across dozens of tissues and 8 species — maps **unseen** species and cell
  types into a coherent universal embedding space zero-shot. Crucially, UCE tokenizes each
  gene by its **ESM-2 protein-sequence embedding**: a built-in cross-modality bridge
  (protein-LM → cell space), making it the one model where a cross-modality convergence test is
  possible on this substrate.
  The strongest "it works" claim in the field, and the paper that prompted this literature
  pass. Checkpoints: 4-layer (CPU-feasible) and 33-layer (GPU-only, the published model).
- **Cui, Wang, … Wang.** *scGPT: toward building a foundation model for single-cell multi-omics
  using generative AI.* ***Nature Methods* 2024**. Value-aware masked model with a
  specialized attention masking scheme for unordered gene sets; whole-human checkpoint
  ≈ 12 layers, d=512, 8 heads, ~50M params. A core model under test.
- **Theodoris, … Ellinor.** *Transfer learning enables predictions in network biology
  (Geneformer).* ***Nature* 2023**. Rank-value encoding (no expression magnitudes),
  BERT-style. V1 6-layer (~10M) / 12-layer (~38M); V2-316M is 18-layer, d=1152, 18 heads.
  A core model under test.
- Context / reviews: *Single-cell foundation models: bringing artificial intelligence into
  cell biology* (PMC12586647, 2026); *Transformer-based Single-Cell Language Model: A Survey*
  (arXiv:2407.13205); *Large-scale foundation model on single-cell transcriptomics*
  (scFoundation, *Nature Methods* 2024) as an additional scale point if a fourth arm is ever
  wanted.

---

## 3. Interpretability wave (2026) — the nearest methodological neighbours to differentiate from

This substrate has a young but fast-moving interpretability literature. None runs *our*
unified probe-vs-geometry-vs-baseline-convergence protocol, but each must be cited and
distinguished so the contribution is not overclaimed.

- **arXiv:2602.22247** — *Multi-Dimensional Spectral Geometry of Biological Knowledge in
  Single-Cell Transformer Representations.* Residual-stream **spectral geometry** across depth
  (e.g. B-cell regulators BATF/BACH2 converging toward the PAX5 identity anchor). Closest
  geometry study; uses functional-map / eigenvalue tools, **not** k-NN-purity/LVR geometry
  paired with supervised decodability against a linear baseline.
- **arXiv:2603.01752** — *Causal Circuit Tracing Reveals Distinct Computational Architectures
  in Single-Cell Foundation Models: … Cross-Model Convergence.* Geneformer (cooperative,
  ~80/20 inhibitory/excitatory, chromatin/RNA-processing hubs) vs scGPT (competitive, ~65/35,
  mitochondrial hubs). Mechanistic circuit tracing + a cross-model convergence claim — but via
  circuits, not CKA/SVCCA layer-grids or probe accuracy.
- **arXiv:2602.17532** — *Systematic Evaluation of Single-Cell Foundation Model Interpretability
  Reveals Attention Captures Co-Expression Rather Than Unique Regulatory Signal.* Attention
  edge-scores carry no regulatory information beyond gene-level features (survives expression
  residualization, degree-preserving nulls, causal ablation). **Load-bearing caveat: treat the
  attention/"edge" arm as diagnostic only, never a regulatory-signal headline.**
- **arXiv:2603.02952** — *Sparse autoencoders reveal organized biological knowledge but minimal
  regulatory logic in single-cell foundation models: a comparative atlas of Geneformer and
  scGPT.* The SAE "phase-2" analogue (cf. Gujral 2025 PNAS in proteins) — a candidate optional
  stage 06.
- **arXiv:2509.14723** — *Transcoder-based Circuit Analysis for Interpretable Single-Cell
  Foundation Models.* Earlier mechanistic-interpretability entry in the same space.

---

## 4. Adjacent / benchmark infrastructure

- **Luecken et al.** *Benchmarking atlas-level data integration in single-cell genomics.*
  ***Nature Methods* 2021** (scIB; doi:10.1038/s41592-021-01336-8). Source of the **Immune**
  (33,506 human immune cells, 10 batches) and **Pancreas** (16,382 cells, 9 technology
  batches) datasets — both with `cell_type` + `batch` labels — and the ARI/NMI/ASW/kBET metric
  family our chance-corrected k-NN-purity / LVR geometry mirrors. h5ad on Figshare
  (`12420968`, `25953868`) and `github.com/theislab/scib`.
- **CZ CELLxGENE Census** — standardized datasets/metrics used by Souza & Mehta (2026); our
  comparability cross-check.
- **Tabula Sapiens** — the UCE consortium's multi-organ human atlas; the natural scale-up set.
- **scVI** (Lopez et al. 2018, *Nature Methods*) and **Harmony** (Korsunsky et al. 2019,
  *Nature Methods*) — the non-FM baselines the skeptics use and we adopt.

---

## 5. Novelty residual — what this project adds

We claim novelty of the **analysis**, not the phenomenon. What is new is that no single work below does **all** of the following *together,
on one matched model set, with one protocol, aimed squarely at the open dispute*:

1. **A depth × linearity law indexed by locality.** Per-gene (local) targets peaking early and
   stored *linearly*; per-cell (global) targets (cell type, tissue, disease) peaking late and
   *non-linearly* (the XGB−linear gap); best extraction layer tracking locality. The
   single-cell instance of a locality-indexed depth-and-linearity law.
2. **Geometry-vs-decodability divergence as a named phenomenon.** The skeptics' "UMAPs don't
   cluster" is a *geometry* claim; we ask whether *decodability* also fails, or whether biology
   is present-but-not-geometrically-organised — quantifying where unsupervised k-NN/LVR/UMAP
   and supervised probe accuracy disagree. The sharpest contribution on this substrate.
3. **Convergence-toward-a-linear-baseline.** Recasting "one PCA still rules them all" as a
   *measurable* CKA/SVCCA convergence of each FM's layers toward HVG-PCA — and, with UCE, a
   cross-modality convergence test. A shared FM map that HVG-PCA does *not* sit inside is the
   cleanest evidence *against* the skeptics; strong FM→PCA convergence is the cleanest evidence
   *for* them.
4. **Beyond-highly-expressed-genes.** Expression/HVG-stratified probes and HVG-ablation of the
   input, testing the specific disputed claim and locating the layer (if any) where non-HVG
   signal appears.

All four use one shared similarity toolkit (`scripts/qc_common.py`: CKA/SVCCA/mutual-kNN/
k-NN-purity/LVR) — applied where the answer is genuinely unknown.

**Positioning for the writeup.** Frame the contribution as a *unified matched-protocol
adjudication* sitting **between** the named poles: cite the optimist pole (UCE, *Nature* 2026)
and the skeptic pole (Kedzierska *Genome Biology* 2025; Boiarsky; Souza & Mehta
arXiv:2602.16696) as the anchors we adjudicate, and the 2026 interpretability wave
(2602.22247 / 2603.01752 / 2602.17532 / 2603.02952) as the nearest methodological neighbours
we extend — do not claim a void.

---

## 6. Honest caveats

- **Not a new phenomenon, and not an empty field.** The 2026 interpretability wave exists;
  the contribution is the unified protocol and its use to adjudicate, not "nobody has looked."
- **Attention-as-edges is weaker here than in proteins.** In proteins contacts genuinely live
  in attention (Rao 2021); in single-cell, attention encodes co-expression, not regulation
  (arXiv:2602.17532). Keep the edge arm diagnostic.
- **Batch effects are a feature, not a bug.** Batch as a *nuisance* target is part of the
  story — good embeddings should *not* organise by batch — report cell-type purity and batch
  purity side by side.
- **CPU + tokenization cost.** scGPT/Geneformer are batch-1 on CPU; UCE needs precomputed
  ESM-2 gene embeddings and is added only after the core pair. Use light checkpoints.
- **Redundancy / leakage.** Use donor/dataset-grouped train/test splits so probe accuracy is
  not memorised composition — the analogue of the protein port's chain-grouped splits.

---

### Provenance
References located via web search on 2026-07-10 (model specs, skeptic/optimist poles,
interpretability wave, benchmark datasets). Verify author lists, venues, DOIs and citation
counts before submission — some venue/year fields are from secondary sources, and the
single-cell interpretability preprints are recent and fast-moving.
