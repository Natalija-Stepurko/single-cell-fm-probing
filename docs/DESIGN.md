# Design — do foundation-model cell states stratify cancer patients better than a linear baseline?

*The literature and novelty analysis live in `../research/literature.md`. This document is the
experimental design: what is measured, against what, and in what order. Sections 1–10 were
written before any data were seen. Section 11 was added on 2026-09-30, after the primary run.*

## 1. The question

Single-cell foundation models (scGPT, Geneformer) are trained to reconstruct masked gene
expression across tens of millions of cells. Whether the cell states they learn carry biology
beyond what a linear method on ~2,000 highly-variable genes (HVGs) already captures is disputed:
one pole reports strong zero-shot generalisation, the other that HVG + PCA matches or beats them
(`literature.md` §1–2).

That dispute has been argued on embedding geometry and cell-type clustering. This study asks it
on a **clinical endpoint**: take the cell states each representation finds in a tumour
single-cell atlas, translate them into gene signatures, score those signatures in bulk tumour
RNA-seq from an independent patient cohort, and ask whether they **stratify patient survival** —
and whether the foundation-model signatures add prognostic value over the linear baseline's.

**Indication: breast cancer.** Largest TCGA cohort with RNA-seq and survival (~1,100 patients);
a well-annotated public single-cell atlas; and a canonical published reference (PAM50 intrinsic
subtypes) that gives the ladder a fixed rung the field already accepts.

## 2. Why this is not just another probe

Every number in this study is read against a ladder, because a hazard ratio or a concordance
index means nothing on its own:

| Rung | What it is | What it establishes |
|---|---|---|
| **Null** | survival outcomes permuted across patients | what any signature scores by chance |
| **Floor** | 200 random gene sets matched to each signature for size and mean expression | what *any* gene set of that shape scores — the honest zero |
| **Baseline** | cell-state signatures derived from **HVG-PCA** clusters | what the linear method the sceptics defend achieves |
| **Test** | cell-state signatures derived from **scGPT** and **Geneformer** clusters | the claim under test |
| **Reference** | PAM50 subtype call and a published prognostic signature | where the field already is |

The comparison that decides the question is **Test against Baseline**, both read as their
distance above Floor. Test against Reference says whether either is clinically interesting.

**The shared-input confound, carried over from the protein study.** Every representation here
— scGPT, Geneformer, HVG-PCA — is built from the *same* expression matrix, and the bulk scores
are computed on the *same* genes. Agreement between them, or prognostic signal in all of them,
can therefore come from the shared input and not from anything a model learned. The
Baseline rung is the control for this: it is what the shared input yields with no learned model
at all. Any FM claim is stated as the margin over that rung, never as an absolute.

## 3. Falsifiable predictions

- **P1.** FM-derived signatures reach a concordance index above the Floor. *If not, the FM cell
  states carry no prognostic biology at all.*
- **P2.** FM-derived signatures exceed HVG-PCA-derived signatures by a margin whose bootstrap
  interval excludes zero. *If not, "one PCA rules them all" holds on this endpoint.*
- **P3.** At least one FM-derived state is prognostic **within** a PAM50 subtype. *If so, it is
  biology the reference does not already encode; if not, the FM has rediscovered PAM50.*
- **P4.** FM signatures that pass P2 are enriched for genes *outside* the HVG set. *If so, the
  "beyond highly-expressed genes" claim has a concrete instance.*

## 4. Data

| Arm | Source | What is used |
|---|---|---|
| Single-cell tumour atlas | CELLxGENE Census, `disease == "breast cancer"`, primary data only | raw counts, `cell_type`, `donor_id`, `dataset_id` |
| Bulk tumour RNA-seq | TCGA-BRCA via UCSC Xena, `HiSeqV2` (log2 RSEM+1) | ~1,100 tumour samples |
| Clinical | Xena `BRCA_clinicalMatrix` + TCGA-CDR survival table | age, stage, PAM50, OS / OS.time, PFI / PFI.time |

Gene vocabulary is harmonised to HGNC symbols across all three. Single-cell and bulk cohorts
are **disjoint patients** by construction, so every bulk result is an out-of-cohort test.

## 5. Pipeline

| Stage | Script | Output |
|---|---|---|
| 01 | `01_data.py` | atlas `.h5ad` (QC'd, HVGs flagged), bulk expression matrix, clinical table, gene map |
| 02 | `02_embed.py` | per-cell embeddings for `scgpt`, `geneformer`, `hvg_pca`, one schema |
| 03 | `03_states.py` | Leiden clusters per representation → marker-gene signatures (`signatures.json`) |
| 04 | `04_translate.py` | signature scores in bulk; Cox models (age + stage adjusted); C-index; permutation null; matched-random floor |
| 05 | `05_ladder.py` | the ladder assembled, bootstrap intervals over patients, P1–P4 tested |
| 06 | `06_report.py` | figures, candidate shortlist, validation-strategy template |
| — | `run.py` | orchestrator: runs stages as tools, writes a JSON run log |

Every stage writes `params.json` beside its outputs (arguments, command, git commit, library
versions, timestamp) via `qc_common.record_params`.

## 6. Statistical conventions

Carried over from the protein study, where each was learned the hard way:

- **The patient is the unit of independence.** Intervals are bootstraps over patients, never over
  genes or cells. Cell-level clustering is done once; its stability is a separate diagnostic.
- **Every score is read above its null and its floor.** The floor is matched on size and mean
  expression because prognostic signal correlates with both.
- **Multiple testing by family-wise permutation.** For each representation the null is the
  distribution of the *best* signature's C-index under permuted outcomes, so "best of k" is
  controlled at the ladder level, not per signature.
- **Adjustment, not stratification, for age and stage** in the primary model; subtype-stratified
  models are P3 and reported separately.
- **Nothing is compared across cohorts of different size** without saying so; C-index is
  budget-sensitive.

## 7. What a positive and a negative result each look like

- **FM wins.** Test > Baseline above Floor, interval excludes zero, at least one state prognostic
  within a PAM50 subtype, signature enriched outside HVGs. Shortlist those states' marker genes.
- **Sceptics win.** Test ≈ Baseline, both ≪ Reference. Report it as the finding: on a clinical
  endpoint, the linear method is sufficient. This is a publishable negative.
- **Both fail.** Test ≈ Baseline ≈ Floor. The cell-state → bulk translation loses the signal;
  the study cannot separate the poles on this design.

## 8. Validation strategy (what stage 06 templates)

For each shortlisted state: (1) replicate the Cox result in METABRIC, an independent bulk
cohort with survival; (2) check the marker genes against the state's cell-type annotation in the
atlas; (3) name the orthogonal assay that would confirm the state in tissue.

## 9. Compute

CPU-only. Light checkpoints (scGPT whole-human, Geneformer 6-layer V1). The atlas is subsampled
to ≤50,000 cells stratified by cell type. Embedding runs in bfloat16 on length-sorted batches;
setup and measured throughput are in the README. The bulk arm is cheap.

## 10. Risks

- **Signature-derivation choices dominate.** Cluster resolution and marker-gene count set the
  signatures; both are swept, and the ladder is reported at each setting.
- **Bulk deconvolution is approximate.** A signature score in bulk is a proxy for state abundance;
  stated as such throughout.
- **TCGA-BRCA survival is right-censored with ~15% events.** PFI is the secondary endpoint for
  power; both reported.
- **Atlas composition ≠ TCGA composition.** Subtype mix differs; P3 addresses this directly.

## 11. Sensitivity analysis: direction of effect (added after the primary run)

The pre-registered ladder reads the concordance index in one direction: a signature earns a high
C-index only when a higher score means worse survival. The design did not say whether protective
signatures, where a higher score means longer survival, should count. The primary run showed that
they are common and include the largest departures from 0.5 (C down to 0.40).

The sensitivity analysis reads every signature in the direction it acts on the full cohort:

- a protective signature is scored as 1 − C, and each of its matched random sets is scored the
  same way, so the floor asks whether the signature beats random genes in its own direction;
- the family-wise null is the best max(C, 1 − C) across a representation's signatures, computed on
  the same 500 outcome permutations as the primary null;
- bootstrap intervals and the within-subtype test (P3) keep each signature's direction fixed.

The primary ladder and predictions are unchanged and remain the pre-registered result. The
sensitivity ladder is reported beside them and labelled as post hoc.

## 13. Independent replication in METABRIC (specified before the run)

*Written and committed on 2026-10-01, after the TCGA analysis and before any frozen signature was scored
against METABRIC outcomes. The stage refuses to run on real outcomes unless the frozen file below is committed
and unmodified.*

**Cohort.** METABRIC (Curtis et al. 2012; Pereira et al. 2016) from the cBioPortal datahub, `brca_metabric`, pinned to
datahub commit `dca75cb3f32b82d54a6f78bf0a6323e5b975aca1`, files checked by sha256. Illumina HT-12 log2 intensities,
1,980 patients with expression. Duplicate gene symbols collapse to the row with the highest mean; missing values take
the gene's mean; every gene is z-scored across patients. METABRIC is ODbL-licensed and downloaded at run time; only
aggregate statistics are written to `results/`.

**Endpoints.** Overall survival (OS) is primary, as in TCGA. Disease-specific survival (DSS) is secondary: deaths
from other causes are censored at death, because METABRIC's long follow-up makes many deaths unrelated to the cancer.

**What is frozen.** `results/replicate/frozen_signatures.json` lists every signature with its genes and the direction
it was read in TCGA. Nothing is re-selected or re-oriented in METABRIC. Signatures with identical gene sets are one
test.

| Signature | Direction | Genes | TCGA C (in its direction) | Role |
|---|---|---|---|---|
| `hvg_pca:0.3/c21_k50` | risk | 50 | 0.547 | primary pick |
| `scgpt:1.0/c11_k50` | risk | 50 | 0.576 | primary pick, family maximum (primary) |
| `geneformer:0.3/c14_k100` | risk | 100 | 0.546 | primary pick |
| `hvg_pca:1.0/c17_k100:protective` | protective | 100 | 0.591 | sensitivity pick, P3 lead (within LumA) |
| `scgpt:1.0/c8_k50:protective` | protective | 50 | 0.583 | sensitivity pick, P3 lead (within LumA) |
| `geneformer:1.0/c14_k25:protective` | protective | 25 | 0.578 | sensitivity pick, P3 lead (within LumB) |
| `hvg_pca:1.0/c50_k50` | risk | 50 | 0.571 | family maximum (primary) |
| `geneformer:0.5/c13_k25` | risk | 25 | 0.557 | family maximum (primary) |
| `hvg_pca:0.3/c14_k50:protective` = 0.5/c21_k50, 1.0/c23_k50 | protective | 50 | 0.600 | family maximum (sensitivity) |
| `scgpt:1.0/c5_k50:protective` | protective | 50 | 0.593 | family maximum (sensitivity) |
| `geneformer:1.0/c14_k100:protective` | protective | 100 | 0.590 | family maximum (sensitivity) |
| `reference:proliferation` | risk | 11 | 0.569 | reference |

**Statistics, per signature and endpoint.** A signature is scored only if at least 80% of its genes are measured.
The concordance index is read in the TCGA direction. A floor of 200 random gene sets matched on size and
expression bin is drawn in METABRIC and read the same way. An age-adjusted Cox model, stratified by METABRIC
cohort, gives the hazard ratio per standard deviation of the score in the TCGA direction.

**Replication criterion.** A signature replicates on an endpoint if its oriented C exceeds its METABRIC floor's
95th percentile *and* its age-adjusted, cohort-stratified hazard ratio is in the TCGA direction with p < 0.05.

**P3 leads.** Every sensitivity-analysis pick is tested within the PAM50 subtype where its TCGA within-subtype C was
highest (HVG-PCA and scGPT in luminal A, Geneformer in luminal B), whatever its TCGA P3 result: oriented C and a
within-subtype outcome-permutation p (1,000 permutations). METABRIC's subtype call is PAM50 plus claudin-low; the
named subtypes are taken as they are.

**Foundation model against baseline.** For the primary and the sensitivity picks, the difference in floor-adjusted
C between each foundation model's pick and HVG-PCA's pick, with a paired patient bootstrap (1,000 resamples; floors
held at their METABRIC values). A foundation-model advantage requires the interval to lie above zero.

**How the outcomes will be read.**
- A replicated signature is evidence that its programme is prognostic beyond one cohort and one platform; it is not
  evidence for a foundation-model advantage unless the margin interval also lies above zero.
- If the protective picks of the baseline and of scGPT both replicate, that is one luminal programme found twice,
  reported as such.
- A pick that fails the criterion is reported as not replicated, with its numbers; no other signature is tested in
  its place.

Seeds: floor 53, P3 permutations 54, margin bootstrap 55.
