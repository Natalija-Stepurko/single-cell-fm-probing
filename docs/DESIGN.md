# Design — do foundation-model cell states stratify cancer patients better than a linear baseline?

*Related work is in `../research/literature.md`. This document is the experimental design: what is
measured, against what, and in what order. The discovery cohort is TCGA-BRCA (The Cancer Genome Atlas
breast cohort). Sections 1–10 were committed on 2026-09-28 (commit 8d2c2d0), before any data were
downloaded; this commit is what "pre-specified" refers to (there is no external registry). The later
changes to §1–§10 are one wording edit on 2026-09-29 (5af292a, before any data) and an update of §9
Compute on 2026-09-30 (b86c51d, after the primary run; no prediction, statistic or threshold changed).
Section 11 was added on 2026-09-30, after the primary run. Section 12, added on 2026-10-01, lists where
the first run departed from §1–§10 and what is reported now. Section 13 was committed on 2026-10-01
(2ebb73e), before any signature was scored against outcomes in METABRIC (Molecular Taxonomy of Breast
Cancer International Consortium).*

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

## 12. Deviations from the design and corrections

*Added on 2026-10-01, after the TCGA analysis. A review of the code against §1–§10 found places where
the first run computed something other than what the design describes, or did not compute it. Each entry
gives what the design said, what the first run did, and what is reported now. The as-coded values stay in
`results/` beside the corrected ones (`results/ladder/ladder.json`: `predictions`, `verdicts`, `deviations`;
`results/stratify/cv_summary.json`: `deviations`), and the corrected analyses draw from their own random
streams, so no pre-specified number changed. Numbers are OS in TCGA-BRCA, 1,094 patients and 151 deaths,
unless stated.*

### 12.1 Statistics

**Reference rung (§2).**
- *Design:* the PAM50 subtype call and a published prognostic signature.
- *First run:* one Cox model of PAM50 subtype + age + stage, fitted and scored in sample on the 821
  patients with all three (116 deaths), reported as "PAM50" (C = 0.745) beside unadjusted single-score
  C-indices on 1,094 patients. No published signature.
- *Now:* named references on the same 821 patients, each with its in-sample and cross-validated C
  (5 × 5-fold): age + stage 0.738 in sample, 0.732 cross-validated; PAM50 subtype alone 0.598 and 0.525;
  PAM50 + age + stage 0.745 and 0.723. Age + stage on all 1,072 patients with both: 0.755 and 0.751. The
  published signature is the PAM50 11-gene proliferation score (Nielsen et al. 2010): C = 0.569, above its
  own matched-random floor (mean 0.480), permutation p = 0.0085 (10,000 permutations). The cross-validated
  C is the figure to compare; the lifelines penalty (0.01) and the unpenalised in-sample C are recorded for
  each reference, because they differ for PAM50 alone (0.598 penalised, 0.557 unpenalised).

**Covariate adjustment (§6).**
- *Design:* Cox models adjusted for age and stage in the primary model.
- *First run:* every ladder C-index is Harrell's C of the unadjusted score; the age- and stage-adjusted
  Cox model gave only the hazard ratio and p columns of `results/translate/scores.csv` (1,072 patients).
- *Now:* the ladder stays unadjusted, as run, and says so. Added value over age + stage is reported
  separately (`results/ladder/added_value.csv`): hazard ratio per standard deviation, likelihood-ratio p
  and the change in cross-validated C, once with the signature fixed and once with the selection and
  orientation redone inside every training fold (nested). Every pick was selected on the outcomes of these
  patients, so only the nested change and the published proliferation reference are free of selection.
  Nested, no foundation-model pick or family maximum increases C over age + stage (base C 0.752;
  the changes run from −0.0082 to −0.0001); the proliferation score adds 0.008.

**P1 statistic (§3, §6).**
- *Design:* the family-wise null is the distribution of the best signature's C-index under permuted
  outcomes.
- *First run:* the signature with the largest margin above its floor was compared with the 95th
  percentile of that null. The null describes the family maximum, so this pairs a statistic with a null
  built for a different one.
- *Now:* the family maximum with its family-wise p, and P1 passes if p < 0.05 and the maximising signature
  is above its floor. Pre-specified (risk) direction: p = 0.160 HVG-PCA, 0.126 scGPT, 0.367 Geneformer;
  P1 fails for both foundation models as coded and for all three representations corrected. In the
  post-hoc sensitivity reading (§11), the first run found all three margin-selected picks just below the
  null's 95th percentile; with the matched statistic HVG-PCA (p = 0.024) and scGPT (p = 0.034, borderline) clear it and Geneformer
  (p = 0.052, borderline) does not. The margin-selected pick's p against the same null is also given.

**P2 margin (§2, §3).**
- *Design:* Test against Baseline, both read as their distance above Floor; P2 passes if the bootstrap
  interval of the margin excludes zero.
- *First run:* raw C of the foundation-model pick minus raw C of the baseline pick, both picks held fixed
  in every resample, floors not used: +0.030 for scGPT (interval −0.049 to +0.102), −0.002 for Geneformer.
  These intervals hold the picks fixed, so they leave out the selection over 201–369 signatures.
- *Now:* the above-floor margin with a patient bootstrap (1,000 resamples) that re-selects each pick in
  every resample, floors held at their full-cohort values; P2 passes if the lower 2.5% bound is above zero.
  scGPT −0.021 (−0.078 to +0.012), Geneformer −0.039 (−0.084 to −0.001); in the sensitivity reading
  −0.010 and −0.013, both intervals spanning zero. P2 fails for both, in both readings. The fixed-pick
  intervals are kept as coded and labelled conditional on the selection.

**P3 (§3).**
- *Design:* at least one foundation-model state is prognostic within a PAM50 subtype.
- *First run:* the within-subtype C of the one margin-selected signature per foundation model, passing if
  it exceeds 0.6 in any subtype with at least 10 deaths. The 0.6 and the 10 were in the code before any
  data were downloaded but not in this text. There was no null, and the baseline was not tested. As
  coded, P3 failed for both foundation models in the primary analysis and passed for both in the
  sensitivity analysis.
- *Now:* the maximum within-subtype C over Basal, Her2, LumA and LumB (819 patients), against
  within-subtype outcome permutations, for every representation including the baseline. The run used
  10,000 permutations (`p3_corrected.n_perm`; the `deviations` text in `ladder.json` says 1,000). Primary
  picks: p = 0.38 to 0.40, all fail. Sensitivity picks: HVG-PCA 0.702 in LumA, p = 0.016; scGPT 0.668 in
  LumA, p = 0.050 (borderline); Geneformer 0.627 in LumB, p = 0.149. The pick and its direction were fixed
  on the full cohort, which contains these patients.

**P4 (§3).**
- *Design:* foundation-model signatures that pass P2 are enriched for genes outside the HVG set.
- *First run:* passes if fewer than half of the pick's genes are HVGs, evaluated whether or not P2 passed,
  with no enrichment test. Only 1,526 of the 16,297 genes that markers are drawn from are HVGs, so most
  signatures of every representation, the baseline's included, pass that cut.
- *Now:* not applicable, because P2 fails. Had P2 passed, a hypergeometric test of HVG depletion against
  the scored-gene universe would decide it.

**Monte-Carlo error.**
- *Design:* decisions at α = 0.05.
- *First run:* permutation p-values from 500 or 1,000 permutations, without their Monte-Carlo error.
- *Now:* every permutation p carries its Monte-Carlo standard error and is flagged borderline within two
  standard errors of 0.05. The family-wise null keeps its pre-specified 500 permutations, and a
  10,000-permutation null is reported beside it (§12.4).

**Ladder per setting (§10).**
- *Design:* resolution and marker count are swept, and the ladder is reported at each setting.
- *First run:* one pick per representation across all nine settings, chosen by its margin above the floor,
  which uses the outcomes. The family-wise null pools the same nine settings, so the choice is paid for.
- *Now:* `results/ladder/ladder_by_setting.csv` gives the pick at every setting with that setting's own
  family-wise p. At each of the nine settings, both foundation models' best above-floor margin is below
  the baseline's, in either reading of direction, and in the pre-specified direction no setting's best
  signature clears its own null (smallest p = 0.076, scGPT at resolution 1.0 with 50 genes).

**Clustering stability (§6).**
- *Design:* cell-level clustering is done once; its stability is a separate diagnostic.
- *First run:* the diagnostic was not run. Leiden uses a fixed seed and reproduces the stored cluster
  sizes exactly; states are reported at three resolutions and nine settings.
- *Now:* run post hoc under donor subsampling (§12.4).

**Secondary endpoint (§10).**
- *Design:* the progression-free interval (PFI) is the secondary endpoint, for power; both reported.
- *First run:* only OS was analysed. PFI has fewer events in this cohort (145) than OS (151), so the
  reason given for it does not hold.
- *Now:* PFI runs through translate and ladder with the same statistics (`results/translate_pfi`,
  `results/ladder_pfi`). Pre-specified direction: no representation clears its family-wise null
  (p = 0.178 HVG-PCA, 0.218 scGPT, 0.156 Geneformer).

**Multivariable stratification (not in the design).**
- Added on 2026-10-01: cross-validated Cox models of all of a representation's states with age + stage
  (ridge and gradient-boosted), `results/stratify/`. As first run, the ridge penalty also shrank age and
  stage, and the inner cross-validation chose the smallest penalty in every fold (ridge scGPT 0.7456,
  HVG-PCA 0.7525); the penalty now applies to the state features only. Differences of a few thousandths
  between models are within the spread across fold seeds. The gradient-boosted models are compared only
  with the baseline's gradient-boosted model, because at fixed hyperparameters added features mainly
  dilute age and stage.

### 12.2 Data and states

**Atlas query (§4, §9).**
- *Design:* CELLxGENE Census, `disease == "breast cancer"`, primary data; at most 50,000 cells,
  stratified by cell type.
- *As run:* any of nine breast-carcinoma disease labels (`config.BREAST_CANCER_LABELS`), primary data,
  droplet-based assays from both 3′ and 5′ chemistries, with no batch correction; 50,002 cells from 152 donors in 8
  datasets, 48 cell types, drawn in proportion to cell type (rounding takes the total over 50,000). Four of
  the 8 datasets are one collection (Chen et al. 2026), which supplies 31,552 cells.

**Duplicate cells (not in the design).** 577 atlas cells (1.2%) have a raw count vector identical to
another atlas cell from the same donor in another dataset: Chen et al.'s Global Atlas (`de5416ef`) repeats
cells of its epithelial, immune and stromal compartment datasets (343, 193 and 41 pairs). They were not
removed, so the atlas holds 49,425 unique cells; the cell count, the dataset framing and the state
compositions include the duplicates. `results/report/atlas_duplicates.json` records the count, the
dataset pairs and the duplicated cells in each pick (0 to 78 per pick).

**Gene vocabulary (§4).**
- *Design:* gene vocabulary harmonised to HGNC symbols across atlas, bulk and clinical data.
- *As run:* genes are matched by exact symbol between the Census feature names and TCGA `HiSeqV2`; no alias
  mapping. Marker genes are drawn only from the 16,297 shared genes, so every signature is scored in bulk
  with all its genes. Bulk genes under symbols absent from the atlas are not scored.

**State labels.**
- *First run:* each state was labelled by its most frequent annotated cell type. The baseline's
  pre-specified pick (`hvg_pca` 0.3/c21_k50, "exhausted T cell") comes from one dataset and is 83% one
  donor; its markers (WDR76, TKTL1, CEACAM7, FAM9C, INSL3, …) are not T-cell genes.
- *Now:* `results/states/state_composition.csv` gives each state's donors, top-donor share, datasets,
  chemistry and first and second cell types, and flags a state as single-donor when one donor supplies at
  least 80% of its cells (50 of 123 HVG-PCA states, 6 of 67 scGPT, 10 of 71 Geneformer).

**Cohort separation (§4).** The design states that the cohorts are disjoint by construction. The `data`
stage now asserts that no atlas donor carries a TCGA barcode.

### 12.3 Pipeline and reporting

**Stages (§5).** The scripts `01_data.py` … `06_report.py` and `run.py` became the `scfm` package with the
same logic: stages `data`, `embed`, `states`, `translate`, `ladder`, `report`, plus `stratify` and
`replicate`, run through `scfm run`. The restructure reproduced the tracked results byte for byte from
stage `states` onward.

**Provenance (§5).** The first run's `params.json` recorded a subset of library versions, and some stages
ran from a working tree with uncommitted changes. `params.json` now records the versions of every library
that sets a result, the git commit with a dirty flag, and the pinned revisions of both model
checkpoints; every file download (Xena, METABRIC) is checked against a sha256, and the atlas is pinned to
its Census release. The tracked embeddings were computed (fa3449b, 2026-09-29) before the revisions were
pinned in the configuration, and their `params.json` records no revision; the pinned revisions are the only
snapshots of each model in the download cache that run used.

**Validation template (§8).** Stage 06 writes one three-step template and a shortlist of foundation-model
states above their floor and their family-wise null; the shortlist is empty. Step 1 of the template, the
METABRIC replication, was specified in §13 before it ran and covers the picks, the family maxima and the
proliferation reference.

**Reading against §7.** The primary TCGA result matches none of the three outcome shapes exactly. Test ≈
Baseline, as in the second and third shapes. Every pick lies above its own random-gene floor, which the
third shape (Test ≈ Baseline ≈ Floor) excludes, but no representation's best signature clears its
family-wise null, a rung the outcome table does not name. On the reference the data fit the second shape:
every pick's C-index (0.546–0.591, in sample and selected on these outcomes) is far below age + stage
(cross-validated C 0.751). METABRIC (§13), where every foundation-model margin over the baseline is
negative, points to the second shape.

**METABRIC (§13).** The replication ran as specified. The cBioPortal files were downloaded and checked
for size, coverage and outcome counts before §13 was written; no signature was scored against METABRIC
outcomes before commit `2ebb73e`, and the stage was first tested on permuted outcomes
(`--shuffle-outcomes`). One patient with non-positive OS time is excluded, and one more from disease-specific survival (DSS)
with unknown status. Geneformer's protective pick is scored on 22 of its 25 genes (coverage 0.88, above the 0.8
minimum).

### 12.4 Post hoc additions (after the METABRIC replication)

Added on 2026-10-05 in response to an external review. None of them changes a pre-specified number or
verdict; each draws from its own random stream or reads no outcomes.

**Family-wise null with 10,000 permutations (§6).**
- *Design:* 500 outcome permutations.
- *Added:* `family_wise_10k` in `results/ladder/ladder.json` and `results/ladder_pfi/ladder.json`: one
  set of 10,000 permutations shared by the three representations (seed 57), with the vectorised C-index
  checked against lifelines on 36 permutations (largest difference 0). It is the family-wise p reported on
  the project page and in the README, because its Monte-Carlo standard error near 0.05 is about 0.002
  against about 0.010 with 500. The 500-permutation values stay in `family_wise`, and the two agree on
  every verdict (`verdicts_changed` is empty for OS and PFI). OS, pre-specified direction: p = 0.199
  HVG-PCA, 0.098 scGPT, 0.370 Geneformer. Either direction: 0.025, 0.041 and 0.055; Geneformer's p lies
  2.0 Monte-Carlo standard errors above 0.05, so it does not clear and is reported as borderline. PFI,
  either direction: Geneformer clears (p = 0.049, borderline), as with 500.

**PAM50 subtype association of the picks (outcome-free).** `stages/subtype.py`, called by `ladder`:
for every pick, family maximum and the proliferation reference, the score's median and quartiles per
subtype, Kruskal-Wallis H with ε² = H / (n − 1), and one-vs-rest AUCs; TCGA with the PAM50 call (842
patients), METABRIC with CLAUDIN_SUBTYPE (1,974 patients, NC excluded, no outcomes read). Four of the six
picks track subtype in both cohorts (scGPT's risk pick HER2-enriched, AUC 0.88 in TCGA; Geneformer's risk
pick basal-like, 0.89; the luminal protective picks of HVG-PCA and scGPT low in basal-like tumours, 0.06
and 0.24); the baseline's single-dataset pick and Geneformer's stress state barely differ by subtype
(ε² 0.03 and 0.04). The survival association given PAM50 is the age + stage + PAM50 row of
`added_value.csv`. Outputs: `results/ladder/subtype_association.{csv,json}`.

**Donor mixing of the states (outcome-free composition, margins from translate).** `stages/donors.py`,
called by `ladder`: per representation, the kept states (each once, through its 50-gene signature;
three resolutions pooled, so counts and the Fisher p are descriptive), the number dominated by one donor
(largest donor ≥ 80% of the cells), the median top-donor share, and the best above-floor margin among
multi-donor states (largest donor < 50%) in both readings. Single-donor states: 50 of 123 HVG-PCA, 6 of
67 scGPT, 10 of 71 Geneformer. Best multi-donor margin, risk reading: 0.075, 0.074, 0.035; either
direction: 0.081, 0.091, 0.083. The foundation models produce fewer donor-specific states, and their best
multi-donor states are not above the baseline's in the risk reading and 0.010 and 0.002 above it in the
either-direction reading (single selected values, no interval). Outputs:
`results/ladder/donor_mixing.{csv,json}`.

**Cluster stability under donor subsampling.** New stage `stability` (after `ladder`, not part of
`all`): 20 repeats, each keeping 80% of the 152 donors (122; one draw per repeat shared by the
representations), the 15-NN graph and Leiden clustering rebuilt with the settings of `states`. Adjusted
Rand index against the full clustering restricted to the kept cells: 0.80–0.93 (means per
representation and resolution). Per pick, the best-matching subsample cluster's cell Jaccard averages
0.44 (scGPT's luminal ER pick) to 0.90 (the baseline's single-dataset pick); single-donor picks are
reproduced almost exactly whenever their donor is kept. Outputs: `results/states/stability.csv`,
`stability_picks.csv`, `stability_repeats.csv`, `stability_params.json`.

**Project page.** One page: a summary (question, benchmark, headline result, added value, METABRIC,
the states, novelty, scope, how the analysis was checked) and a technical report below it, with the
deviations of §12 as one table.

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
