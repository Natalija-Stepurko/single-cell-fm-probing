"""Study configuration: one indication, its data sources, and the sweep grid.

Everything a stage needs to know about *which* study this is lives here, so a second indication
is a second config block, not a second pipeline.
"""
from pathlib import Path

ROOT = Path("/data/scfm")
DATA = ROOT / "data"
RESULTS = ROOT / "results"

INDICATION = "BRCA"

# ---- single-cell atlas: CELLxGENE Census, filtered at query time --------------------------
CENSUS_VERSION = "stable"
CENSUS_ORGANISM = "Homo sapiens"
# Census obs filter syntax; validated on first run (field names follow the CELLxGENE schema).
CENSUS_OBS_FILTER = (
    "disease == 'breast cancer' and is_primary_data == True "
    "and suspension_type == 'cell' and tissue_general == 'breast'"
)
CENSUS_OBS_COLUMNS = ["cell_type", "donor_id", "dataset_id", "assay", "tissue", "disease",
                      "development_stage", "sex"]
MAX_CELLS = 50_000            # stratified by cell_type; CPU budget for the embedding stage
MIN_GENES_PER_CELL = 200
MAX_MT_FRACTION = 0.20
N_HVG = 2_000

# ---- bulk cohort: TCGA-BRCA via UCSC Xena hubs (all three verified live 2026-09-28) --------
XENA_EXPR = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz"
XENA_CLIN = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
XENA_SURV = ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
             "Survival_SupplementalTable_S1_20171025_xena_sp")
PRIMARY_TUMOUR_SUFFIX = "-01"  # TCGA sample-type code for primary solid tumour
ENDPOINTS = {"OS": ("OS", "OS.time"), "PFI": ("PFI", "PFI.time")}
PRIMARY_ENDPOINT = "OS"
COVARIATES = ["age", "stage_ord"]        # adjusted for, not stratified on, in the primary model
SUBTYPE_COLUMN = "PAM50Call_RNAseq"      # the Reference rung; also the P3 strata

# ---- representations -----------------------------------------------------------------
MODELS = {
    "hvg_pca":    {"kind": "baseline", "n_comps": 50},
    "scgpt":      {"kind": "fm", "checkpoint": "scgpt-whole-human", "layers": "last"},
    "geneformer": {"kind": "fm", "checkpoint": "ctheodoris/Geneformer", "variant": "gf-6L-30M-i2048",
                   "layers": "last"},
}

# ---- signature derivation sweep --------------------------------------------------------
LEIDEN_RESOLUTIONS = [0.3, 0.5, 1.0]
MARKER_TOPK = [25, 50, 100]
MIN_CELLS_PER_STATE = 200

# ---- ladder ----------------------------------------------------------------------------
N_PERMUTATIONS = 500          # outcome permutations for the family-wise null
N_FLOOR_SETS = 200            # matched random gene sets per signature
N_BOOTSTRAP = 200             # patient bootstraps for intervals
EXPRESSION_BINS = 10          # for matching random sets on mean expression
SEED = 42
