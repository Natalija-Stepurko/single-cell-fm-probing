"""Study configuration: one indication, its data sources, and the sweep grid.

Everything a stage needs to know about *which* study this is lives here, so a second indication
is a second config block, not a second pipeline.

Paths resolve against ROOT: the repository checkout, or $SCFM_ROOT when set.
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ["SCFM_ROOT"]) if os.environ.get("SCFM_ROOT") else REPO
# SCFM_SMOKE=1 runs every stage end to end on a small atlas and short control loops,
# writing to smoke/ so it never touches the real outputs.
SMOKE = os.environ.get("SCFM_SMOKE") == "1"
DATA = ROOT / ("smoke/data" if SMOKE else "data")
RESULTS = ROOT / ("smoke/results" if SMOKE else "results")

INDICATION = "BRCA"

# ---- single-cell atlas: CELLxGENE Census, filtered at query time --------------------------
CENSUS_VERSION = "2025-11-08"
CENSUS_ORGANISM = "Homo sapiens"
# Every breast-carcinoma label in this release (checked 2026-09-29); TCGA-BRCA spans all of them.
BREAST_CANCER_LABELS = [
    "breast cancer", "invasive ductal breast carcinoma", "triple-negative breast carcinoma",
    "estrogen-receptor positive breast cancer", "invasive lobular breast carcinoma",
    "HER2 positive breast carcinoma",
    "invasive tubular breast carcinoma || invasive lobular breast carcinoma",
    "breast carcinoma", "breast mucinous carcinoma",
]
CENSUS_OBS_FILTER = (
    f"disease in {BREAST_CANCER_LABELS} and is_primary_data == True "
    "and suspension_type == 'cell' and tissue_general == 'breast'"
)
ASSAY_PREFIX = "10x"          # one chemistry family, so counts are on one scale
CENSUS_OBS_COLUMNS = ["soma_joinid", "cell_type", "donor_id", "dataset_id", "assay", "tissue",
                      "disease", "development_stage", "sex", "nnz"]
MAX_CELLS = 2_000 if SMOKE else 50_000   # stratified by cell_type; CPU budget for embedding
MIN_GENES_PER_CELL = 200
MAX_MT_FRACTION = 0.20
N_HVG = 2_000

# ---- bulk cohort: TCGA-BRCA via UCSC Xena hubs (all three verified live 2026-09-28) --------
XENA_EXPR = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz"
XENA_CLIN = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
XENA_SURV = ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
             "Survival_SupplementalTable_S1_20171025_xena_sp")
# sha256 of each download, keyed by the file name it is saved under in data/raw/
XENA_SHA256 = {
    "HiSeqV2.gz": "263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746",
    "BRCA_clinicalMatrix": "39eb3be0fb86e6a577bd2cc01502a7fa5a271e1e1cba294e9dc644ad99580d7f",
    "TCGA_CDR_survival.tsv": "a5e704158bb5c51cded8a368accd999dfb259d428e88bbb0c4386078c5df9617",
}
PRIMARY_TUMOUR_SUFFIX = "-01"  # TCGA sample-type code for primary solid tumour
ENDPOINTS = {"OS": ("OS", "OS.time"), "PFI": ("PFI", "PFI.time")}
PRIMARY_ENDPOINT = "OS"
COVARIATES = ["age", "stage_ord"]        # adjusted for, not stratified on, in the primary model
SUBTYPE_COLUMN = "PAM50Call_RNAseq"      # the Reference rung; also the P3 strata

# ---- representations -----------------------------------------------------------------
MODELS = {
    "hvg_pca":    {"kind": "baseline", "n_comps": 50},
    # whole-human checkpoint as released by the authors' lab; CLS embedding, L2-normalised
    "scgpt":      {"kind": "fm", "hf_repo": "wanglab/scGPT-human",
                   "revision": "a24c237737a40f3720f75abb555489e9fe753be6", "max_len": 1200},
    # 6-layer V1 model (pretrained on ~30M cells); needs the V1 (gc30M) dictionaries
    "geneformer": {"kind": "fm", "hf_repo": "ctheodoris/Geneformer",
                   "revision": "1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5", "variant": "Geneformer-V1-10M",
                   "dict_dir": "geneformer/gene_dictionaries_30m", "max_len": 2048},
}
# scGPT pins torchtext and an older torch, so it runs in its own environment (envs/scgpt).
SCGPT_PYTHON = os.environ.get("SCGPT_PYTHON", str(REPO / "envs/scgpt/.venv/bin/python"))

# ---- embedding runtime (CPU) -------------------------------------------------------------
TORCH_THREADS = 4             # physical cores; hyperthreads do not help GEMMs
EMBED_DTYPE = "bf16"          # AMX-BF16 on this CPU; checked against fp32 in the smoke run
TOKENS_PER_BATCH = 32_768     # cells are sorted by length and packed to this budget
SHARD_CELLS = 250 if SMOKE else 2_000   # results written per shard; a restart resumes


# ---- signature derivation sweep --------------------------------------------------------
LEIDEN_RESOLUTIONS = [0.3, 0.5, 1.0]
MARKER_TOPK = [25, 50, 100]
MIN_CELLS_PER_STATE = 20 if SMOKE else 200

# ---- ladder ----------------------------------------------------------------------------
N_PERMUTATIONS = 20 if SMOKE else 500   # outcome permutations for the family-wise null
N_FLOOR_SETS = 20 if SMOKE else 200     # matched random gene sets per signature
N_BOOTSTRAP = 20 if SMOKE else 200      # patient bootstraps for intervals
N_JOBS = int(os.environ.get("SCFM_N_JOBS", "8"))   # worker processes; results do not depend on it
EXPRESSION_BINS = 10          # for matching random sets on mean expression
SEED = 42
# Analyses added after the pre-registered run each draw from their own generator, so the
# pre-registered streams (SEED: floors, nulls, bootstrap; SEED + 1: sensitivity bootstrap)
# and every number they produce stay unchanged.
SEED_REF_CV = SEED + 4        # cross-validation folds of the reference Cox models
SEED_REF_FLOOR = SEED + 5     # matched-random floor of the proliferation reference
SEED_P2_BOOT = SEED + 6       # selection-aware patient bootstrap (P2)
SEED_P3_PERM = SEED + 7       # within-subtype outcome permutations (P3)
SEED_ADDED_VALUE = SEED + 8   # cross-validation folds of the added-value models
SEED_STRATIFY = SEED + 9      # cross-validation folds of the multivariable stratification
SEED_REF_PERM = SEED + 10     # outcome permutations of the proliferation reference
N_CV_SPLITS = 5
N_CV_REPEATS = 5
N_REF_PERM = 20 if SMOKE else 500
N_P3_PERM = 20 if SMOKE else 1_000
P3_MIN_EVENTS = 10            # subtypes with fewer events are not tested within
ALPHA = 0.05

# PAM50 11-gene proliferation score (Nielsen et al. 2010; ROR-P, Parker et al. 2009), as published;
# aliases are resolved against the bulk symbols at run time
PROLIFERATION_GENES = {
    "BIRC5": [], "CCNB1": [], "CDC20": [], "CDCA1": ["NUF2"], "CEP55": [], "KNTC2": ["NDC80"],
    "MKI67": [], "PTTG1": [], "RRM2": [], "TYMS": [], "UBE2C": [],
}

# ---- prediction thresholds ---------------------------------------------------------------
P3_WITHIN_SUBTYPE_CINDEX = 0.6   # P3 as coded: passes if the C-index exceeds this within any subtype
P4_MAX_HVG_FRAC = 0.5            # P4 as coded: passes if fewer than this fraction of signature genes are HVGs
SINGLE_DONOR_FRAC = 0.8          # flagged single-donor when one donor gives this share of a state's cells
