"""Study configuration: one indication, its data sources, and the sweep grid.

Everything a stage needs to know about *which* study this is lives here, so a second indication
is a second config block, not a second pipeline.
"""
import os
from pathlib import Path

ROOT = Path("/data/scfm")
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
PRIMARY_TUMOUR_SUFFIX = "-01"  # TCGA sample-type code for primary solid tumour
ENDPOINTS = {"OS": ("OS", "OS.time"), "PFI": ("PFI", "PFI.time")}
PRIMARY_ENDPOINT = "OS"
COVARIATES = ["age", "stage_ord"]        # adjusted for, not stratified on, in the primary model
SUBTYPE_COLUMN = "PAM50Call_RNAseq"      # the Reference rung; also the P3 strata

# ---- representations -----------------------------------------------------------------
MODELS = {
    "hvg_pca":    {"kind": "baseline", "n_comps": 50},
    # whole-human checkpoint as released by the authors' lab; CLS embedding, L2-normalised
    "scgpt":      {"kind": "fm", "hf_repo": "wanglab/scGPT-human", "max_len": 1200},
    # 6-layer V1 model (pretrained on ~30M cells); needs the V1 (gc30M) dictionaries
    "geneformer": {"kind": "fm", "hf_repo": "ctheodoris/Geneformer", "variant": "Geneformer-V1-10M",
                   "dict_dir": "geneformer/gene_dictionaries_30m", "max_len": 2048,
                   "layers": "last"},
}
# scGPT pins torchtext and old scvi-tools, so it runs in its own environment.
SCGPT_PYTHON = os.environ.get("SCGPT_PYTHON", "/scratch/.venv-scgpt/bin/python")

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
N_JOBS = 8
EXPRESSION_BINS = 10          # for matching random sets on mean expression
SEED = 42
