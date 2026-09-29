"""Stage 01 — assemble the two cohorts and harmonise their gene vocabulary.

Single-cell: a breast-cancer tumour atlas pulled from the CELLxGENE Census at query time, QC'd,
log-normalised, HVGs flagged, subsampled to a CPU-sized set stratified by cell type.
Bulk: TCGA-BRCA expression, clinical covariates and TCGA-CDR survival from the UCSC Xena hubs,
restricted to primary tumours and joined on patient barcode.

The two cohorts share no patients, so everything downstream is an out-of-cohort test.
"""
import argparse
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import config as C
import qc_common as qc
from surv_common import stage_to_ordinal


def fetch(url: str, dest: Path) -> Path:
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url.split('/')[-1]} …", flush=True)
    urllib.request.urlretrieve(url, dest)
    return dest


# ───────────────────────────────────────────────────────────────── single-cell arm
def build_atlas(out: Path, max_cells: int, seed: int):
    import cellxgene_census
    import scanpy as sc

    print("  querying Census …", flush=True)
    with cellxgene_census.open_soma(census_version=C.CENSUS_VERSION) as census:
        adata = cellxgene_census.get_anndata(
            census, organism=C.CENSUS_ORGANISM,
            obs_value_filter=C.CENSUS_OBS_FILTER,
            column_names={"obs": C.CENSUS_OBS_COLUMNS, "var": ["feature_id", "feature_name"]})
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes before QC", flush=True)
    adata.var_names = adata.var["feature_name"].astype(str).values
    adata.var_names_make_unique()

    # QC: drop low-complexity and high-mitochondrial cells
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    keep = (adata.obs["n_genes_by_counts"] >= C.MIN_GENES_PER_CELL) & \
           (adata.obs["pct_counts_mt"] <= 100 * C.MAX_MT_FRACTION)
    adata = adata[keep].copy()

    # subsample stratified by cell type so rare states survive the CPU budget
    rng = np.random.default_rng(seed)
    if adata.n_obs > max_cells:
        frac = max_cells / adata.n_obs
        idx = (adata.obs.groupby("cell_type", observed=True, group_keys=False)
               .apply(lambda g: g.sample(frac=frac, random_state=int(rng.integers(1 << 31))))
               .index)
        adata = adata[idx].copy()

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=C.N_HVG, flavor="seurat_v3", layer="counts")
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes after QC; "
          f"{int(adata.var['highly_variable'].sum())} HVGs", flush=True)
    adata.write_h5ad(out / "atlas.h5ad")
    return adata


# ────────────────────────────────────────────────────────────────────── bulk arm
def build_bulk(out: Path, raw: Path):
    expr_p = fetch(C.XENA_EXPR, raw / "HiSeqV2.gz")
    clin_p = fetch(C.XENA_CLIN, raw / "BRCA_clinicalMatrix")
    surv_p = fetch(C.XENA_SURV, raw / "TCGA_CDR_survival.tsv")

    with gzip.open(expr_p, "rt") as f:
        expr = pd.read_csv(f, sep="\t", index_col=0)          # genes × samples, log2(RSEM+1)
    primary = [s for s in expr.columns if s.endswith(C.PRIMARY_TUMOUR_SUFFIX)]
    expr = expr[primary]
    expr.columns = [s[:12] for s in expr.columns]              # patient barcode
    expr = expr.loc[:, ~expr.columns.duplicated()]

    clin = pd.read_csv(clin_p, sep="\t", index_col=0)
    clin.index = [s[:12] for s in clin.index]
    clin = clin[~clin.index.duplicated()]
    surv = pd.read_csv(surv_p, sep="\t", index_col=1)         # index = _PATIENT (barcode)
    surv = surv[surv["cancer type abbreviation"] == C.INDICATION]
    surv = surv[~surv.index.duplicated()]

    tbl = pd.DataFrame(index=expr.columns)
    tbl["age"] = pd.to_numeric(clin.reindex(tbl.index)["age_at_initial_pathologic_diagnosis"],
                               errors="coerce")
    tbl["stage_ord"] = stage_to_ordinal(clin.reindex(tbl.index)["pathologic_stage"])
    tbl["subtype"] = clin.reindex(tbl.index).get(C.SUBTYPE_COLUMN)
    for name, (ev, tm) in C.ENDPOINTS.items():
        tbl[ev] = pd.to_numeric(surv.reindex(tbl.index)[ev], errors="coerce")
        tbl[tm] = pd.to_numeric(surv.reindex(tbl.index)[tm], errors="coerce")
    tbl = tbl.dropna(subset=[C.ENDPOINTS[C.PRIMARY_ENDPOINT][1]])
    expr = expr[tbl.index]

    expr.to_parquet(out / "bulk_expr.parquet")
    tbl.to_csv(out / "bulk_clinical.csv")
    ev, tm = C.ENDPOINTS[C.PRIMARY_ENDPOINT]
    print(f"  bulk: {expr.shape[1]:,} primary tumours × {expr.shape[0]:,} genes; "
          f"{int(tbl[ev].sum())} {C.PRIMARY_ENDPOINT} events", flush=True)
    return expr, tbl


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--max-cells", type=int, default=C.MAX_CELLS)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--skip-atlas", action="store_true")
    ap.add_argument("--skip-bulk", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and exit")
    args = ap.parse_args()
    out = Path(args.data_dir); raw = out / "raw"
    out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(json.dumps({"atlas": {"filter": C.CENSUS_OBS_FILTER, "max_cells": args.max_cells},
                          "bulk": {"expr": C.XENA_EXPR, "clin": C.XENA_CLIN, "surv": C.XENA_SURV},
                          "out": str(out)}, indent=2))
        return

    adata = None if args.skip_atlas else build_atlas(out, args.max_cells, args.seed)
    expr, _ = (None, None) if args.skip_bulk else build_bulk(out, raw)

    if adata is not None and expr is not None:
        shared = sorted(set(adata.var_names) & set(expr.index))
        json.dump({"n_shared_genes": len(shared), "genes": shared},
                  open(out / "gene_map.json", "w"))
        print(f"  {len(shared):,} genes shared between atlas and bulk", flush=True)

    qc.record_params(out, args, extra={"indication": C.INDICATION})


if __name__ == "__main__":
    main()
