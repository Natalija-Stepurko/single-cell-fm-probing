"""Stage 01 — assemble the two cohorts and harmonise their gene vocabulary.

Single-cell: a breast-cancer tumour atlas pulled from the CELLxGENE Census at query time, QC'd,
log-normalised, HVGs flagged, subsampled to a CPU-sized set stratified by cell type.
Bulk: TCGA-BRCA expression, clinical covariates and TCGA-CDR survival from the UCSC Xena hubs,
restricted to primary tumours and joined on patient barcode.

The two cohorts share no patients, so everything downstream is an out-of-cohort test.
"""
import argparse
import gzip
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params
from scfm.survival import stage_to_ordinal


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, want: str | None = None) -> Path:
    """Download to <dest>.part and move it into place only once its sha256 matches the pin
    (`want`, else the Xena pin for this file name)."""
    want = want or C.XENA_SHA256.get(dest.name)
    if dest.exists():
        if want is not None and (got := sha256(dest)) != want:
            raise RuntimeError(f"{dest}: sha256 {got} does not match the pinned {want}; "
                               f"delete {dest} and rerun")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    print(f"  downloading {url.split('/')[-1]} …", flush=True)
    urllib.request.urlretrieve(url, part)
    if want is not None and (got := sha256(part)) != want:
        raise RuntimeError(f"{part}: sha256 {got} does not match the pinned {want}; the upstream file "
                           f"changed or the download is corrupt; delete {part} and rerun")
    part.replace(dest)
    return dest


def check_disjoint(adata):
    """The atlas and the bulk cohort must share no patients: no atlas donor may be a TCGA case."""
    tcga = adata.obs["donor_id"].astype(str).str.contains(r"TCGA-")
    assert not tcga.any(), f"{int(tcga.sum())} atlas cells come from TCGA donors"


def atlas_facts(adata) -> dict:
    return {"n_cells": int(adata.n_obs), "n_donors": int(adata.obs["donor_id"].nunique()),
            "n_datasets": int(adata.obs["dataset_id"].nunique()),
            "n_cell_types": int(adata.obs["cell_type"].nunique())}


# ───────────────────────────────────────────────────────────────── single-cell arm
def build_atlas(out: Path, max_cells: int, seed: int):
    """Metadata first, then only the chosen cells' counts: the full query is ~1M cells."""
    import cellxgene_census
    import scanpy as sc

    rng = np.random.default_rng(seed)
    print("  querying Census metadata …", flush=True)
    with cellxgene_census.open_soma(census_version=C.CENSUS_VERSION) as census:
        exp = census["census_data"]["homo_sapiens"]
        obs = exp.obs.read(value_filter=C.CENSUS_OBS_FILTER,
                           column_names=C.CENSUS_OBS_COLUMNS).concat().to_pandas()
        obs = obs[obs["assay"].astype(str).str.startswith(C.ASSAY_PREFIX)
                  & (obs["nnz"] >= C.MIN_GENES_PER_CELL)]
        print(f"  {len(obs):,} candidate cells, {obs['donor_id'].nunique()} donors, "
              f"{obs['dataset_id'].nunique()} datasets", flush=True)
        # oversample by 25% so the mitochondrial filter below can still reach max_cells
        take = min(len(obs), int(max_cells * 1.25))
        pick = (obs.groupby("cell_type", observed=True, group_keys=False)
                .apply(lambda g: g.sample(frac=take / len(obs),
                                          random_state=int(rng.integers(1 << 31)))))
        coords = np.sort(pick["soma_joinid"].to_numpy())
        print(f"  downloading counts for {len(coords):,} cells …", flush=True)
        adata = cellxgene_census.get_anndata(
            census, organism=C.CENSUS_ORGANISM, obs_coords=coords,
            obs_column_names=C.CENSUS_OBS_COLUMNS,
            var_column_names=["feature_id", "feature_name"])
    adata.obs_names = adata.obs_names.astype(str)
    adata.var["ensembl_id"] = adata.var["feature_id"].astype(str).values
    adata.var_names = adata.var["feature_name"].astype(str).values
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    keep = (adata.obs["n_genes_by_counts"] >= C.MIN_GENES_PER_CELL) & \
           (adata.obs["pct_counts_mt"] <= 100 * C.MAX_MT_FRACTION)
    adata = adata[keep].copy()
    if adata.n_obs > max_cells:
        frac = max_cells / adata.n_obs
        idx = (adata.obs.groupby("cell_type", observed=True, group_keys=False)
               .apply(lambda g: g.sample(frac=frac, random_state=int(rng.integers(1 << 31))))
               .index)
        adata = adata[idx].copy()
    adata = finish_atlas(adata)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes after QC; "
          f"{int(adata.var['highly_variable'].sum())} HVGs; "
          f"{adata.obs['donor_id'].nunique()} donors", flush=True)
    check_disjoint(adata)
    adata.write_h5ad(out / "atlas.h5ad")
    # the cell list plus the pinned Census release rebuilds this exact atlas
    C.RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"obs_name": adata.obs_names, "soma_joinid": adata.obs["soma_joinid"].values,
                  "cell_type": adata.obs["cell_type"].astype(str).values,
                  "donor_id": adata.obs["donor_id"].astype(str).values,
                  "dataset_id": adata.obs["dataset_id"].astype(str).values}
                 ).to_csv(C.RESULTS / "atlas_cells.csv", index=False)
    return adata


def finish_atlas(adata):
    """Everything after the final cell set is fixed: gene filter, counts layer, normalisation, HVGs."""
    import scanpy as sc

    # drop genes no retained cell expresses; keeps every downstream matrix smaller
    sc.pp.filter_genes(adata, min_cells=1)

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()   # Geneformer needs it
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=C.N_HVG, flavor="seurat_v3", layer="counts")
    return adata


def atlas_from_cell_list(out: Path, cell_list: Path):
    """Rebuild the atlas from a saved cell list: the listed cells, in the listed order, are the atlas.

    The list already reflects QC and subsampling, so neither is repeated. Counts, HVGs and per-cell
    QC metrics match the original build; the per-gene QC columns in var describe the listed cells
    only, not the oversampled pool they were drawn from (nothing downstream reads them).
    """
    import cellxgene_census
    import scanpy as sc

    cells = pd.read_csv(cell_list, dtype={"obs_name": str})
    coords = cells["soma_joinid"].to_numpy()
    print(f"  downloading counts for {len(coords):,} listed cells …", flush=True)
    with cellxgene_census.open_soma(census_version=C.CENSUS_VERSION) as census:
        adata = cellxgene_census.get_anndata(
            census, organism=C.CENSUS_ORGANISM, obs_coords=coords,
            obs_column_names=C.CENSUS_OBS_COLUMNS,
            var_column_names=["feature_id", "feature_name"])
    pos = pd.Index(adata.obs["soma_joinid"].to_numpy()).get_indexer(coords)
    assert (pos >= 0).all(), f"{int((pos < 0).sum())} listed cells are missing from the Census release"
    adata = adata[pos].copy()
    adata.obs_names = cells["obs_name"].astype(str).to_numpy()
    adata.var["ensembl_id"] = adata.var["feature_id"].astype(str).values
    adata.var_names = adata.var["feature_name"].astype(str).values
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = finish_atlas(adata)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes; "
          f"{int(adata.var['highly_variable'].sum())} HVGs; "
          f"{adata.obs['donor_id'].nunique()} donors", flush=True)
    check_disjoint(adata)
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
    for ev, tm in C.ENDPOINTS.values():
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--max-cells", type=int, default=C.MAX_CELLS)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--cell-list", default=None,
                    help="rebuild the atlas from this cell list (results/atlas_cells.csv); "
                         "no sampling or QC")
    ap.add_argument("--skip-atlas", action="store_true")
    ap.add_argument("--skip-bulk", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and exit")
    args = ap.parse_args(argv)
    out = Path(args.data_dir); raw = out / "raw"
    out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(json.dumps({"atlas": {"filter": C.CENSUS_OBS_FILTER, "max_cells": args.max_cells,
                                    "cell_list": args.cell_list},
                          "bulk": {"expr": C.XENA_EXPR, "clin": C.XENA_CLIN, "surv": C.XENA_SURV},
                          "out": str(out)}, indent=2))
        return

    if args.skip_atlas:
        adata = None
    elif args.cell_list:
        adata = atlas_from_cell_list(out, Path(args.cell_list))
    else:
        adata = build_atlas(out, args.max_cells, args.seed)
    expr, _ = (None, None) if args.skip_bulk else build_bulk(out, raw)

    if adata is not None and expr is not None:
        shared = sorted(set(adata.var_names) & set(expr.index))
        json.dump({"n_shared_genes": len(shared), "genes": shared},
                  open(out / "gene_map.json", "w"))
        print(f"  {len(shared):,} genes shared between atlas and bulk", flush=True)

    extra = {"indication": C.INDICATION}
    if adata is not None:
        extra.update(atlas_facts(adata))
    record_params(out, args, extra=extra)


if __name__ == "__main__":
    main()
