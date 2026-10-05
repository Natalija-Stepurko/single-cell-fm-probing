import pandas as pd

from scfm import config as C
from scfm.stages.states import composition, marker_shares
from scfm.stages.translate import endpoint_dir, resolve_proliferation


def test_endpoint_dir_keeps_primary_outputs_apart():
    assert endpoint_dir("translate", C.PRIMARY_ENDPOINT) == "translate"
    assert endpoint_dir("ladder", "PFI") == "ladder_pfi"


def test_proliferation_aliases_resolve_to_bulk_symbols():
    used = resolve_proliferation(["BIRC5", "NUF2", "KNTC2", "NDC80", "MKI67"])
    assert used["CDCA1"] == "NUF2"
    assert used["KNTC2"] == "KNTC2"          # the published symbol wins when present
    assert used["BIRC5"] == "BIRC5" and "UBE2C" not in used


def test_marker_shares():
    s = marker_shares(["RPL13", "RPS6", "MT-CO1", "MALAT1", "CD3E", "RPLP0"])
    assert s["ribo"] == 2 / 6                # RPLP0 is not ^RP[SL]\d
    assert s["mito"] == 1 / 6 and s["lnc"] == 1 / 6


def test_composition_flags_single_donor_states():
    obs = pd.DataFrame({
        "leiden": ["0"] * 10 + ["1"] * 4,
        "donor_id": ["d1"] * 9 + ["d2"] + ["d1", "d2", "d3", "d4"],
        "dataset_id": ["s1"] * 14, "assay": ["10x 3' v3"] * 14,
        "cell_type": ["malignant cell"] * 8 + ["T cell"] * 2 + ["B cell"] * 4,
        "n_genes_by_counts": range(14),
    })
    rows = {r["cluster"]: r for r in composition(obs, "leiden", min_cells=5)}
    assert rows["0"]["single_donor"] and rows["0"]["top_donor_frac"] == 0.9
    assert rows["0"]["second_cell_type"] == "T cell" and rows["0"]["kept"]
    assert not rows["1"]["single_donor"] and not rows["1"]["kept"]
    assert rows["1"]["donor_entropy"] == 2.0
