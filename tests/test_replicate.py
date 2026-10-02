import numpy as np
import pandas as pd
import pytest

from scfm import config as C
from scfm.stages.replicate import analyse, load_metabric, shuffle_outcomes, usable

PATIENT = """#Patient Identifier\tAge\tCohort\tOS months\tOS status\tSubtype\tVital
#desc\tdesc\tdesc\tdesc\tdesc\tdesc\tdesc
#STRING\tNUMBER\tSTRING\tNUMBER\tSTRING\tSTRING\tSTRING
#1\t1\t1\t1\t1\t1\t1
PATIENT_ID\tAGE_AT_DIAGNOSIS\tCOHORT\tOS_MONTHS\tOS_STATUS\tCLAUDIN_SUBTYPE\tVITAL_STATUS
MB-1\t50.5\t1\t10\t1:DECEASED\tLumA\tDied of Disease
MB-2\t61\t1\t20\t0:LIVING\tBasal\tLiving
MB-3\t70\t2\t30\t1:DECEASED\tLumB\tDied of Other Causes
MB-4\t45\t2\t\t\tHer2\t
MB-5\t55\t3\t5\t1:DECEASED\tclaudin-low\tDied of Disease
"""
SAMPLE = """#Patient Identifier\tSample Identifier
#d\td
#STRING\tSTRING
#1\t1
PATIENT_ID\tSAMPLE_ID\tTUMOR_STAGE
MB-1\tMB-1\t2
MB-2\tMB-2\t1
MB-3\tMB-3\t
MB-4\tMB-4\t3
MB-5\tMB-5\t2
"""
EXPR = """Hugo_Symbol\tEntrez_Gene_Id\tMB-1\tMB-2\tMB-3\tMB-4\tMB-6
G1\t1\t1.0\t2.0\t3.0\t4.0\t9.0
G1\t11\t5.0\t6.0\t7.0\tNA\t9.0
G2\t2\t2.0\tNA\t4.0\t6.0\t1.0
\t3\t1.0\t1.0\t1.0\t1.0\t1.0
G3\t4\t1.0\t1.0\t1.0\t1.0\t1.0
"""


@pytest.fixture
def files(tmp_path):
    out = {}
    for name, text in [("data_clinical_patient.txt", PATIENT), ("data_clinical_sample.txt", SAMPLE),
                       ("data_mrna_illumina_microarray.txt", EXPR)]:
        (tmp_path / name).write_text(text)
        out[name] = tmp_path / name
    return out


def test_loader_on_a_tiny_synthetic_file(files):
    z, gene_mean, clin = load_metabric(files)
    assert list(z.columns) == ["MB-1", "MB-2", "MB-3", "MB-4"]       # MB-6 has no clinical row
    assert sorted(z.index) == ["G1", "G2", "G3"]                     # unnamed row dropped
    # the duplicate G1 with the higher mean (second row) is kept; its NaN gets that row's mean
    assert gene_mean["G1"] == pytest.approx(np.mean([5, 6, 7, 6]))
    assert gene_mean["G2"] == pytest.approx(4.0)                    # NaN filled with the row mean
    np.testing.assert_allclose(z.loc["G1"].mean(), 0, atol=1e-12)
    assert (z.loc["G3"] == 0).all()                                 # constant gene
    assert clin.loc["MB-1", "time"] == pytest.approx(10 * 30.4375)
    assert clin["os_event"].tolist()[:3] == [1.0, 0.0, 1.0]
    assert clin["dss_event"].tolist()[:3] == [1.0, 0.0, 0.0]        # other-cause death censored
    assert np.isnan(clin.loc["MB-4", "os_event"]) and np.isnan(clin.loc["MB-4", "dss_event"])
    assert clin.loc["MB-1", "age"] == 50.5 and clin.loc["MB-3", "cohort"] == "2"
    assert usable(clin, "os_event").tolist() == [True, True, True, False]


def test_loader_rejects_sample_ids_that_differ(files):
    files["data_clinical_sample.txt"].write_text(SAMPLE.replace("MB-2\tMB-2", "MB-2\tMB-2-T"))
    with pytest.raises(ValueError):
        load_metabric(files)


def test_shuffle_moves_outcomes_jointly_and_nothing_else(files):
    _, _, clin = load_metabric(files)
    s = shuffle_outcomes(clin, 1)
    assert s[["age", "subtype", "cohort"]].equals(clin[["age", "subtype", "cohort"]])
    pairs = lambda d: sorted(map(tuple, d[["time", "os_event", "dss_event"]].fillna(-1).to_numpy()))  # noqa: E731
    assert pairs(s) == pairs(clin)


def _synthetic(n=300, n_genes=60, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    X = rng.normal(size=(n_genes, n))
    z = pd.DataFrame(X, index=genes, columns=[f"P{i}" for i in range(n)])
    risk = X[:5].mean(axis=0)
    t = rng.exponential(np.exp(-risk))
    clin = pd.DataFrame({"time": t * 1000, "os_event": (rng.random(n) < 0.7).astype(float),
                         "age": rng.normal(60, 10, n), "subtype": rng.choice(["LumA", "Basal"], n),
                         "cohort": rng.choice(["1", "2", "3"], n)}, index=z.columns)
    clin["dss_event"] = clin["os_event"] * (rng.random(n) < 0.8)
    return z, pd.Series(rng.normal(8, 1, n_genes), index=genes), clin


def _spec(sigs, comparisons=()):
    return {"signatures": sigs, "comparisons": list(comparisons), "min_gene_coverage": 0.8,
            "n_floor_sets": 30, "n_perm_p3": 49, "n_boot_margin": 30}


def test_analyse_end_to_end_on_synthetic_data():
    z, gm, clin = _synthetic()
    sigs = [{"id": "fm", "model": "m", "genes": ["G0", "G1", "G2", "G3", "G4"], "direction": 1,
             "analyses": ["primary"], "p3_subtypes": ["LumA"]},
            {"id": "base", "model": "b", "genes": ["G10", "G11", "G12", "G13", "G14"], "direction": -1,
             "analyses": ["primary"]},
            {"id": "sparse", "model": "m", "genes": ["G20", "G21", "X1", "X2", "X3"], "direction": 1,
             "analyses": ["primary"]}]
    res = analyse(_spec(sigs, [{"analysis": "primary", "fm": "fm", "baseline": "base"}]), z, gm, clin)
    by = {r["id"]: r for r in res["signatures"]}
    assert not by["sparse"]["scored"] and by["sparse"]["coverage"] == pytest.approx(0.4)
    os_fm = by["fm"]["endpoints"]["OS"]
    assert os_fm["cindex_oriented"] > 0.6 and os_fm["hr_per_sd"] > 1 and os_fm["replicates"]
    assert os_fm["n"] == 300 and os_fm["n_cox"] == 300
    ob = by["base"]["endpoints"]["OS"]
    assert ob["cindex_oriented"] == pytest.approx(1 - ob["cindex"])
    assert by["fm"]["p3"]["LumA"]["OS"]["n"] == int((clin["subtype"] == "LumA").sum())
    assert by["fm"]["p3"]["LumA"]["OS"]["perm_p"] == pytest.approx(1 / 50)
    comp = [c for c in res["comparisons"] if c["endpoint"] == "OS"][0]
    assert comp["margin"] == pytest.approx(os_fm["above_floor"] - ob["above_floor"])
    assert comp["ci"][0] > 0 and comp["n_boot"] == 30


def test_analyse_on_shuffled_outcomes_does_not_replicate():
    z, gm, clin = _synthetic(seed=3)
    sig = [{"id": "fm", "model": "m", "genes": ["G0", "G1", "G2", "G3", "G4"], "direction": 1,
            "analyses": ["primary"]}]
    res = analyse(_spec(sig), z, gm, shuffle_outcomes(clin, C.SEED_REPLICATE_SHUFFLE))
    e = res["signatures"][0]["endpoints"]["OS"]
    assert abs(e["cindex_oriented"] - 0.5) < 0.06


def test_freeze_reads_picks_maximisers_leads_and_reference(tmp_path):
    import json

    from scfm.stages.replicate import freeze
    sig = lambda g: {"topk": 25, "top_cell_type": "malignant cell", "genes": g}  # noqa: E731
    sigs = {"hvg_pca": {"0.3": {"a": sig(["A", "B", "C"]), "b": sig(["D", "E", "F"])}},
            "scgpt": {"1.0": {"c": sig(["G", "H", "I"]), "d": sig(["J", "K", "L"])}}}
    scores = pd.DataFrame({"model": ["hvg_pca", "hvg_pca", "scgpt", "scgpt"],
                           "resolution": [0.3, 0.3, 1.0, 1.0], "signature": list("abcd"),
                           "cindex": [0.55, 0.40, 0.56, 0.41], "direction": [1, -1, 1, -1]})
    ladder = {"ladder": [{"model": "hvg_pca", "resolution": 0.3, "signature": "a"},
                         {"model": "scgpt", "resolution": 1.0, "signature": "c"}],
              "sensitivity": {"ladder": [{"model": "hvg_pca", "resolution": 0.3, "signature": "b"},
                                         {"model": "scgpt", "resolution": 1.0, "signature": "d"}]},
              "family_wise": {"primary": {"scgpt": {"argmax": ["1.0/c"]}},
                              "sensitivity": {"hvg_pca": {"argmax": ["0.3/b"]},
                                              "scgpt": {"argmax": ["1.0/c"]}}},
              "p3_corrected": {"picks": {"sensitivity/scgpt": {"P3_pass": True, "argmax_subtype": "LumA"},
                                         "primary/scgpt": {"P3_pass": False, "argmax_subtype": "Basal"},
                                         "sensitivity/hvg_pca": {"P3_pass": False, "borderline": True,
                                                                 "argmax_subtype": "Her2"}}}}
    nulls = {"reference": {"proliferation": {"genes_used": {"CDCA1": "NUF2", "MKI67": "MKI67"},
                                             "cindex": 0.57}}}
    paths = {}
    for name, obj in [("ladder.json", ladder), ("signatures.json", sigs), ("nulls.json", nulls)]:
        paths[name] = tmp_path / name
        paths[name].write_text(json.dumps(obj))
    paths["scores.csv"] = tmp_path / "scores.csv"
    scores.to_csv(paths["scores.csv"], index=False)
    spec = freeze(paths["ladder.json"], paths["signatures.json"], paths["scores.csv"], paths["nulls.json"])
    by = {s["id"]: s for s in spec["signatures"]}
    assert set(by) == {"hvg_pca:0.3/a", "scgpt:1.0/c", "hvg_pca:0.3/b:protective",
                       "scgpt:1.0/d:protective", "reference:proliferation"}
    assert by["scgpt:1.0/c"]["analyses"] == ["primary", "family_max_primary", "family_max_sensitivity"]
    assert by["scgpt:1.0/d:protective"]["analyses"] == ["sensitivity", "p3_lead_sensitivity"]
    assert by["scgpt:1.0/d:protective"]["p3_subtypes"] == ["LumA"]
    assert by["scgpt:1.0/d:protective"]["direction"] == -1
    assert by["hvg_pca:0.3/b:protective"]["p3_subtypes"] == ["Her2"]
    assert "p3_subtypes" not in by["scgpt:1.0/c"]
    assert by["reference:proliferation"]["genes"] == ["NUF2", "MKI67"]
    assert spec["comparisons"] == [
        {"analysis": "primary", "fm": "scgpt:1.0/c", "baseline": "hvg_pca:0.3/a"},
        {"analysis": "sensitivity", "fm": "scgpt:1.0/d:protective", "baseline": "hvg_pca:0.3/b:protective"}]
