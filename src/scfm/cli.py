"""scfm — the pipeline as a set of named stages, with a run log.

`scfm list` prints the stages; `scfm run <stage> [--dry-run] [-- args]` runs one; `scfm run all`
runs them in order, stopping at the first failure. Every real run appends to results/run_log.json:
the command, start time, duration, exit code and the params.json the stage wrote, so a run is
reconstructible from the log alone. Dry runs print each stage's plan and are not logged.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from scfm import config as C
from scfm.provenance import relpath

STAGES = [
    ("data",      "assemble the single-cell atlas and the bulk cohort; harmonise genes"),
    ("embed",     "one embedding per cell for hvg_pca, scgpt, geneformer"),
    ("states",    "Leiden clusters per representation -> marker-gene signatures"),
    ("stability", "clusters and picks under donor subsampling (needs the ladder; not part of 'all')"),
    ("translate", "score signatures in bulk; Cox; permutation null; matched-random floor; references"),
    ("ladder",    "assemble the ladder; patient bootstraps; test P1-P4, as coded and corrected"),
    ("stratify",  "cross-validated multivariable stratification: states vs age + stage"),
    ("report",    "figures, candidate shortlist, validation template"),
    ("replicate", "frozen signatures in METABRIC (run on its own; not part of 'all')"),
]
NOT_IN_ALL = {"replicate", "stability"}
SRC = Path(__file__).resolve().parents[1]


def log_path() -> Path:
    return C.RESULTS / "run_log.json"


def env() -> dict:
    e = dict(os.environ)
    e.setdefault("TORCHDYNAMO_DISABLE", "1")
    e["PYTHONPATH"] = os.pathsep.join([str(SRC)] + ([e["PYTHONPATH"]] if e.get("PYTHONPATH") else []))
    return e


def log_append(rec: dict):
    log = log_path()
    log.parent.mkdir(parents=True, exist_ok=True)
    hist = json.loads(log.read_text()) if log.exists() else []
    hist.append(rec)
    log.write_text(json.dumps(hist, indent=1))


def run_stage(name: str, extra_args: list[str], dry: bool) -> int:
    module = f"scfm.stages.{name}"
    cmd = [sys.executable, "-m", module] + extra_args + (["--dry-run"] if dry else [])
    t0 = time.time()
    print(f"\n▶ {name}  {' '.join(extra_args)}{'  (dry run)' if dry else ''}", flush=True)
    rc = subprocess.call(cmd, env=env())
    if dry:
        return rc
    args = [relpath(a) for a in extra_args]
    rec = {"stage": name, "cmd": ["python", "-m", module] + args,
           "started": t0, "seconds": round(time.time() - t0, 1), "rc": rc, "dry_run": False}
    # the params.json the stage wrote, if any (the data stage writes into data/), so the log
    # points at its provenance
    found = [p for d in (C.RESULTS, C.DATA) if d.exists() for p in d.rglob("params.json")]
    for p in sorted(found, key=os.path.getmtime, reverse=True)[:1]:
        if os.path.getmtime(p) >= t0:
            rec["params"] = relpath(str(p))
    log_append(rec)
    return rc


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    extra = []
    if "--" in argv:
        i = argv.index("--")
        argv, extra = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser(prog="scfm", description=__doc__.split("\n")[0],
                                 epilog="arguments after '--' are passed to the stage")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="print the stages")
    r = sub.add_parser("run", help="run one stage, or all of them in order")
    r.add_argument("stage", help="a stage name, or 'all'")
    r.add_argument("--dry-run", action="store_true", help="each stage prints its plan only")
    a = ap.parse_args(argv)

    if a.cmd == "list":
        print("stages:")
        for n, d in STAGES:
            print(f"  {n:<10} {d}")
        print(f"\npython: {sys.executable}\nroot:   {C.ROOT}\nlog:    {log_path()}")
        return 0
    names = [n for n, _ in STAGES if n not in NOT_IN_ALL] if a.stage == "all" else [a.stage]
    for n in names:
        if n not in dict(STAGES):
            sys.exit(f"unknown stage {n!r}; try: scfm list")
        rc = run_stage(n, extra, a.dry_run)
        if rc != 0:
            sys.exit(f"\n✗ {n} exited {rc}; stopping.")
    print("\n✓ done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
