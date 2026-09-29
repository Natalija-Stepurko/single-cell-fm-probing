"""Orchestrator — the pipeline as a set of callable tools, with a run log.

Each stage is a tool an agent (or a person) can call by name. `run.py list` prints the tools
with their arguments; `run.py <stage> [-- args]` runs one; `run.py all` runs them in order,
stopping at the first failure. Every call appends to results/run_log.json: the command, start
and end time, exit code, and the params.json the stage wrote — so a run is reconstructible from
the log alone.

Every analysis step has a name, an argument list, a provenance record and an exit code, so an
agent or a person can drive the pipeline step by step.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import config as C

HERE = Path(__file__).parent
VENV_PY = Path(os.environ.get("SCFM_PYTHON", "/scratch/.venv-scfm/bin/python"))
LOG = C.RESULTS / "run_log.json"

STAGES = [
    ("01_data",      "assemble the single-cell atlas and the bulk cohort; harmonise genes"),
    ("02_embed",     "one embedding per cell for hvg_pca, scgpt, geneformer"),
    ("03_states",    "Leiden clusters per representation -> marker-gene signatures"),
    ("04_translate", "score signatures in bulk; Cox; permutation null; matched-random floor"),
    ("05_ladder",    "assemble the ladder; patient bootstraps; test P1-P4"),
    ("06_report",    "figures, candidate shortlist, validation template"),
]


def env():
    e = dict(os.environ)
    # never inherit the machine-wide venv pointer; this project has its own
    e["UV_PROJECT_ENVIRONMENT"] = "/scratch/.venv-scfm"
    e.setdefault("HF_HOME", "/scratch/.hf-cache")
    e.setdefault("TORCH_HOME", "/scratch/.torch-hub")
    e.setdefault("TORCHDYNAMO_DISABLE", "1")
    return e


def log_append(rec):
    LOG.parent.mkdir(parents=True, exist_ok=True)
    hist = json.load(open(LOG)) if LOG.exists() else []
    hist.append(rec)
    json.dump(hist, open(LOG, "w"), indent=1)


def run_stage(name, extra_args, dry):
    script = HERE / f"{name}.py"
    cmd = [str(VENV_PY), str(script)] + extra_args + (["--dry-run"] if dry else [])
    t0 = time.time()
    print(f"\n▶ {name}  {' '.join(extra_args)}{'  (dry run)' if dry else ''}", flush=True)
    rc = subprocess.call(cmd, env=env())
    rec = {"stage": name, "cmd": cmd, "started": t0, "seconds": round(time.time() - t0, 1),
           "rc": rc, "dry_run": dry}
    # find the params.json the stage wrote, if any, so the log points at its provenance
    for p in sorted(C.RESULTS.rglob("params.json"), key=os.path.getmtime, reverse=True)[:1]:
        if os.path.getmtime(p) >= t0:
            rec["params"] = str(p)
    if not dry:
        log_append(rec)
    return rc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", help="'list', 'all', or a stage name")
    ap.add_argument("--dry-run", action="store_true", help="each stage prints its plan only")
    ap.add_argument("args", nargs=argparse.REMAINDER, help="passed through after '--'")
    a = ap.parse_args()
    extra = [x for x in a.args if x != "--"]

    if a.what == "list":
        print("tools:")
        for n, d in STAGES:
            print(f"  {n:<13} {d}")
        print(f"\npython: {VENV_PY}  exists={VENV_PY.exists()}\nlog:    {LOG}")
        return
    names = [n for n, _ in STAGES] if a.what == "all" else [a.what]
    for n in names:
        if n not in dict(STAGES):
            sys.exit(f"unknown stage {n!r}; try: run.py list")
        rc = run_stage(n, extra, a.dry_run)
        if rc != 0:
            sys.exit(f"\n✗ {n} exited {rc}; stopping.")
    print("\n✓ done")


if __name__ == "__main__":
    main()
