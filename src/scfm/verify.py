"""Compare regenerated results with results/MANIFEST.sha256, or rewrite the manifest.

The manifest covers every tracked file under results/ except params.json and run_log.json, which
carry timestamps and commits. It lives in the checkout; the files it is checked against live under
ROOT, so `SCFM_ROOT=<elsewhere> python -m scfm.verify` checks a separate run. Only a fresh atlas
build writes results/atlas_cells.csv, so a run that starts later (stage 03, or `make reproduce`
into an empty root) reports it missing; pass --skip-missing for such partial reruns.
"""
import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

from scfm import config as C

MANIFEST = C.REPO / "results" / "MANIFEST.sha256"
EXCLUDE = {"params.json", "run_log.json", MANIFEST.name}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tracked_results() -> list[str]:
    out = subprocess.run(["git", "-C", str(C.REPO), "ls-files", "results"],
                         capture_output=True, text=True, check=True).stdout.split()
    return sorted(p for p in out if Path(p).name not in EXCLUDE)


def read_manifest() -> dict[str, str]:
    rows = [ln.split(maxsplit=1) for ln in MANIFEST.read_text().splitlines() if ln.strip()]
    return {path: digest for digest, path in rows}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--write", action="store_true",
                    help="rewrite the manifest from the tracked files in the checkout")
    ap.add_argument("--skip-missing", action="store_true",
                    help="ignore manifest entries with no file under ROOT (a partial rerun)")
    a = ap.parse_args(argv)

    if a.write:
        files = tracked_results()
        MANIFEST.write_text("".join(f"{sha256(C.REPO / p)}  {p}\n" for p in files))
        print(f"wrote {len(files)} entries to {MANIFEST}")
        return 0

    bad = 0
    for path, digest in read_manifest().items():
        f = C.ROOT / path
        if not f.exists():
            status = "missing"
            bad += not a.skip_missing
        elif sha256(f) == digest:
            status = "ok"
        else:
            status = "MISMATCH"
            bad += 1
        print(f"{status:<9} {path}")
    print(f"\n{'all match' if not bad else f'{bad} file(s) differ or are missing'} (root: {C.ROOT})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
