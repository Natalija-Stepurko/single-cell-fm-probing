"""Provenance: every stage writes params.json beside its outputs, recording what produced them."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

from scfm import config as C

VERSIONED = ("numpy", "scipy", "pandas", "sklearn", "scanpy", "anndata", "lifelines", "leidenalg",
             "igraph", "joblib", "xgboost", "matplotlib", "torch", "transformers", "huggingface_hub",
             "cellxgene_census")


def _version(mod: str) -> str | None:
    """Installed version of an importable module, read from package metadata without importing it."""
    from importlib.metadata import PackageNotFoundError, version
    from importlib.util import find_spec
    dist = {"sklearn": "scikit-learn", "huggingface_hub": "huggingface-hub",
            "cellxgene_census": "cellxgene-census", "igraph": "igraph"}.get(mod, mod)
    try:
        if find_spec(mod) is None:
            return None
        return version(dist)
    except (ImportError, ValueError, PackageNotFoundError):
        return None


def relpath(s: str) -> str:
    """A path under ROOT (or the checkout), written relative to it so the record is portable."""
    if not s.startswith(("/", ".")):
        return s
    p = Path(os.path.abspath(s))
    for base in (C.ROOT, C.REPO):
        for b in {Path(os.path.abspath(base)), base.resolve()}:
            if p.is_relative_to(b):
                return str(p.relative_to(b))
    return s


def _jsonable(v):
    if isinstance(v, Path):
        return relpath(str(v))
    if isinstance(v, str):
        return relpath(v)
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (int, float, bool)) or v is None:
        return v
    return repr(v)


def _git(*cmd, default=None):
    try:
        r = subprocess.run(("git", "-C", str(Path(__file__).resolve().parent)) + cmd,
                           capture_output=True, text=True, timeout=10)
        return r.stdout.strip() if r.returncode == 0 else default
    except Exception:
        return default


def record_params(out_dir, args=None, extra=None, filename="params.json") -> dict:
    """Write `params.json` beside a stage's outputs.

    Captures the resolved argument namespace, the command line, the git commit (flagged dirty if
    the tree has uncommitted changes), the versions of libraries whose numerics affect results,
    the pinned model revisions and a UTC timestamp. Paths under ROOT are written relative to it.

    Args:
        out_dir:  directory the stage writes into; created if absent.
        args:     the argparse namespace (or any object with __dict__), optional.
        extra:    dict of stage-specific facts worth pinning.
        filename: override when several stages share one output directory.

    Returns the dict it wrote.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    versions = {"python": sys.version.split()[0]}
    for mod in VERSIONED:
        v = _version(mod)
        if v is not None:
            versions[mod] = v

    payload = {
        "script": Path(sys.argv[0]).name,
        "timestamp_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "command": " ".join(_jsonable(list(sys.argv))),
        "args": {k: _jsonable(v) for k, v in vars(args).items()} if args is not None else {},
        "git_commit": _git("rev-parse", "HEAD", default="unknown"),
        "git_dirty": bool(_git("status", "--porcelain", default="")),
        "versions": versions,
        "model_revisions": {m: s["revision"] for m, s in C.MODELS.items() if "revision" in s},
    }
    if extra:
        payload["extra"] = {k: _jsonable(v) for k, v in extra.items()}

    (out / filename).write_text(json.dumps(payload, indent=2) + "\n")
    return payload
