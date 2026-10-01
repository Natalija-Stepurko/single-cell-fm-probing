import os
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"


def scfm(*args, root):
    env = {k: v for k, v in os.environ.items() if k not in ("SCFM_SMOKE",)}
    env["SCFM_ROOT"] = str(root)
    env["PYTHONPATH"] = os.pathsep.join([str(SRC)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    return subprocess.run([sys.executable, "-m", "scfm", *args], env=env, capture_output=True, text=True,
                          timeout=300)


def test_list(tmp_path):
    r = scfm("list", root=tmp_path)
    assert r.returncode == 0, r.stderr
    for stage in ("data", "embed", "states", "translate", "ladder", "report"):
        assert stage in r.stdout
    assert str(tmp_path) in r.stdout


def test_run_all_dry_run_needs_no_data(tmp_path):
    r = scfm("run", "all", "--dry-run", root=tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "done" in r.stdout
    assert not (tmp_path / "results" / "run_log.json").exists()
    assert not list(tmp_path.rglob("params.json"))


def test_unknown_stage_fails(tmp_path):
    r = scfm("run", "nope", "--dry-run", root=tmp_path)
    assert r.returncode != 0
