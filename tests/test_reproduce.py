"""Slow check: re-run every part from its code and confirm the saved results.json reproduces.

Run with:  pytest -m slow
"""
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

PARTS = {
    "human_mortality": ["run_all.py"],
    "dog_lifespan": ["run_all.py"],
    "damage_repair": ["evolution.py", "saturating_removal.py", "clocks.py", "cancer.py", "resilience.py"],
    "interventions": ["translation_math.py", "evidence_tables.py", "monte_carlo.py", "trial_design.py"],
    "escape_velocity": ["run_all.py"],
}
SKIP_KEYS = ("runtime", "seconds", "elapsed", "time_s")


def differences(old, new, path=""):
    """Yield paths where two JSON trees differ beyond floating-point noise."""
    if isinstance(old, dict) and isinstance(new, dict):
        for key in sorted(set(old) | set(new)):
            if any(s in str(key).lower() for s in SKIP_KEYS):
                continue
            if key not in old or key not in new:
                yield f"{path}/{key} (missing on one side)"
            else:
                yield from differences(old[key], new[key], f"{path}/{key}")
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            yield f"{path} (length {len(old)} vs {len(new)})"
        else:
            for i, (a, b) in enumerate(zip(old, new)):
                yield from differences(a, b, f"{path}[{i}]")
    elif isinstance(old, (int, float)) and isinstance(new, (int, float)) \
            and not isinstance(old, bool) and not isinstance(new, bool):
        if math.isnan(old) and math.isnan(new):
            return
        if not math.isclose(old, new, rel_tol=1e-6, abs_tol=1e-9):
            yield f"{path}: {old} vs {new}"
    elif old != new:
        yield f"{path}: {old!r} vs {new!r}"


@pytest.mark.slow
@pytest.mark.parametrize("part", list(PARTS))
def test_part_reproduces(part, tmp_path):
    saved_path = ROOT / part / "results.json"
    backup = tmp_path / "results.json"
    shutil.copy(saved_path, backup)
    with open(backup) as fh:
        saved = json.load(fh)
    try:
        for script in PARTS[part]:
            out = subprocess.run([sys.executable, script], cwd=ROOT / part,
                                 capture_output=True, text=True, timeout=1200)
            assert out.returncode == 0, f"{part}/{script} failed:\n{out.stderr[-2000:]}"
        with open(saved_path) as fh:
            fresh = json.load(fh)
    finally:
        shutil.copy(backup, saved_path)
    diffs = list(differences(saved, fresh))
    assert not diffs, f"{part}: {len(diffs)} values changed, e.g. {diffs[:5]}"
