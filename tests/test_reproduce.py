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
    "experiments/001_trial_fingerprint": ["fingerprint.py"],
    "experiments/002_dog_cell_model": ["dog_cell_model.py"],
    "experiments/003_fading_effect": ["fading.py"],
    "experiments/004_old_dogs": ["old_dogs.py"],
    "experiments/006_stop_and_watch": ["two_period.py"],
    "experiments/007_breed_differences": ["breeds.py"],
    "experiments/008_breed_levels": ["breed_levels.py"],
    "experiments/010_trial_that_settles_it": ["design_grid.py"],
    "experiments/011_realistic_trial": ["realistic_trial.py"],
}
# Experiments count outcomes over simulated trials and stop optimisers at a tolerance, so a different
# machine can move a count by one or two trials. They get a looser comparison.
# damage_repair contains one three-setting refit that is poorly pinned down by the data; its optimiser
# stops at slightly different points on different processors (seen: 4 parts in a million).
LOOSE = {"damage_repair": (1e-4, 1e-9),
         "experiments/001_trial_fingerprint": (0.05, 0.011), "experiments/002_dog_cell_model": (0.02, 1e-3),
         "experiments/003_fading_effect": (0.05, 0.011), "experiments/004_old_dogs": (0.05, 0.011),
         "experiments/006_stop_and_watch": (0.05, 0.011), "experiments/007_breed_differences": (0.05, 0.011),
         "experiments/008_breed_levels": (0.05, 0.011), "experiments/010_trial_that_settles_it": (0.05, 0.011),
         "experiments/011_realistic_trial": (0.05, 0.011)}
SKIP_KEYS = ("runtime", "seconds", "elapsed", "time_s")


def differences(old, new, path="", rel=1e-6, abs_=1e-9):
    """Yield paths where two JSON trees differ beyond floating-point noise."""
    if isinstance(old, dict) and isinstance(new, dict):
        for key in sorted(set(old) | set(new)):
            if any(s in str(key).lower() for s in SKIP_KEYS):
                continue
            if key not in old or key not in new:
                yield f"{path}/{key} (missing on one side)"
            else:
                yield from differences(old[key], new[key], f"{path}/{key}", rel, abs_)
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            yield f"{path} (length {len(old)} vs {len(new)})"
        else:
            for i, (a, b) in enumerate(zip(old, new)):
                yield from differences(a, b, f"{path}[{i}]", rel, abs_)
    elif isinstance(old, (int, float)) and isinstance(new, (int, float)) \
            and not isinstance(old, bool) and not isinstance(new, bool):
        if math.isnan(old) and math.isnan(new):
            return
        if not math.isclose(old, new, rel_tol=rel, abs_tol=abs_):
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
    rel, abs_ = LOOSE.get(part, (1e-6, 1e-9))
    diffs = list(differences(saved, fresh, rel=rel, abs_=abs_))
    assert not diffs, f"{part}: {len(diffs)} values changed, e.g. {diffs[:5]}"
