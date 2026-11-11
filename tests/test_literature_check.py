"""Checks on the reconciliation with the rapamycin dog trial's published design inputs."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "interventions"))
import triad_published_inputs as tpi  # noqa: E402


def recorded(case_start):
    r = json.loads((ROOT / "interventions" / "triad_published_inputs.json").read_text())
    return next(c for c in r["cases"] if c["case"].startswith(case_start))


def test_recorded_reconciliation():
    # the published inputs reproduce the paper's "about one third of placebo dogs die in 3 years"
    assert 0.27 < recorded("no effect")["placebo_share_dead_in_window"] < 0.34
    assert 0.02 < recorded("no effect")["power"] < 0.08
    # a 20% risk cut is well short of 80% power in 3 years; about a 38% cut is needed
    assert 0.2 < recorded("20% risk cut for all")["power"] < 0.4
    assert 0.75 < recorded("38% risk cut")["power"] < 0.87
    # an effect that lasts only for the dosing year, or a slower rise of risk, is close to invisible
    assert recorded("20% risk cut during the dosing year")["power"] < 0.12
    assert recorded("20% slower rise")["power"] < 0.15


def test_reconciliation_runs():
    rng = np.random.default_rng(1)
    power, dead = tpi.run_case(3, {"r": 0.62}, 400, rng)
    assert 0.7 < power < 0.9 and 0.27 < dead < 0.34
    # year-1-only effect: survival is continuous at the one-year mark and returns to the normal hazard
    x0 = np.full(200000, 8.0)
    t = tpi.death_times(x0, rng, r_year1=0.8)
    h0 = tpi.A * np.exp(tpi.B * 8.0)
    expect_1y = np.exp(-0.8 * h0 * (np.exp(tpi.B) - 1) / tpi.B)
    assert abs((t > 1).mean() - expect_1y) < 0.005
