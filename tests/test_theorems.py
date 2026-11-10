"""Fast checks: every theorem in the theory paper, and the headline numbers in the saved results."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(part):
    with open(ROOT / part / "results.json") as fh:
        return json.load(fh)


def test_all_theorems_pass():
    """theory/verify_theorems.py asserts each theorem numerically or symbolically."""
    out = subprocess.run([sys.executable, "verify_theorems.py"], cwd=ROOT / "theory",
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    assert "ALL CHECKS PASSED" in out.stdout


def test_human_mortality_headline():
    h = load("human_mortality")["headline"]
    assert h["male"]["e0_period_2023"] == pytest.approx(75.79, abs=0.01)
    assert h["female"]["e0_period_2023"] == pytest.approx(81.06, abs=0.01)
    # risk of death doubles every 7 to 9 years, whichever fitting method is used
    for sex in ("male", "female"):
        for key in ("GM_doubling_time_years", "Gompertz_40_95_doubling_time_years"):
            assert 6.5 < h[sex][key] < 9.0
    # if risk stayed at its age-30 level
    assert h["male"]["no_aging_life_expectancy_hazard_frozen_at_30"] == pytest.approx(479, abs=1)


def test_dog_mortality_headline():
    d = load("dog_lifespan")["part1_mortality_fits"]["RECOMMENDED_dog_parameters"]
    assert d["b_per_year"] == pytest.approx(0.3195, abs=0.001)
    assert 2.0 < d["MRDT_years"] < 2.3
    human = d["human_comparator"]
    ratio = d["b_per_year"] / human["b_per_year"]
    assert 3.3 < ratio < 3.8          # dogs age about 3.5 times faster


def test_interventions_use_fitted_parameters():
    r = load("interventions")
    params = r["translation_math"]["mortality_parameters"]
    dog_b = load("dog_lifespan")["part1_mortality_fits"]["RECOMMENDED_dog_parameters"]["b_per_year"]
    assert params["dog"]["b"] == pytest.approx(dog_b, rel=1e-9)
    assert params["human"]["b"] == pytest.approx(0.0851, abs=0.0005)


def test_interventions_headline():
    mc = load("interventions")["monte_carlo"]
    human = mc["human"]["combo_all_overlapping_pathways"]
    dog = mc["dog"]["combo_all_overlapping_pathways"]
    assert human["mean"] == pytest.approx(4.96, abs=0.05)
    assert 1.5 < human["p5"] < 2.5 and 9.0 < human["p95"] < 10.0
    assert dog["mean"] == pytest.approx(1.18, abs=0.05)
    # everything known today is a small fraction of what stopping aging would give
    assert mc["solving_ageing_comparison"]["human"]["fraction_of_gap_closed_mean"] < 0.05


def test_escape_velocity_threshold_is_the_rate_of_aging():
    """Remaining years for a 6-year-old dog: constant-age result at rho = b exceeds half pace."""
    c = load("escape_velocity")["escape_velocity"]["cohort_remaining_le"]["dog age 6 (b=0.32)"]
    assert c["rho=0"] < c["rho=0.5b"] < c["rho=0.9b"] < c["rho=b"]
