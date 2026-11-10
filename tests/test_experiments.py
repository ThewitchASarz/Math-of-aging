"""Fast checks on the recorded experiment results, plus a short live run of experiment 001."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXP1 = ROOT / "experiments" / "001_trial_fingerprint"
EXP2 = ROOT / "experiments" / "002_dog_cell_model"
EXP3 = ROOT / "experiments" / "003_fading_effect"
EXP4 = ROOT / "experiments" / "004_old_dogs"
EXP6 = ROOT / "experiments" / "006_stop_and_watch"
EXP7 = ROOT / "experiments" / "007_breed_differences"


def case(results, design_start, effect):
    for c in results["cases"]:
        if c["design"].startswith(design_start) and c["effect"] == effect:
            return c
    raise KeyError((design_start, effect))


def test_experiment_001_recorded_result():
    r = json.loads((EXP1 / "results.json").read_text())
    assert r["gradient_check_max_relative_error"] < 1e-5
    assert max(c["fits_not_converged"] for c in r["cases"]) < 0.02
    # P1, P2: short trials rarely recognise a true 20% slowing of aging
    assert case(r, "580 dogs, 3 years (like", "aging slowed 20%")["slope_change_flagged"] < 0.30
    assert case(r, "1,300 dogs, 4 years", "aging slowed 20%")["slope_change_flagged"] < 0.50
    # P3: length beats size
    long_ = case(r, "1,300 dogs, 10 years", "aging slowed 20%")["slope_change_flagged"]
    big = case(r, "5,000 dogs, 3 years", "aging slowed 20%")["slope_change_flagged"]
    assert long_ > 0.8 > big
    # P4: a pure risk cut is not misread as a slope change
    false_alarms = [c["slope_change_flagged"] for c in r["cases"] if c["kind"] in ("none", "level")]
    assert 0.02 < min(false_alarms) and max(false_alarms) < 0.09


def test_experiment_001_runs():
    out = subprocess.run([sys.executable, "fingerprint.py", "--quick"], cwd=EXP1,
                         capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr[-2000:]
    quick = EXP1 / "results_quick.json"
    r = json.loads(quick.read_text())
    quick.unlink()
    assert len(r["cases"]) == 21


def test_experiment_002_recorded_result():
    r = json.loads((EXP2 / "results.json").read_text())
    for table in r["tables"].values():
        fits = table["fits"]
        assert fits["A"]["rmse_log_death_rate"] > 1.5                  # production alone fails
        for name in ("B", "C", "D"):
            assert fits[name]["rmse_log_death_rate"] < 0.3             # one more setting fixes most of it
            assert fits[name]["rmse_log_death_rate"] > 0.15            # but not to within 15% (P2 refuted)
        assert fits["B"]["beta_ratio_to_human"] < 0.6                   # dogs clear less in that fit
    check = json.loads((EXP2 / "sde_check.json").read_text())
    assert 0.1 < check["max_abs_log_error"] < 0.3                       # the shortcut is only approximate


def test_experiment_003_recorded_result():
    r = json.loads((EXP3 / "results.json").read_text())
    short = "580 dogs, 3 years"
    assert case(r, short, "aging slowed 20% for 1 year, benefit kept")["detected_by_usual_test"] < 0.10
    assert 0.40 < case(r, short, "aging stopped for 1 year, benefit kept")["detected_by_usual_test"] < 0.75
    fade = case(r, short, "risk cut 30% for 1 year, then fades")
    assert fade["slope_change_flagged"] > 0.08
    assert fade["flagged_as_faster_aging"] > 0.9 * fade["slope_change_flagged"]   # read as faster aging
    assert case(r, "5,000 dogs", "risk cut 30% for 1 year, then fades")["slope_change_flagged"] > 0.25


def test_experiment_004_recorded_result():
    r = json.loads((EXP4 / "results.json").read_text())
    p1 = r["part1_uk_fit"]
    assert p1["likelihood_ratio"] > 50 and 0.05 < p1["with_frailty"]["frailty_variance"] < 0.2
    ratio = r["part2_forecast_old_ages"]["forecast_over_observed"]
    assert ratio["14"] < 1.1 and ratio["15"] < 1.1          # the doubling law holds to about 15
    assert ratio["16"] > 1.15 and ratio["18"] > 1.5         # and overstates risk after that
    rows = {(c["design"], c["effect"]): c for c in r["part4_trials_with_frailty"]["cases"]}
    assert rows[("1,300 dogs, 10 years", "risk cut 30%")]["slope_change_flagged"] < 0.10
    assert rows[("1,300 dogs, 10 years", "aging slowed 20%")]["slope_change_flagged"] > 0.8


def test_experiment_006_recorded_result():
    r = json.loads((EXP6 / "results.json").read_text())
    assert r["gradient_check_max_relative_error"] < 1e-5
    for c in r["cases"]:
        if "fades" in c["effect"] or c["effect"] == "no effect":
            assert 0.02 < c["lasting_benefit_reported"] < 0.09          # not fooled by a fading benefit
    short = "580 dogs, 3 years"
    assert case(r, short, "risk cut 30%, persists after treatment")["lasting_benefit_reported"] > 0.5
    assert case(r, short, "aging slowed 20% for 1 year, benefit kept")["lasting_benefit_reported"] < 0.10
    assert case(r, "580 dogs, 10 years", "aging stopped for 1 year, benefit kept")["lasting_benefit_reported"] > 0.85


def test_experiment_007_recorded_result():
    r = json.loads((EXP7 / "results.json").read_text())
    p1 = r["part1_spread_of_risk"]
    assert p1["variance_between_breeds_uk"] > 0.2 > p1["variance_between_size_classes_us"]
    assert p1["shortest_lived"][0] == "French Bulldog"
    rows = {(c["design"], c["effect"], c["analysis"]): c for c in r["part2_trials_with_mixed_sizes"]["cases"]}
    for (design, effect, analysis), c in rows.items():
        if analysis == "size-aware" and effect != "aging slowed 20%":
            assert 0.02 < c["slope_change_flagged"] < 0.09               # size-aware analysis keeps its size
    assert rows[("1,300 dogs, 10 years", "aging slowed 20%", "size-aware")]["slope_change_flagged"] > 0.8
    assert rows[("1,300 dogs, 10 years", "risk cut 50%", "ignoring size")]["slope_change_flagged"] < 0.10
