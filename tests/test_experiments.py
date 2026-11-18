"""Fast checks on the recorded experiment results, plus a short live run of experiment 001."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXP1 = ROOT / "experiments" / "001_trial_fingerprint"
EXP2 = ROOT / "experiments" / "002_dog_cell_model"
EXP5 = ROOT / "experiments" / "005_cell_model_direct"
EXP3 = ROOT / "experiments" / "003_fading_effect"
EXP4 = ROOT / "experiments" / "004_old_dogs"
EXP6 = ROOT / "experiments" / "006_stop_and_watch"
EXP7 = ROOT / "experiments" / "007_breed_differences"
EXP8 = ROOT / "experiments" / "008_breed_levels"
EXP9 = ROOT / "experiments" / "009_partial_fade"
EXP10 = ROOT / "experiments" / "010_trial_that_settles_it"
EXP11 = ROOT / "experiments" / "011_realistic_trial"
EXP12 = ROOT / "experiments" / "012_hidden_differences_test"
EXP13 = ROOT / "experiments" / "013_breed_slopes"
EXP14 = ROOT / "experiments" / "014_health_linked_dropout"
EXP15 = ROOT / "experiments" / "015_delayed_partial_effect"
EXP16 = ROOT / "experiments" / "016_real_trial_designs"
EXP19 = ROOT / "experiments" / "019_rebound_after_stopping"


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


def test_experiment_008_recorded_result():
    r = json.loads((EXP8 / "results.json").read_text())
    assert r["gradient_check_max_relative_error"] < 1e-5
    lists = r["breed_lists"]
    assert 0.65 < lists["all 19 groups"]["spread_at_birth"] < 0.75
    assert 0.25 < lists["without the French Bulldog"]["spread_at_birth"] < 0.33
    assert 0.12 < lists["without the four shortest-lived breeds"]["spread_at_birth"] < 0.18
    for c in r["cases"]:
        if c["effect"] != "aging slowed 20%":
            assert 0.03 < c["slope_change_flagged"] < 0.08               # breed levels do not break the test
        elif c["design"] == "1,300 dogs, 10 years":
            assert c["slope_change_flagged"] > 0.8                       # and power is kept
    # widely spread hidden differences of the same nominal size still break it
    assert r["added_after_the_result"]["gamma_spread_0_3_risk_cut_50_slope_change_flagged"] > 0.15


def test_experiment_010_recorded_result():
    r = json.loads((EXP10 / "results.json").read_text())
    cell = {(c["dogs"], c["years"], c["effect"]): c["slowing_recognised"] for c in r["cases"]}
    slow20, slow30 = "aging slowed 20%", "aging slowed 30%"
    assert all(cell[(n, 4.0, slow20)] < 0.8 for n in (300, 600, 1000, 1300, 2000, 5000))   # 4 years is never enough
    assert cell[(1300, 10.0, slow20)] > 0.8 > cell[(600, 10.0, slow20)]
    assert cell[(600, 10.0, slow30)] > 0.8
    # more dogs and more years never hurt, allowing for simulation noise
    for n in (300, 600, 1000, 1300, 2000, 5000):
        assert cell[(n, 10.0, slow20)] >= cell[(n, 4.0, slow20)]


def test_experiment_011_recorded_result():
    r = json.loads((EXP11 / "results.json").read_text())
    assert 0.09 < r["frailty_variance_from_experiment_004"] < 0.11
    cell = {(c["scenario"], c["dogs"], c["effect"]): c for c in r["cases"]}
    slow, real, harsh = "aging slowed 20%", "realistic (bending, 5% dropout)", "harsh (bending, 10% dropout)"
    flag = lambda *k: cell[k]["slope_change_flagged"]
    assert flag("size mix + bending", 1300, slow) > 0.8                  # bending and size mix cost little
    assert 0.70 < flag("size mix + 10% dropout", 1300, slow) < 0.85      # dropout costs more
    assert 0.74 < flag(real, 1300, slow) < 0.85                          # 1,300 dogs sit near the line
    assert 0.68 < flag(harsh, 1300, slow) < 0.80                         # and below it with heavy dropout
    assert flag(real, 2000, slow) > 0.85 and flag(harsh, 2000, slow) > 0.85   # 2,000 dogs give a margin
    for sc in (real, harsh):
        for effect in ("no effect", "risk cut 50%"):
            assert 0.03 < flag(sc, 1300, effect) < 0.09                  # the test stays honest
        assert flag(sc, 2600, slow) >= flag(sc, 1300, slow)              # more dogs never hurt
    assert 0.35 < cell[(harsh, 1300, slow)]["share_of_dogs_dropped_out"] < 0.42


def test_experiment_012_recorded_result():
    r = json.loads((EXP12 / "results.json").read_text())
    assert r["gradient_check_max_relative_error"] < 1e-5
    cell = {(c["world"][0], c["effect"], c["analysis"]): c for c in r["cases"]}
    flag = lambda *k: cell[k]["slope_change_flagged"]
    cut, slow = "risk cut 50%", "aging slowed 20%"
    assert flag("C", cut, "plain") > 0.15 and 0.03 < flag("C", cut, "new") < 0.08      # P1: broken, then repaired
    assert cell[("C", cut, "plain")]["flagged_as_slower_aging"] < 0.02                  # misread as faster aging
    for w in "ABC":
        for effect in ("no effect", cut):
            assert 0.03 < flag(w, effect, "new") < 0.08                                 # P2: honest
    assert flag("A", slow, "plain") > 0.85                                              # P3 refuted: no cost
    assert abs(flag("A", slow, "new") - flag("A", slow, "plain")) < 0.03
    assert flag("B", slow, "new") > 0.75 and flag("C", slow, "new") > 0.75
    assert 0.08 < flag("D", cut, "plain") < 0.15 and flag("D", cut, "new") < 0.10       # P4: half and half
    assert 0.2 < cell[("C", cut, "new")]["median_spread_estimate"] < 0.4


def test_experiment_016_recorded_result():
    r = json.loads((EXP16 / "results.json").read_text())
    assert r["triad_level_per_year"] == 0.0066 and r["triad_slope_per_year"] == 0.29      # published inputs
    cell = {(c["design"], c["effect"]): c for c in r["cases"]}
    A = "rapamycin trial as published"
    slowed, fades = "aging slowed 20% in the dosing year, gain kept", "risk cut 50% in the dosing year, then fades"
    assert 0.28 < cell[(A, "no effect")]["placebo_share_dead"] < 0.32
    assert cell[(A, slowed)]["any_difference_detected"] < 0.10                           # P1 confirmed
    assert cell[(A, slowed)]["slope_change_flagged"] < 0.10
    assert 0.10 < cell[(A, fades)]["any_difference_detected"] < 0.25                     # P2 first part refuted
    assert cell[(A, fades)]["slope_change_flagged"] > 0.20                               # P2 second part confirmed
    assert cell[(A, fades)]["flagged_as_slower_aging"] < 0.02                            # misread as faster aging
    assert cell[(A, fades)]["lasting_benefit_reported"] < 0.09                           # two-period test honest
    assert 0.72 < cell[(A, "risk cut 38%, lasts all 3 years")]["any_difference_detected"] < 0.86
    for ages in ("even 10-14", "population 10-16"):                                      # P3 confirmed
        B = f"STAY-like, entry ages {ages}"
        assert 0.60 < cell[(B, "no effect")]["placebo_share_dead"] < 0.85
        assert cell[(B, "risk cut 20%")]["any_difference_detected"] > 0.85
        assert 0.22 < cell[(B, "aging slowed 20%")]["slope_change_flagged"] < 0.40
        assert 0.02 < cell[(B, "risk cut 20%")]["slope_change_flagged"] < 0.09
    C = "rapamycin trial's dogs, treated throughout, 10 years, {} dogs"
    assert 0.70 < cell[(C.format(1300), "aging slowed 20%")]["slope_change_flagged"] < 0.84   # P4: just short of 80%
    assert cell[(C.format(580), "aging slowed 20%")]["slope_change_flagged"] < 0.65
    assert 0.5 < r["side_calculations"]["stay_like_risk_ratio_that_adds_one_year"]["even 10-14"] < 0.6


def test_experiment_019_recorded_result():
    r = json.loads((EXP19 / "results.json").read_text())
    assert r["triad_level_per_year"] == 0.0066 and r["triad_slope_per_year"] == 0.29      # published inputs
    assert r["n_dogs"] == 580 and r["follow_up_years"] == 3.0 and r["dosing_years"] == 1.0
    cell = {(c["effect"], c["hidden_variance"]): c for c in r["cases"]}
    assert 0.28 < cell[("no effect", 0.0)]["placebo_share_dead"] < 0.32
    reb = cell[("lasting rebound", 0.0)]
    assert 0.20 < reb["any_difference_detected"] < 0.45                                   # P1 confirmed
    assert reb["any_difference_with_more_treated_deaths"] >= 0.9 * reb["any_difference_detected"]
    assert 0.60 < reb["reported_worse_after_stopping"] < 0.85                             # P2 confirmed
    short = cell[("short rebound", 0.0)]
    assert short["reported_worse_after_stopping"] < 0.35
    assert short["any_difference_detected"] < 0.09                                        # benefit and rebound cancel
    half = cell[("half the benefit lasts", 0.0)]
    assert 0.50 < half["any_difference_detected"] < 0.75                                  # P3 confirmed
    assert 0.20 < half["reported_lasting_benefit"] < 0.45
    for v in (0.0, 0.10, 0.30):                                                           # P4 confirmed
        for effect in ("no effect", "fades"):
            assert 0.02 < cell[(effect, v)]["after_dosing_difference_reported"] < 0.08
    side = r["survivor_risk_ratio_after_dosing"]["rapamycin design, fades (0.5 then 1)"]["hidden variance 0.3"]
    assert 1.0 < side["average_after_dosing"] < side["pooled_deaths_per_dog_year"] < 1.03
    brute = r["survivor_risk_ratio_after_dosing"]["brute-force check, 2 million dogs, rapamycin design, hidden variance 0.3"]
    assert abs(brute - side["pooled_deaths_per_dog_year"]) < 0.01


def test_experiment_014_recorded_result():
    r = json.loads((EXP14 / "results.json").read_text())
    assert r["dogs"] == 1300 and r["follow_up_years"] == 10.0 and r["background_dropout_per_year"] == 0.05
    cell = {(c["pattern"], c["effect"]): c for c in r["cases"]}
    none, slow, cut = "no effect", "aging slowed 20%", "risk cut 50%"
    for p in ("treated dogs leave more", "sick dogs leave, equally"):                     # P1 (see borderline)
        assert 0.03 < cell[(p, none)]["slope_test_fires"] < 0.09
        assert 0.03 < cell[(p, cut)]["slope_test_fires"] < 0.08
    unequal = cell[("sick treated dogs leave more", none)]                                # P2 confirmed
    assert unequal["risk_test_fires"] > 0.5
    assert unequal["risk_test_fires_pointing_to_benefit"] == unequal["risk_test_fires"]
    assert 0.78 < unequal["median_risk_ratio"] < 0.88
    assert 0.03 < unequal["slope_test_fires"] < 0.08
    assert cell[("sick untreated dogs leave more", none)]["median_risk_ratio"] > 1.12
    assert 0.06 < cell[("missed deaths drift, to 15%", none)]["slope_test_fires"] < 0.15   # P3 confirmed
    drift = cell[("missed deaths drift, to 30%", none)]
    assert 0.15 < drift["slope_test_fires"] < 0.40
    assert drift["slope_test_fires_pointing_to_slowed_aging"] >= 0.9 * drift["slope_test_fires"]
    for p in ("treated dogs leave more", "sick dogs leave, equally", "sick treated dogs leave more"):
        assert 0.65 < cell[(p, slow)]["slope_test_fires"] < 0.82                          # P4 confirmed
    assert cell[("reference", slow)]["slope_test_fires"] > 0.75
    b = json.loads((EXP14 / "borderline.json").read_text())                               # edge cases, 2,500 trials
    assert b["trials_per_case"] == 2500
    for c in b["cases"]:
        assert c["effect"] == none and 0.04 < c["slope_test_fires"] < 0.08


def test_experiment_005_recorded_result():
    r = json.loads((EXP5 / "results.json").read_text())
    old = {"UK": dict(B=0.283, C=0.230, D=0.204), "US": dict(B=0.196, C=0.168, D=0.155)}
    for label, t in r["tables"].items():
        f = t["fits"]
        assert f["A"]["rmse_log_death_rate"] > 1.0                                         # P1 confirmed
        for name in ("B", "C", "D"):                                                       # P2 confirmed
            e = f[name]["rmse_log_death_rate"]
            assert 0.12 < e and abs(e - old[label[:2]][name]) < 0.08
            assert abs(f[name]["experiment_002"]["rmse_log_death_rate"] - old[label[:2]][name]) < 0.001
            assert max(f[name]["rmse_from_each_start"]) - min(f[name]["rmse_from_each_start"]) < 0.002
        assert 0.40 < f["B"]["beta_ratio_to_human"] < 0.65                                 # P3 confirmed
        assert 2.8 < f["B"]["eta_ratio_to_human"] < 4.6
        assert t["spread_of_BCD_errors"] < 0.10 and t["order_best_to_worst_BCD"] == ["D", "C", "B"]   # P4
        res = t["residuals_fit_B_log_model_minus_log_observed"]                            # shape of the miss
        assert res[0] < -0.35 and max(res[3:6]) > 0.15
        for name in ("B", "C", "D"):                                                       # extra: ages 5 to 15
            assert t["extra_ages_5_to_15"]["fits"][name]["rmse_log_death_rate"] < 0.07
    g = r["solver_guard"]
    assert g["max_abs_z"] < 3.0 and abs(g["mean_log_ratio"]) < 0.02
    assert g["max_abs_log_gap_grid_refinement"] < 0.002
    assert -0.09 < g["experiment_002_simulation"]["mean_log_ratio"] < -0.03
    assert 0.02 < r["human_check"]["max_abs_log_gap"] < 0.05


def test_experiment_009_recorded_result():
    r = json.loads((EXP9 / "results.json").read_text())
    assert r["triad_level_per_year"] == 0.0066 and r["triad_slope_per_year"] == 0.29      # published inputs
    assert r["dosing_years"] == 1.0
    cell = {(c["design"][:2], c["effect"]): c for c in r["cases"]}
    assert len(cell) == 21
    assert cell[("D1", "no effect")]["n_dogs"] == 580 and cell[("D1", "no effect")]["follow_up_years"] == 3.0
    assert 0.28 < cell[("D1", "no effect")]["placebo_share_dead"] < 0.32
    assert cell[("D2", "no effect")]["placebo_share_dead"] > 0.95
    for d in ("D1", "D2", "D3"):                                                          # false alarms near 5%
        for effect in ("no effect", "fades at once"):
            assert 0.02 < cell[(d, effect)]["after_dosing_difference_reported"] < 0.08
    half = "half lasts"                                                                   # P1 confirmed
    assert 0.25 < cell[("D1", half)]["reported_lasting_benefit"] < 0.45
    assert cell[("D2", half)]["reported_lasting_benefit"] > 0.85
    assert cell[("D3", half)]["reported_lasting_benefit"] > 0.95
    fifth = "a fifth lasts"                                                               # P2 confirmed
    assert cell[("D1", fifth)]["reported_lasting_benefit"] < 0.15
    assert cell[("D2", fifth)]["reported_lasting_benefit"] < 0.50
    assert cell[("D3", fifth)]["reported_lasting_benefit"] < 0.50
    slow = "slow fade (halves every 2 years)"                                             # P3 confirmed
    assert 0.45 < cell[("D1", slow)]["reported_lasting_benefit"] < 0.70
    assert 0.60 < cell[("D1", slow)]["median_ratio_after_dosing"] < 0.72
    assert 0.35 < cell[("D2", slow)]["reported_lasting_benefit"] < 0.70
    assert 0.78 < cell[("D2", slow)]["median_ratio_after_dosing"] < 0.92
    quick = "quick fade (halves every 6 months)"                                          # P4: half and half
    assert cell[("D1", quick)]["reported_lasting_benefit"] < 0.20
    assert 0.80 < cell[("D3", quick)]["reported_lasting_benefit"] < 0.92                  # on the 85% line
    assert cell[("D2", quick)]["reported_lasting_benefit"] < 0.12                         # washes out in 10 years
    small = "half of a smaller benefit lasts"                                             # look the same in D1
    assert abs(cell[("D1", quick)]["median_ratio_after_dosing"] - cell[("D1", small)]["median_ratio_after_dosing"]) < 0.03
    assert cell[("D3", quick)]["true_multiplier_at_end_of_follow_up"] > 0.96


def test_experiment_013_recorded_result():
    r = json.loads((EXP13 / "results.json").read_text())
    fits = json.loads((ROOT / "dog_lifespan" / "results.json").read_text())
    assert r["dog_b_per_year"] == fits["part1_mortality_fits"]["RECOMMENDED_dog_parameters"]["b_per_year"]
    assert r["trials_per_case"] == 1000 and r["gradient_check_max_relative_error"] < 1e-5
    main = r["worlds"][0]
    assert main["breed_list"] == "without the French Bulldog" and main["share_on_rate"] == 1.0
    assert main["analysis_groups"] == 12 and len(main["breeds"]) == 18
    assert 0.27 < main["lowest_rate"] < 0.28 and 0.55 < main["highest_rate"] < 0.58
    assert 0.035 < main["rate_sd_at_entry"] < 0.045            # rates are close together among dogs at entry
    assert all(abs(b["level"] / main["reference_level"] - 1) < 1e-6 for b in main["breeds"])  # all on the rate
    cell = {(c["breed_list"], c["share_on_rate"], c["design"][-8:].strip(), c["effect"], c["analysis"]): c
            for c in r["cases"]}
    assert len(cell) == 36 and all(c["trials_analysed"] == 1000 for c in r["cases"])
    for key, c in cell.items():
        if key[3] != "aging slowed 20%":                        # P1 and P3: no true change of slope
            assert 0.03 < c["slope_change_flagged"] < 0.08, key
    m = ("without the French Bulldog", 1.0, "10 years")
    cut = "risk cut 50%"                                        # P2 refuted: under 8 percent, simple fix enough
    assert cell[m + (cut, "ignoring breed")]["slope_change_flagged"] < 0.08
    assert cell[m + (cut, "level per breed")]["slope_change_flagged"] < 0.08
    assert 1.0 < cell[m + (cut, "ignoring breed")]["median_slope_estimate"] < 1.05
    slow = "aging slowed 20%"                                   # P4 confirmed
    assert 0.75 < cell[m + (slow, "level and rate per breed")]["slope_change_flagged"] < 0.90
    four = ("without the French Bulldog", 1.0, "4 years")
    assert 0.15 < cell[four + (slow, "level and rate per breed")]["slope_change_flagged"] < 0.35
    for w in (m, ("without the French Bulldog", 0.5, "10 years"), ("all 19 groups", 1.0, "10 years")):
        best = cell[w + (slow, "level per breed")]["slope_change_flagged"]      # a level per breed has the best power
        assert best > 0.83
        assert best > cell[w + (slow, "ignoring breed")]["slope_change_flagged"] + 0.04
        assert best > cell[w + (slow, "level and rate per breed")]["slope_change_flagged"]


def test_experiment_015_recorded_result():
    r = json.loads((EXP15 / "results.json").read_text())
    r11 = json.loads((EXP11 / "results.json").read_text())
    assert r["rate_of_aging_per_year"] == r11["rate_of_aging_per_year"]          # same trial as experiment 011
    assert r["frailty_variance_from_experiment_004"] == r11["frailty_variance_from_experiment_004"]
    assert r["trials_per_case"] == 800 and r["dropout_per_year"] == 0.05
    assert r["grid_solver_max_error_years"] < 1e-5
    cell = {(c["world"], c["dogs"]): c for c in r["cases"]}
    assert len(cell) == 17
    flag = lambda w, n=1300: cell[(w, n)]["slope_change_flagged"]
    base = flag("slowed 20%, at once")
    assert 0.74 < base < 0.84                                                    # experiment 011 found 79%
    assert abs(flag("slowed 20%, builds over 1 year") - base) < 0.05             # P1 confirmed: delay costs nothing
    assert abs(flag("slowed 20%, builds over 2 years") - base) < 0.05
    assert 0.64 < flag("slowed 20%, 10% never given") < 0.76                     # P2: half and half
    assert 0.66 < flag("slowed 20%, 20% never given") < 0.72                     # above the predicted 52 to 66
    assert 0.56 < flag("slowed 20%, 30% never given") < 0.66                     # above the predicted 40 to 56
    assert 0.45 < flag("slowed 20%, 10% stop each year") < 0.52
    combined = "slowed 20%, builds over 1 year, 20% never given"
    assert flag(combined, 1600) < 0.80 <= flag(combined, 2000) < flag(combined, 2600)
    assert r["smallest_grid_size_reaching_80_percent_combined_case"] == 2000
    two = cell[("risk cut 50%, builds over 2 years", 1300)]                      # P3 confirmed: misread as slowing
    assert 0.20 < two["slope_change_flagged"] < 0.28 and 0.90 < two["median_slope_estimate"] < 0.93
    assert 0.08 < flag("risk cut 50%, builds over 1 year") < 0.12
    for w in ("no effect", "risk cut 50%, at once", "risk cut 50%, 20% never given", "risk cut 50%, 30% never given"):
        assert 0.03 < flag(w) < 0.10, w
    assert 0.03 < two["from_year_2_slope_change_flagged"] < 0.08                 # P4 confirmed: the repair works
    assert 0.03 < cell[("no effect", 1300)]["from_year_2_slope_change_flagged"] < 0.08
    assert 0.45 < cell[("slowed 20%, at once", 1300)]["from_year_2_slope_change_flagged"] < 0.58   # at a price
