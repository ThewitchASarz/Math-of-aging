"""Experiment 014: dropout that depends on health or on treatment.

The trial is experiment 011's "realistic" one: 1,300 dogs, half treated throughout, entry ages 7 to 12,
10 years, five size classes, a hidden personal risk factor z (gamma, variance V = 0.102 among dogs alive at
age 2, from experiment 004), and a background dropout of 5 percent a year.

Model for one dog of size class c, entering at age x0, t years after enrolment:
    hazard(t) = z * a * exp(b_c * x0) * r * exp(b_c * s * t)        (r = s = 1 if untreated)

Dropout patterns (ALL ARE ASSUMPTIONS; no published dropout pattern for a dog trial was found):
    rate    yearly rate of leaving, unrelated to health:  rate_control, rate_treated
    frail   extra yearly rate for treated dogs = frail * z   (frail dogs leave the treated group more)
    sick    share of dogs about to die that are withdrawn first, at death time minus Uniform(0, LAG) years
            (not before enrolment): sick_control, sick_treated
    drift   treated dogs only: the share withdrawn before death is drift * (death time / 10 years)
A withdrawn dog is censored: it counts as alive up to the day it left, and its death is not recorded.

Analysis (does not model z or any dropout): maximum likelihood with a baseline (a, b) for each of four size
groups, from experiment 007's code.
    slope test: "risk cut and slope change" against "risk cut only"   (likelihood ratio, 5 percent level)
    risk test:  "risk cut only" against "no difference"

Run:  python3 dropout.py            (about 5 minutes on 2 cores)
      python3 dropout.py --quick
"""
import importlib.util
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from multiprocessing import Pool

import numpy as np
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, folder, script):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", folder, script))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


br = _load("breeds", "007_breed_differences", "breeds.py")          # size classes, stratified likelihood, fit
rt = _load("realistic", "011_realistic_trial", "realistic_trial.py")  # frailty variance, constants

SEED = 20261118
N = 1300
AGE0, ENTRY, FOLLOW, K = rt.AGE0, rt.ENTRY, rt.FOLLOW, rt.K
BACKGROUND = 0.05            # ASSUMPTION, as in experiment 011's realistic scenario
LAG = 0.5                    # ASSUMPTION: sick dogs are withdrawn up to 6 months before death
EFFECTS = {"no effect": (1.0, 1.0), "risk cut 50%": (0.5, 1.0), "aging slowed 20%": (1.0, 0.8)}
# name: dict(rate=(control, treated), frail=..., sick=(control, treated), drift=...)   ALL ASSUMPTIONS
PATTERNS = {
    "reference": dict(),
    "treated dogs leave more": dict(rate=(0.05, 0.10)),
    "sick dogs leave, equally": dict(sick=(0.30, 0.30)),
    "sick treated dogs leave more": dict(sick=(0.15, 0.30)),
    "sick untreated dogs leave more": dict(sick=(0.30, 0.15)),
    "frail treated dogs leave": dict(frail=0.05),
    "missed deaths drift, to 15%": dict(drift=0.15),
    "missed deaths drift, to 30%": dict(drift=0.30),
}
P = list(PATTERNS)
CASES = ([(p, "no effect") for p in P]
         + [(P[i], "aging slowed 20%") for i in (0, 1, 2, 3, 5)]
         + [(P[i], "risk cut 50%") for i in (1, 2, 3, 5)])


def simulate(n, r, s, a, b_class, shares, v, pattern, rng):
    rate = pattern.get("rate", (BACKGROUND, BACKGROUND))
    sick = pattern.get("sick", (0.0, 0.0))
    frail, drift = pattern.get("frail", 0.0), pattern.get("drift", 0.0)
    window = pattern.get("window", LAG)
    cls = rng.choice(len(shares), size=n, p=shares)
    x0 = rng.uniform(ENTRY[0], ENTRY[1], n)
    arm = (np.arange(n) % 2).astype(float)
    b = b_class[cls]
    H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
    z = rng.gamma(1.0 / v, v / (1.0 + v * H0))               # hidden factor among dogs alive at entry
    h0 = z * a * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    death = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    leave_rate = np.where(arm == 1, rate[1] + frail * z, rate[0])
    leave = rng.exponential(1.0, n) / leave_rate
    p_sick = np.where(arm == 1, sick[1] + drift * np.minimum(death, FOLLOW) / FOLLOW, sick[0])
    withdrawn = rng.random(n) < p_sick
    lag = rng.uniform(0.0, window, n)
    leave = np.where(withdrawn, np.minimum(leave, np.maximum(death - lag, 0.0)), leave)
    end = np.minimum(leave, FOLLOW)
    died = (death <= end).astype(float)
    t = np.minimum(death, end)
    missed = (death <= FOLLOW) & (died == 0)                 # died inside the 10 years but not recorded
    return x0, arm, t, died, cls, missed, (leave < np.minimum(death, FOLLOW))


def analyse(data, cls, a, b0):
    """Returns the two likelihood-ratio statistics, the risk ratio (slope held equal) and the slope ratio."""
    nul = br.fit("null", data, cls, K, a, b0)
    lvl = br.fit("level", data, cls, K, a, b0)
    f1 = br.fit("full", data, cls, K, a, b0, extra=[lvl[1][2 * K], 0.0])
    f2 = br.fit("full", data, cls, K, a, b0, extra=[0.0, np.log(0.8)])
    full = f1 if f1[0] >= f2[0] else f2
    return 2 * (lvl[0] - nul[0]), 2 * (full[0] - lvl[0]), float(np.exp(lvl[1][2 * K])), float(np.exp(full[1][-1]))


def one_case(args):
    name, effect, nsim, a, b_class, shares, v, seed = args
    r, s = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    b0 = float(np.average(b_class, weights=shares))
    risk_flag = risk_benefit = slope_flag = slope_slow = 0
    r_hat, s_hat, seen, missed_c, missed_t, left_c, left_t = [], [], [], [], [], [], []
    for _ in range(nsim):
        *data, cls5, missed, left = simulate(N, r, s, a, b_class, shares, v, PATTERNS[name], rng)
        st_r, st_s, rh, sh = analyse(data, br.STRATUM[cls5], a, b0)
        risk_flag += st_r > c1
        risk_benefit += (st_r > c1) and (rh < 1)
        slope_flag += st_s > c1
        slope_slow += (st_s > c1) and (sh < 1)
        r_hat.append(rh)
        s_hat.append(sh)
        arm = data[1]
        seen.append(float(data[3].mean()))
        missed_c.append(float(missed[arm == 0].mean()))
        missed_t.append(float(missed[arm == 1].mean()))
        left_c.append(float(left[arm == 0].mean()))
        left_t.append(float(left[arm == 1].mean()))
    return dict(pattern=name, effect=effect, n_trials=nsim, dogs=N,
                slope_test_fires=slope_flag / nsim, slope_test_fires_pointing_to_slowed_aging=slope_slow / nsim,
                risk_test_fires=risk_flag / nsim, risk_test_fires_pointing_to_benefit=risk_benefit / nsim,
                median_slope_ratio=float(np.median(s_hat)), median_risk_ratio=float(np.median(r_hat)),
                share_of_dogs_with_death_recorded=float(np.mean(seen)),
                share_left_untreated=float(np.mean(left_c)), share_left_treated=float(np.mean(left_t)),
                share_died_unrecorded_untreated=float(np.mean(missed_c)),
                share_died_unrecorded_treated=float(np.mean(missed_t)))


def main():
    quick = "--quick" in sys.argv
    nsim = 20 if quick else 800
    a, b_class, shares, names = br.size_classes()
    v = rt.frailty_variance()
    cases = [(p, e, nsim, a, b_class, shares, v, SEED + 1000 * i) for i, (p, e) in enumerate(CASES)]
    if quick:
        cases = cases[:6]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases, chunksize=1)
    out = dict(
        note="Simulated 10-year trials of 1,300 dogs, entry ages 7 to 12, half treated throughout, size mix and "
             "hidden risk factor as in experiment 011's realistic scenario. Each number is a share of simulated "
             "trials. ALL DROPOUT PATTERNS ARE ASSUMPTIONS: background 5 percent a year; the 10 percent rate; the "
             "15 and 30 percent shares of sick dogs withdrawn; the 6-month window; the frail-dog rate; the drifts.",
        seed=SEED, trials_per_case=nsim, dogs=N, follow_up_years=FOLLOW, entry_ages=list(ENTRY),
        background_dropout_per_year=BACKGROUND, withdrawal_window_years=LAG,
        frailty_variance_from_experiment_004=v, size_classes=names,
        rate_of_aging_per_year=[float(x) for x in b_class], population_share=[float(x) for x in shares],
        common_level_a=a, patterns={k: {kk: list(vv) if isinstance(vv, tuple) else vv for kk, vv in d.items()}
                                    for k, d in PATTERNS.items()}, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['pattern']:31s} | {r['effect']:16s} | slope {r['slope_test_fires']:.3f} "
              f"(slow {r['slope_test_fires_pointing_to_slowed_aging']:.3f}) s-hat {r['median_slope_ratio']:.3f} | "
              f"risk {r['risk_test_fires']:.3f} (benefit {r['risk_test_fires_pointing_to_benefit']:.3f}) "
              f"r-hat {r['median_risk_ratio']:.3f} | seen {r['share_of_dogs_with_death_recorded']:.2f} | "
              f"left {r['share_left_untreated']:.2f}/{r['share_left_treated']:.2f} | "
              f"missed {r['share_died_unrecorded_untreated']:.2f}/{r['share_died_unrecorded_treated']:.2f}")


if __name__ == "__main__":
    main()
