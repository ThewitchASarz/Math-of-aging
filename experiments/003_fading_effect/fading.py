"""Experiment 003: what can a trial detect when treatment lasts one year and then stops?

Treated dogs follow two phases (t = years since enrolment, tau = length of treatment):
    t <  tau :  hazard = h0 * r1 * exp(b * s1 * t)
    t >= tau :  hazard = h0 * exp(b * s1 * tau) * r2 * exp(b * (t - tau))
with h0 = a * exp(b * entry age). Placebo dogs have r1 = s1 = r2 = 1.

    slowed aging, benefit kept : r1 = 1, s1 < 1, r2 = 1   (the dog stays (1 - s1) * tau years "younger")
    risk cut that fades        : r1 < 1, s1 = 1, r2 = 1   (risk returns to the untreated level)
    risk cut that persists     : r1 < 1, s1 = 1, r2 = r1

The analysis is the one from experiment 001, which assumes a constant effect. Here it is deliberately
applied to data that break that assumption.

Run:  python3 fading.py            (about 3 minutes on 2 cores)
      python3 fading.py --quick
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
_spec = importlib.util.spec_from_file_location(
    "fingerprint", os.path.join(HERE, "..", "001_trial_fingerprint", "fingerprint.py"))
fp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fp)

SEED = 20261111
TAU = 1.0                                         # years of treatment
DESIGNS = {
    "580 dogs, 3 years (the rapamycin trial's shape)": (580, 7.0, 12.0, 3.0),
    "580 dogs, 10 years": (580, 7.0, 12.0, 10.0),
    "5,000 dogs, 3 years": (5000, 7.0, 12.0, 3.0),
}
# name: (r1, s1, r2)
EFFECTS = {
    "no effect": (1.0, 1.0, 1.0),
    "aging slowed 20% for 1 year, benefit kept": (1.0, 0.8, 1.0),
    "aging slowed 50% for 1 year, benefit kept": (1.0, 0.5, 1.0),
    "aging stopped for 1 year, benefit kept": (1.0, 1e-6, 1.0),
    "risk cut 30% for 1 year, then fades": (0.7, 1.0, 1.0),
    "risk cut 50% for 1 year, then fades": (0.5, 1.0, 1.0),
    "risk cut 30%, persists after treatment": (0.7, 1.0, 0.7),
}


def simulate(n, lo, hi, follow, r1, s1, r2, a, b, rng):
    x0 = rng.uniform(lo, hi, n)
    arm = (np.arange(n) % 2).astype(float)
    treated = arm == 1
    R1 = np.where(treated, r1, 1.0)
    S1 = np.where(treated, s1, 1.0)
    R2 = np.where(treated, r2, 1.0)
    h0 = a * np.exp(b * x0)
    e = -np.log(rng.random(n))                    # unit-exponential draw
    bs1 = b * S1
    H1 = h0 * R1 * np.expm1(bs1 * TAU) / bs1      # cumulative hazard over the treatment year
    t_first = np.log1p(e * bs1 / (h0 * R1)) / bs1
    h2 = h0 * np.exp(bs1 * TAU) * R2
    t_second = TAU + np.log1p(np.maximum(e - H1, 0.0) * b / h2) / b
    t = np.where(e < H1, t_first, t_second)
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def one_case(args):
    design, effect, nsim, a, b, seed = args
    n, lo, hi, follow = DESIGNS[design]
    r1, s1, r2 = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    any_level = flagged = flagged_faster = 0
    s_hat, r_hat, deaths = [], [], []
    for _ in range(nsim):
        data = simulate(n, lo, hi, follow, r1, s1, r2, a, b, rng)
        ll = fp.fit_all(data, a, b)
        any_level += 2 * (ll["level"][0] - ll["null"][0]) > c1
        flag = 2 * (ll["full"][0] - ll["level"][0]) > c1
        s = float(np.exp(ll["full"][1][3]))
        flagged += flag
        flagged_faster += flag and s > 1.0
        s_hat.append(s)
        r_hat.append(float(np.exp(ll["level"][1][2])))
        deaths.append(float(data[3].sum()))
    return dict(design=design, effect=effect, n_trials=nsim, mean_deaths=float(np.mean(deaths)),
                detected_by_usual_test=any_level / nsim, slope_change_flagged=flagged / nsim,
                flagged_as_faster_aging=flagged_faster / nsim,
                median_slope_estimate=float(np.median(s_hat)),
                median_hazard_ratio_estimate=float(np.median(r_hat)))


def main():
    quick = "--quick" in sys.argv
    nsim = 40 if quick else 400
    a, b = fp.dog_parameters()
    cases = [(d, e, nsim, a, b, SEED + 1000 * i + j)
             for i, d in enumerate(DESIGNS) for j, e in enumerate(EFFECTS)]
    if quick:
        cases = [c for c in cases if c[0].startswith("580 dogs, 3")]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases)
    out = dict(note="Simulated trials with one year of treatment. Analysed with the constant-effect models of "
                    "experiment 001. Entry ages 7 to 12 are an assumption.",
               dog_a_per_year=a, dog_b_per_year=b, treatment_years=TAU, seed=SEED, trials_per_case=nsim,
               lasting_risk_cut_from_kept_benefit={k: float(np.exp(-b * (1 - min(max(s1, 0.0), 1.0)) * TAU))
                                                   for k, (r1, s1, r2) in EFFECTS.items() if "kept" in k},
               cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['design'][:20]:20s} | {r['effect']:42s} | usual test {r['detected_by_usual_test']:.2f} | "
              f"slope flagged {r['slope_change_flagged']:.2f} (as faster {r['flagged_as_faster_aging']:.2f}) | "
              f"s-hat {r['median_slope_estimate']:.2f} | HR-hat {r['median_hazard_ratio_estimate']:.2f}")


if __name__ == "__main__":
    main()
