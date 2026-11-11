"""Experiment 011: the 10-year trial of experiment 010 with dropout, bending risk and a mix of sizes.

Model for one dog of size class c, entering at age x0, t years after enrolment:
    hazard(t) = z * a * exp(b_c * x0) * r * exp(b_c * s * t)
  a, b_c  size-class calibration of this repository (dog_lifespan/data/size_class_gompertz_calibration.csv,
          scenario A), read through experiment 007's code.
  z       hidden personal risk factor: gamma with mean 1 and variance V among dogs alive at age 2, so among
          dogs alive at entry age x0 it is gamma with shape 1/V and scale V / (1 + V * H_c(x0)), where
          H_c(x0) = (a / b_c) * (exp(b_c * x0) - exp(2 * b_c)). V = 0 switches it off.
          V comes from experiment 004's fit to the UK table.
          ASSUMPTION: the same V inside every size class.
  r, s    treatment: r < 1 cuts risk, s < 1 slows aging; both 1 in the control group.
  dropout each dog leaves at a constant yearly rate, same in both groups, unrelated to health.
          ASSUMPTION: 0, 5 or 10 percent a year. No published rate was opened.
ASSUMPTIONS also: entry ages uniform on 7 to 12; size mix among trial dogs = US population shares.

Analysis: slope test (likelihood ratio, 5 percent level) with a separate baseline (a and b) for each of four
size groups, as in experiment 007. It does not model z or dropout.

Run:  python3 realistic_trial.py            (about 4 minutes on 2 cores)
      python3 realistic_trial.py --quick
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


br = _load("breeds", "007_breed_differences", "breeds.py")      # size classes, stratified likelihood and fit

SEED = 20261117
AGE0 = 2.0
ENTRY = (7.0, 12.0)          # ASSUMPTION
FOLLOW = 10.0
K = 4                        # size groups in the analysis (giant analysed with large)
EFFECTS = {"aging slowed 20%": (1.0, 0.8), "no effect": (1.0, 1.0), "risk cut 50%": (0.5, 1.0)}
# scenario name: (bending on?, dropout per year)   dropout rates are ASSUMPTIONS
SCENARIOS = {
    "size mix only": (False, 0.0),
    "size mix + bending": (True, 0.0),
    "size mix + 5% dropout": (False, 0.05),
    "size mix + 10% dropout": (False, 0.10),
    "realistic (bending, 5% dropout)": (True, 0.05),
    "harsh (bending, 10% dropout)": (True, 0.10),
}
REAL, HARSH = "realistic (bending, 5% dropout)", "harsh (bending, 10% dropout)"
SLOW = "aging slowed 20%"
# (scenario, dogs, effect)
CASES = ([(sc, 1300, SLOW) for sc in SCENARIOS]
         + [(sc, n, SLOW) for sc in (REAL, HARSH) for n in (1600, 2000, 2600)]
         + [(sc, 1300, e) for sc in (REAL, HARSH) for e in ("no effect", "risk cut 50%")])


def frailty_variance():
    with open(os.path.join(HERE, "..", "004_old_dogs", "results.json")) as fh:
        return json.load(fh)["part1_uk_fit"]["with_frailty"]["frailty_variance"]


def simulate(n, r, s, a, b_class, shares, v, drop, rng):
    cls = rng.choice(len(shares), size=n, p=shares)
    x0 = rng.uniform(ENTRY[0], ENTRY[1], n)
    arm = (np.arange(n) % 2).astype(float)
    b = b_class[cls]
    if v > 0:
        H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
        z = rng.gamma(1.0 / v, v / (1.0 + v * H0))
    else:
        z = np.ones(n)
    h0 = z * a * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    leave = rng.exponential(1.0 / drop, n) if drop > 0 else np.full(n, np.inf)
    end = np.minimum(leave, FOLLOW)
    died = (t <= end).astype(float)
    return x0, arm, np.minimum(t, end), died, cls, (leave < np.minimum(t, FOLLOW))


def slope_test(data, cls, a, b0):
    """Fit 'risk cut only' and 'risk cut and slope change'; return the likelihood-ratio statistic and s-hat.
    The larger model is started from the smaller one's answer and from a second point; the better is kept."""
    lvl = br.fit("level", data, cls, K, a, b0)
    f1 = br.fit("full", data, cls, K, a, b0, extra=[lvl[1][2 * K], 0.0])
    f2 = br.fit("full", data, cls, K, a, b0, extra=[0.0, np.log(0.8)])
    full = f1 if f1[0] >= f2[0] else f2
    return 2 * (full[0] - lvl[0]), float(np.exp(full[1][-1]))


def one_case(args):
    scenario, n, effect, nsim, a, b_class, shares, v, seed = args
    bend, drop = SCENARIOS[scenario]
    r, s = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    b0 = float(np.average(b_class, weights=shares))
    flagged, s_hat, seen, left = 0, [], [], []
    for _ in range(nsim):
        *data, cls5, dropped = simulate(n, r, s, a, b_class, shares, v if bend else 0.0, drop, rng)
        stat, sh = slope_test(data, br.STRATUM[cls5], a, b0)
        flagged += stat > c1
        s_hat.append(sh)
        seen.append(float(data[3].mean()))
        left.append(float(dropped.mean()))
    return dict(scenario=scenario, dogs=n, effect=effect, n_trials=nsim, bending=bend, dropout_per_year=drop,
                slope_change_flagged=flagged / nsim, median_slope_estimate=float(np.median(s_hat)),
                share_of_deaths_seen=float(np.mean(seen)), share_of_dogs_dropped_out=float(np.mean(left)))


def main():
    quick = "--quick" in sys.argv
    nsim = 20 if quick else 800
    a, b_class, shares, names = br.size_classes()
    v = frailty_variance()
    cases = [(sc, n, e, nsim, a, b_class, shares, v, SEED + 1000 * i) for i, (sc, n, e) in enumerate(CASES)]
    if quick:
        cases = cases[:6]
    order = sorted(range(len(cases)), key=lambda i: -cases[i][1])       # big cases first, to balance the cores
    with Pool(2) as pool:
        done = pool.map(one_case, [cases[i] for i in order], chunksize=1)
    rows = [None] * len(cases)
    for i, row in zip(order, done):
        rows[i] = row
    smallest = {}
    for sc in (REAL, HARSH):
        ok = [r["dogs"] for r in rows if r["scenario"] == sc and r["effect"] == SLOW and r["slope_change_flagged"] >= 0.8]
        smallest[sc] = min(ok) if ok else None
    out = dict(
        note="Simulated 10-year trials, entry ages 7 to 12, half the dogs treated throughout. 'slope_change_flagged' "
             "is the share of trials in which the size-aware slope test reports a change in the rate of aging. "
             "ASSUMPTIONS: dropout rates (0, 5, 10 percent a year, unrelated to health); the hidden risk factor of "
             "experiment 004 applied inside every size class; entry ages; size mix equal to the US population mix. "
             "The scenario names 'realistic' and 'harsh' are labels, not measurements.",
        seed=SEED, trials_per_case=nsim, frailty_variance_from_experiment_004=v, size_classes=names,
        rate_of_aging_per_year=[float(x) for x in b_class], population_share=[float(x) for x in shares],
        common_level_a=a, smallest_grid_size_reaching_80_percent=smallest, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['scenario']:32s} | {r['dogs']:5d} | {r['effect']:17s} | flagged {r['slope_change_flagged']:.3f} | "
              f"s-hat {r['median_slope_estimate']:.2f} | deaths seen {r['share_of_deaths_seen']:.2f} | "
              f"dropped {r['share_of_dogs_dropped_out']:.2f}")
    print("smallest:", smallest)


if __name__ == "__main__":
    main()
