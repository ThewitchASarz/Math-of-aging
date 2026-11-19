"""Experiment 037: the opposite direction (frail dogs helped more) and a weaker treatment.

Model for one dog entering at age x0, t years after enrolment (experiment 033's, code reused unchanged):
    hazard(t) = z * a * exp(b * x0) * r * exp(b * t)
  z   hidden personal risk factor: gamma with mean 1 and variance V = 0.3 (ASSUMPTION) among dogs alive at age 2.
  r   the treatment's risk multiplier for that dog (1 in the control group): 1 - c_sturdy for the half of dogs
      with lower hidden risk among dogs of the same entry age, 1 - c_frail for the other half.
      Every cut size is an ASSUMPTION. No dog's rate of aging changes in any case.
  a, b  the UK fit in dog_lifespan/results.json.

Design: 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years, one size, no dropout.

Analyses: the plain slope test and the hidden-differences slope test, imported unchanged from experiment 017
(tests_one_size). Each is a 5% likelihood-ratio test of "risk cut only" against "risk cut and slope change".

Run:  python3 opposite_weaker.py            (about 4 minutes on 2 cores)
      python3 opposite_weaker.py --quick
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


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


e33 = _load("benefit_size", "033_benefit_difference_size", "benefit_size.py")
h17, fp = e33.h17, e33.fp

SEED = 20261137
V = 0.3                                   # ASSUMPTION: spread of hidden risk
LONG = e33.LONG                           # entry ages 7 to 12: ASSUMPTION; 10 years
# name: (average cut, cut in the sturdier half, cut in the frailer half)   -- all ASSUMPTIONS
CASES = {
    "average cut 50 | 45 and 55": (0.50, 0.45, 0.55),
    "average cut 50 | 40 and 60": (0.50, 0.40, 0.60),
    "average cut 50 | 35 and 65": (0.50, 0.35, 0.65),
    "average cut 50 | 30 and 70": (0.50, 0.30, 0.70),
    "average cut 30 | 30 and 30": (0.30, 0.30, 0.30),
    "average cut 30 | 36 and 24": (0.30, 0.36, 0.24),
    "average cut 30 | 42 and 18": (0.30, 0.42, 0.18),
    "average cut 30 | 18 and 42": (0.30, 0.18, 0.42),
    "average cut 20 | 20 and 20": (0.20, 0.20, 0.20),
    "average cut 20 | 24 and 16": (0.20, 0.24, 0.16),
    "average cut 20 | 28 and 12": (0.20, 0.28, 0.12),
    "average cut 20 | 12 and 28": (0.20, 0.12, 0.28),
}


def one_chunk(args):
    name, nsim, a, b, seed = args
    _, c_sturdy, c_frail = CASES[name]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    out = {k: dict(flag=0, slower=0, s=[]) for k in ("plain", "new")}
    dead = 0.0
    for _ in range(nsim):
        data, _cls = e33.simulate(LONG, False, V, "halves", c_sturdy, c_frail, a, b, None, None, rng)
        dead += data[3].mean()
        for label, (stat, s_hat, _) in h17.tests_one_size(data, a, b).items():
            out[label]["flag"] += int(stat > c1)
            out[label]["slower"] += int(stat > c1 and s_hat < 1.0)
            out[label]["s"].append(s_hat)
    return name, nsim, dead, out


def main():
    quick = "--quick" in sys.argv
    chunks, per = (1, 10) if quick else (2, 300)
    a, b = fp.dog_parameters()
    for name, (avg, cs, cf) in CASES.items():
        assert abs((cs + cf) / 2 - avg) < 1e-12, name
    jobs = [(name, per, a, b, SEED + 1000 * i + k) for i, name in enumerate(CASES) for k in range(chunks)]
    if quick:
        jobs = jobs[:3]
    with Pool(2) as pool:
        parts = pool.map(one_chunk, jobs, chunksize=1)
    rows = []
    for name in CASES:
        mine = [p for p in parts if p[0] == name]
        if not mine:
            continue
        n = sum(p[1] for p in mine)
        avg, cs, cf = CASES[name]
        for label in ("plain", "new"):
            s_all = [x for p in mine for x in p[3][label]["s"]]
            flag = sum(p[3][label]["flag"] for p in mine) / n
            slower = sum(p[3][label]["slower"] for p in mine) / n
            rows.append(dict(case=name, analysis=label, n_trials=n, average_cut=avg, cut_sturdier_half=cs,
                             cut_frailer_half=cf,
                             helped_more="frail" if cf > cs else ("sturdy" if cs > cf else "neither"),
                             share_of_dogs_dead=sum(p[2] for p in mine) / n,
                             slope_change_flagged=flag, flagged_as_slower_aging=slower,
                             flagged_as_faster_aging=flag - slower,
                             median_slope_estimate=float(np.median(s_all))))
    out = dict(note="Simulated 10-year trials of 1,300 dogs, half treated from enrolment, entry ages 7 to 12. No "
                    "dog's rate of aging changes in any case; the treatment cuts risk by one fraction in the "
                    "sturdier half of dogs and another in the frailer half. 'flagged_as_faster_aging' "
                    "('..slower..') is the share of trials in which a 5% likelihood-ratio test rejects 'the slope "
                    "did not change' with an estimated rate above (below) 1. 'plain' ignores hidden differences; "
                    "'new' adds one setting for the spread of a gamma-shaped hidden risk factor. ASSUMPTIONS: the "
                    "spread of hidden risk (0.3), entry ages, and every cut.",
               seed=SEED, uk_a_per_year=a, uk_b_per_year=b, spread=V, years=LONG[2], entry_ages=list(LONG[:2]),
               n_dogs=e33.N, trials_per_case=chunks * per, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case']:28s} | {r['analysis']:5s} | dead {r['share_of_dogs_dead']:.3f}"
              f" | slower {r['flagged_as_slower_aging']:.3f} | faster {r['flagged_as_faster_aging']:.3f}"
              f" | s-hat {r['median_slope_estimate']:.3f}")


if __name__ == "__main__":
    main()
