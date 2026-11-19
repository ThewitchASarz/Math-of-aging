"""Experiment 033: how large a difference in benefit between sturdy and frail dogs breaks the slope test?

Model for one dog of size class c, entering at age x0, t years after enrolment (experiment 017's):
    hazard(t) = z * a * exp(b_c * x0) * r * exp(b_c * t)
  z   hidden personal risk factor: gamma with mean 1 and variance V among dogs alive at age 2, so among dogs
      alive at entry age x0 it is gamma with shape 1/V and scale V / (1 + V * H0), where
      H0 = (a / b_c) * (exp(b_c * x0) - exp(2 * b_c)).   V = 0.3 is an ASSUMPTION; V = 0.1 is about the size
      fitted to the UK life table in experiment 004.
  r   the treatment's risk multiplier for that dog (1 in the control group). It depends on p, the dog's rank
      (0 sturdiest, 1 frailest) in hidden risk among dogs of its entry age. No dog's rate of aging changes.
        halves:  r = 1 - c_sturdy if p < 0.5, else 1 - c_frail
        steady:  r = 1 - (c_top - (c_top - c_bottom) * p)      (the cut falls in a straight line with rank)
      Every cut size and both patterns are ASSUMPTIONS. The average cut is 50 percent in every case.

Designs: 1,300 dogs, half treated from enrolment.
  LONG   entry ages 7 to 12 (ASSUMPTION), 10 years, one size, no dropout; a, b = UK fit (dog_lifespan/results.json)
  SHORT  entry ages 10 to 14 (ASSUMPTION), 4 years, otherwise the same
  SIZES  LONG with five size classes, each with its own rate of aging (read through experiment 007), size mix =
         US population shares (ASSUMPTION), dropout 0.05 a year unrelated to health (ASSUMPTION)

Analyses: the plain slope test and the hidden-differences slope test, both imported unchanged from experiment 017
(tests_one_size, tests_sizes). Each is a 5% likelihood-ratio test of "risk cut only" against "risk cut and slope
change".

Run:  python3 benefit_size.py            (about 5 minutes on 2 cores)
      python3 benefit_size.py --quick
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
from scipy.stats import gamma as gamma_dist

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


h17 = _load("harder", "017_hidden_differences_harder", "harder.py")
fp, br = h17.fp, h17.br

SEED = 20261133
AGE0 = h17.AGE0
N = h17.N
DROP = h17.DROP                           # ASSUMPTION: dropout rate per year in the SIZES design
SHORT, LONG = h17.SHORT, h17.LONG         # entry ages: ASSUMPTION
# name: (design, sizes and dropout?, spread V, pattern, cut at the sturdy end, cut at the frail end)
CASES = {
    "10 years, spread 0.3 | halves 50 and 50": (LONG, False, 0.3, "halves", 0.50, 0.50),
    "10 years, spread 0.3 | halves 55 and 45": (LONG, False, 0.3, "halves", 0.55, 0.45),
    "10 years, spread 0.3 | halves 60 and 40": (LONG, False, 0.3, "halves", 0.60, 0.40),
    "10 years, spread 0.3 | halves 65 and 35": (LONG, False, 0.3, "halves", 0.65, 0.35),
    "10 years, spread 0.3 | halves 70 and 30": (LONG, False, 0.3, "halves", 0.70, 0.30),
    "10 years, spread 0.3 | steady 70 to 30": (LONG, False, 0.3, "steady", 0.70, 0.30),
    "10 years, spread 0.3 | steady 90 to 10": (LONG, False, 0.3, "steady", 0.90, 0.10),
    "10 years, spread 0.1 | halves 60 and 40": (LONG, False, 0.1, "halves", 0.60, 0.40),
    "10 years, spread 0.1 | halves 70 and 30": (LONG, False, 0.1, "halves", 0.70, 0.30),
    "4 years, spread 0.3 | halves 70 and 30": (SHORT, False, 0.3, "halves", 0.70, 0.30),
    "10 years, sizes, dropout, spread 0.3 | halves 70 and 30": (LONG, True, 0.3, "halves", 0.70, 0.30),
}


def risk_multiplier(p, pattern, c_sturdy, c_frail):
    """Risk multiplier of a treated dog whose rank in hidden risk is p (0 sturdiest, 1 frailest)."""
    if pattern == "halves":
        return 1.0 - np.where(p < 0.5, c_sturdy, c_frail)
    return 1.0 - (c_sturdy - (c_sturdy - c_frail) * p)


def simulate(design, mixed, v, pattern, c_sturdy, c_frail, a, b_one, b_class, shares, rng):
    lo, hi, follow = design
    x0 = rng.uniform(lo, hi, N)
    arm = (np.arange(N) % 2).astype(float)
    if mixed:
        cls = rng.choice(len(shares), size=N, p=shares)
        b = b_class[cls]
    else:
        cls = np.zeros(N, int)
        b = np.full(N, b_one)
    H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
    u = rng.gamma(1.0 / v, 1.0, N)                      # the same for every entry age, so its rank is the dog's
    z = u * v / (1.0 + v * H0)                          # rank among dogs of its own entry age
    p = gamma_dist.cdf(u, 1.0 / v)
    r = np.where(arm == 1, risk_multiplier(p, pattern, c_sturdy, c_frail), 1.0)
    h0 = z * a * np.exp(b * x0) * r
    t = np.log1p(-np.log(rng.random(N)) * b / h0) / b
    leave = rng.exponential(1.0 / DROP, N) if mixed else np.full(N, np.inf)
    end = np.minimum(leave, follow)
    died = (t <= end).astype(float)
    return (x0, arm, np.minimum(t, end), died), cls


def one_chunk(args):
    name, nsim, a_uk, b_uk, a_sz, b_class, shares, seed = args
    design, mixed, v, pattern, c_sturdy, c_frail = CASES[name]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    a = a_sz if mixed else a_uk
    b0 = float(np.average(b_class, weights=shares))
    out = {k: dict(flag=0, slower=0, s=[]) for k in ("plain", "new")}
    dead = 0.0
    for _ in range(nsim):
        data, cls5 = simulate(design, mixed, v, pattern, c_sturdy, c_frail, a, b_uk, b_class, shares, rng)
        dead += data[3].mean()
        res = h17.tests_sizes(data, br.STRATUM[cls5], a, b0) if mixed else h17.tests_one_size(data, a, b_uk)
        for label, (stat, s_hat, _) in res.items():
            out[label]["flag"] += int(stat > c1)
            out[label]["slower"] += int(stat > c1 and s_hat < 1.0)
            out[label]["s"].append(s_hat)
    return name, nsim, dead, out


def main():
    quick = "--quick" in sys.argv
    chunks, per = (1, 10) if quick else (2, 300)
    a_uk, b_uk = fp.dog_parameters()
    a_sz, b_class, shares, names = br.size_classes()
    # the average cut is 50 percent in every pattern
    grid = (np.arange(100000) + 0.5) / 100000
    for name, (_, _, _, pattern, cs, cf) in CASES.items():
        assert abs(risk_multiplier(grid, pattern, cs, cf).mean() - 0.5) < 1e-9, name
    jobs = [(name, per, a_uk, b_uk, a_sz, b_class, shares, SEED + 1000 * i + k)
            for i, name in enumerate(CASES) for k in range(chunks)]
    if quick:
        jobs = [j for j in jobs if j[0] in ("10 years, spread 0.3 | halves 60 and 40",
                                            "10 years, spread 0.3 | steady 90 to 10",
                                            "10 years, sizes, dropout, spread 0.3 | halves 70 and 30")]
    order = sorted(range(len(jobs)), key=lambda i: not CASES[jobs[i][0]][1])       # scheduling only: slow case first
    with Pool(2) as pool:
        parts = pool.map(one_chunk, [jobs[i] for i in order], chunksize=1)
    rows = []
    for name in CASES:
        mine = [p for p in parts if p[0] == name]
        if not mine:
            continue
        n = sum(p[1] for p in mine)
        design, mixed, v, pattern, cs, cf = CASES[name]
        for label in ("plain", "new"):
            s_all = [x for p in mine for x in p[3][label]["s"]]
            flag = sum(p[3][label]["flag"] for p in mine) / n
            slower = sum(p[3][label]["slower"] for p in mine) / n
            rows.append(dict(case=name, analysis=label, n_trials=n, years=design[2], entry_ages=list(design[:2]),
                             spread=v, sizes_and_dropout=mixed, pattern=pattern,
                             cut_sturdy_end=cs, cut_frail_end=cf,
                             share_of_dogs_dead=sum(p[2] for p in mine) / n,
                             slope_change_flagged=flag, flagged_as_slower_aging=slower,
                             flagged_as_faster_aging=flag - slower,
                             median_slope_estimate=float(np.median(s_all))))
    out = dict(note="Simulated trials of 1,300 dogs, half treated from enrolment. No dog's rate of aging changes in "
                    "any case; the treatment cuts risk by 50 percent on average, more in sturdy dogs than in frail "
                    "ones. 'flagged_as_slower_aging' is the share of trials in which a 5% likelihood-ratio test "
                    "rejects 'the slope did not change' with an estimated rate below 1. 'plain' ignores hidden "
                    "differences; 'new' adds one setting for the spread of a gamma-shaped hidden risk factor. "
                    "ASSUMPTIONS: the spread, entry ages, dropout 0.05 a year, size mix, and every cut and pattern.",
               seed=SEED, uk_a_per_year=a_uk, uk_b_per_year=b_uk, size_class_a_per_year=a_sz,
               size_class_b_per_year=[float(x) for x in b_class], size_class_shares=[float(x) for x in shares],
               dropout_per_year=DROP, trials_per_case=chunks * per, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case'][:56]:56s} | {r['analysis']:5s} | dead {r['share_of_dogs_dead']:.2f}"
              f" | slower {r['flagged_as_slower_aging']:.3f} | faster {r['flagged_as_faster_aging']:.3f}"
              f" | s-hat {r['median_slope_estimate']:.3f}")


if __name__ == "__main__":
    main()
