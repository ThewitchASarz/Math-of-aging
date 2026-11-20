"""Experiment 038: does a larger trial make the "benefit differs between dogs" trap worse?

Model for one dog entering at age x0, t years after enrolment (that of experiments 017, 033 and 037):
    hazard(t) = z * a * exp(b * x0) * r * exp(s * b * t)
  z   hidden personal risk factor: gamma with mean 1 and variance V = 0.3 (ASSUMPTION) among dogs alive at age 2,
      so among dogs alive at entry age x0 it is gamma with shape 1/V and scale V / (1 + V * H0), where
      H0 = (a / b) * (exp(b * x0) - exp(2 * b)).
  r   the treatment's risk multiplier for that dog (1 in the control group): 1 - c_sturdy for the half of dogs
      with lower hidden risk among dogs of the same entry age, 1 - c_frail for the other half.
      Every cut size is an ASSUMPTION.
  s   the treatment's multiplier on the rate of aging: 1 in every case except "true slowing" (0.8, with r = 1).
  a, b  the UK fit in dog_lifespan/results.json.

Designs: n dogs (600, 1,300, 2,600 or 5,000), half treated from enrolment, one size, no dropout (ASSUMPTION).
  LONG   entry ages 7 to 12 (ASSUMPTION), 10 years (ASSUMPTION)
  SHORT  entry ages 10 to 14 (ASSUMPTION), 4 years (ASSUMPTION)

Analyses: the plain slope test and the hidden-differences slope test, imported unchanged from experiment 017
(tests_one_size). Each is a 5% likelihood-ratio test of "risk cut only" against "risk cut and slope change".

Run:  python3 larger_trial.py               (everything; about 5 minutes on 2 cores)
      python3 larger_trial.py --quick
The same work in four shorter runs (each under 2 minutes), giving the same results.json:
      python3 larger_trial.py --part 1      (and 2, 3, 4; each writes parts/part_N.json)
      python3 larger_trial.py --merge
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
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


h17 = _load("harder", "017_hidden_differences_harder", "harder.py")
fp = h17.fp

SEED = 20261138
AGE0 = h17.AGE0
V = 0.3                                   # ASSUMPTION: spread of hidden risk
SHORT, LONG = h17.SHORT, h17.LONG         # entry ages: ASSUMPTION
DESIGNS = {"10 years": LONG, "4 years": SHORT}
TRIALS = {600: 600, 1300: 600, 2600: 300, 5000: 300}
# pattern name: (cut in the sturdier half, cut in the frailer half, rate multiplier)   -- all ASSUMPTIONS
PATTERNS = {
    "55 and 45": (0.55, 0.45, 1.0),
    "60 and 40": (0.60, 0.40, 1.0),
    "40 and 60": (0.40, 0.60, 1.0),
    "42 and 18": (0.42, 0.18, 1.0),
    "50 and 50": (0.50, 0.50, 1.0),
    "true slowing 20%": (0.0, 0.0, 0.8),
}
# (design, pattern, sizes)
PLAN = [
    ("10 years", "55 and 45", (600, 2600, 5000)),
    ("10 years", "60 and 40", (600, 2600, 5000)),
    ("10 years", "40 and 60", (600, 2600, 5000)),
    ("10 years", "42 and 18", (600, 2600, 5000)),
    ("10 years", "true slowing 20%", (600, 2600, 5000)),
    ("10 years", "50 and 50", (5000,)),
    ("4 years", "55 and 45", (600, 1300, 2600, 5000)),
    ("4 years", "60 and 40", (600, 1300, 2600, 5000)),
    ("4 years", "true slowing 20%", (600, 2600, 5000)),
]
CASES = [(d, p, n) for d, p, sizes in PLAN for n in sizes]          # 27 cases, fixed order
CHUNKS = 2
# The 1,300-dog figures already recorded (share of trials reporting the direction named), used only for the
# square-root rule in --merge. Sources: experiments 012, 017, 033 and 037 results.json, hidden-differences test.
KNOWN_1300 = {
    ("10 years", "55 and 45"): ("slower", 0.098),          # 033
    ("10 years", "60 and 40"): ("slower", 0.287),          # 033
    ("10 years", "40 and 60"): ("faster", 0.203),          # 037
    ("10 years", "42 and 18"): ("slower", 0.285),          # 037
    ("10 years", "true slowing 20%"): ("slower", 0.814),   # 012, world C
    ("4 years", "true slowing 20%"): ("slower", 0.258),    # 017
}


def simulate(n, design, c_sturdy, c_frail, s, a, b, rng):
    lo, hi, follow = design
    x0 = rng.uniform(lo, hi, n)
    arm = (np.arange(n) % 2).astype(float)
    H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
    u = rng.gamma(1.0 / V, 1.0, n)                      # the same for every entry age, so its rank is the dog's
    z = u * V / (1.0 + V * H0)                          # rank among dogs of its own entry age
    frail = u > gamma_dist.ppf(0.5, 1.0 / V)
    r = np.where(arm == 1, 1.0 - np.where(frail, c_frail, c_sturdy), 1.0)
    h0 = z * a * np.exp(b * x0) * r
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def one_chunk(args):
    idx, nsim, a, b, seed = args
    dname, pname, n = CASES[idx]
    c_sturdy, c_frail, s = PATTERNS[pname]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    out = {k: dict(flag=0, slower=0, s=[]) for k in ("plain", "new")}
    dead = 0.0
    for _ in range(nsim):
        data = simulate(n, DESIGNS[dname], c_sturdy, c_frail, s, a, b, rng)
        dead += data[3].mean()
        for label, (stat, s_hat, _) in h17.tests_one_size(data, a, b).items():
            out[label]["flag"] += int(stat > c1)
            out[label]["slower"] += int(stat > c1 and s_hat < 1.0)
            out[label]["s"].append(s_hat)
    return idx, nsim, dead, out


def rows_for(parts):
    rows = []
    for idx, (dname, pname, n) in enumerate(CASES):
        mine = [p for p in parts if p[0] == idx]
        if not mine:
            continue
        nt = sum(p[1] for p in mine)
        cs, cf, s = PATTERNS[pname]
        for label in ("plain", "new"):
            s_all = [x for p in mine for x in p[3][label]["s"]]
            flag = sum(p[3][label]["flag"] for p in mine) / nt
            slower = sum(p[3][label]["slower"] for p in mine) / nt
            rows.append(dict(case=f"{dname} | {pname} | {n} dogs", design=dname, pattern=pname, n_dogs=n,
                             analysis=label, n_trials=nt, cut_sturdier_half=cs, cut_frailer_half=cf,
                             rate_multiplier=s, share_of_dogs_dead=sum(p[2] for p in mine) / nt,
                             slope_change_flagged=flag, flagged_as_slower_aging=slower,
                             flagged_as_faster_aging=flag - slower,
                             median_slope_estimate=float(np.median(s_all)),
                             sd_of_log_slope_estimate=float(np.std(np.log(s_all)))))
    return rows


def root_n_rule(p_1300, n):
    """Share of trials expected at n dogs if the error in the estimate is fixed and the noise shrinks as 1/sqrt(n)."""
    delta = norm.ppf(p_1300) + norm.ppf(0.975)
    return float(norm.cdf(delta * np.sqrt(n / 1300.0) - norm.ppf(0.975)))


def _make_parts(budget):
    """Split the 27 cases into parts of roughly equal work (dogs x trials)."""
    part, load, out = 1, 0.0, {}
    for idx, (_, _, n) in enumerate(CASES):
        w = n * TRIALS[n]
        if load + w > budget and load > 0:
            part, load = part + 1, 0.0
        out[idx] = part
        load += w
    return out


PART_OF = _make_parts(8.2e6)
N_PARTS = max(PART_OF.values())


def write_results(rows, a, b):
    order = {f"{d} | {p} | {n} dogs": i for i, (d, p, n) in enumerate(CASES)}
    rows.sort(key=lambda r: (order[r["case"]], r["analysis"] == "new"))
    assert len(rows) == 2 * len(CASES)
    rule = []
    for (dname, pname), (direction, p0) in KNOWN_1300.items():
        for r in rows:
            if r["design"] == dname and r["pattern"] == pname and r["analysis"] == "new":
                rule.append(dict(design=dname, pattern=pname, n_dogs=r["n_dogs"], direction=direction,
                                 known_at_1300_dogs=p0, square_root_rule=root_n_rule(p0, r["n_dogs"]),
                                 simulated=r[f"flagged_as_{direction}_aging"]))
    out = dict(note="Simulated trials of 600 to 5,000 dogs, half treated from enrolment; 10 years with entry "
                    "ages 7 to 12, or 4 years with entry ages 10 to 14. In every case but 'true slowing 20%' no "
                    "dog's rate of aging changes: the treatment cuts risk by one fraction in the sturdier half "
                    "of dogs and another in the frailer half. 'flagged_as_slower_aging' ('..faster..') is the "
                    "share of trials in which a 5% likelihood-ratio test rejects 'the slope did not change' "
                    "with an estimated rate below (above) 1. 'plain' ignores hidden differences; 'new' adds one "
                    "setting for the spread of a gamma-shaped hidden risk factor. 'square_root_rule' applies "
                    "'fixed error, noise shrinking as 1/sqrt(dogs)' to the figure recorded earlier for 1,300 "
                    "dogs. ASSUMPTIONS: the spread of hidden risk (0.3), entry ages, trial lengths, no dropout, "
                    "and every cut.",
               seed=SEED, uk_a_per_year=a, uk_b_per_year=b, spread=V,
               trials_per_case={str(k): v for k, v in TRIALS.items()}, cases=rows, square_root_rule=rule)
    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case']:42s} | {r['analysis']:5s} | dead {r['share_of_dogs_dead']:.3f}"
              f" | slower {r['flagged_as_slower_aging']:.3f} | faster {r['flagged_as_faster_aging']:.3f}"
              f" | s-hat {r['median_slope_estimate']:.3f}")


def jobs_for(indices, a, b):
    jobs = [(idx, TRIALS[CASES[idx][2]] // CHUNKS, a, b, SEED + 1000 * idx + c) for idx in indices
            for c in range(CHUNKS)]
    return sorted(jobs, key=lambda j: -CASES[j[0]][2] * j[1])     # largest jobs first, for an even split


def main():
    a, b = fp.dog_parameters()
    if "--merge" in sys.argv:
        rows = []
        for k in range(1, N_PARTS + 1):
            with open(os.path.join(HERE, "parts", f"part_{k}.json")) as fh:
                rows += json.load(fh)
        write_results(rows, a, b)
        return
    if "--quick" in sys.argv:
        jobs, name = [(idx, 6, a, b, SEED + 1000 * idx) for idx in (2, 17, 26)], "results_quick.json"
    elif "--part" in sys.argv:
        k = int(sys.argv[sys.argv.index("--part") + 1])
        jobs = jobs_for([i for i in range(len(CASES)) if PART_OF[i] == k], a, b)
        os.makedirs(os.path.join(HERE, "parts"), exist_ok=True)
        name = os.path.join("parts", f"part_{k}.json")
    else:
        jobs, name = jobs_for(range(len(CASES)), a, b), None
    with Pool(2) as pool:
        parts = pool.map(one_chunk, jobs, chunksize=1)
    rows = rows_for(parts)
    if name is None:
        write_results(rows, a, b)
        return
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(rows, fh, indent=1)
    for r in rows:
        print(f"{r['case']:42s} | {r['analysis']:5s} | slower {r['flagged_as_slower_aging']:.3f}"
              f" | faster {r['flagged_as_faster_aging']:.3f} | s-hat {r['median_slope_estimate']:.3f}")


if __name__ == "__main__":
    main()
