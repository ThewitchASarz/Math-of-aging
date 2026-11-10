"""Experiment 010: how large and how long must a dog trial be to recognise slowed aging?

Uses the simulation and slope test of experiment 001 (treatment continuous from enrolment, entry ages 7 to 12,
UK dog mortality, baseline re-estimated in each trial) over a grid of trial sizes and follow-up lengths.

Run:  python3 design_grid.py            (about 5 minutes on 2 cores)
      python3 design_grid.py --quick
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

SEED = 20261115
SIZES = [300, 600, 1000, 1300, 2000, 5000]
YEARS = [4.0, 6.0, 8.0, 10.0]
SLOWING = {"aging slowed 20%": 0.8, "aging slowed 30%": 0.7}
ENTRY = (7.0, 12.0)


def one_case(args):
    n, follow, effect, nsim, a, b, seed = args
    s = SLOWING[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    recognised = detected = 0
    deaths = []
    for _ in range(nsim):
        data = fp.simulate(n, ENTRY[0], ENTRY[1], follow, 1.0, s, a, b, rng)
        ll = fp.fit_all(data, a, b)
        recognised += 2 * (ll["full"][0] - ll["level"][0]) > c1
        detected += 2 * (ll["level"][0] - ll["null"][0]) > c1
        deaths.append(float(data[3].sum()))
    return dict(dogs=n, years=follow, effect=effect, n_trials=nsim, share_of_dogs_dead=float(np.mean(deaths)) / n,
                slowing_recognised=recognised / nsim, detected_by_usual_test=detected / nsim)


def main():
    quick = "--quick" in sys.argv
    nsim = 30 if quick else 300
    a, b = fp.dog_parameters()
    sizes = SIZES[:2] if quick else SIZES
    cases = [(n, y, e, nsim, a, b, SEED + 100000 * i + 1000 * j + k)
             for i, n in enumerate(sizes) for j, y in enumerate(YEARS) for k, e in enumerate(SLOWING)]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases)
    smallest = {}
    for e in SLOWING:
        for y in YEARS:
            ok = [r["dogs"] for r in rows if r["effect"] == e and r["years"] == y and r["slowing_recognised"] >= 0.8]
            smallest[f"{e}, {y:g} years"] = min(ok) if ok else None
    out = dict(note="Share of simulated trials in which the slope test of experiment 001 recognises a true slowing "
                    "of aging. Treatment is continuous; entry ages 7 to 12 are an assumption.",
               seed=SEED, trials_per_case=nsim, smallest_trial_reaching_80_percent=smallest, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for e in SLOWING:
        print(e)
        print("  dogs   " + "".join(f"{y:>8g}y" for y in YEARS))
        for n in sizes:
            vals = [next(r["slowing_recognised"] for r in rows if r["dogs"] == n and r["years"] == y and r["effect"] == e)
                    for y in YEARS]
            print(f"  {n:5d}  " + "".join(f"{v:9.2f}" for v in vals))
    print("smallest:", smallest)


if __name__ == "__main__":
    main()
