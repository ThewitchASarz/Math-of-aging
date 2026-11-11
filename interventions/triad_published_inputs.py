"""Power of a 580-dog, 3-year trial, using the inputs published in the trial's own design paper.

Source of inputs: Coleman et al. 2025, GeroScience 47:2851 (the TRIAD design paper),
https://pmc.ncbi.nlm.nih.gov/articles/PMC12181551/ :
  Gompertz level 0.0066 and slope 0.29 per year; half the dogs enter at age 7-8, a quarter at 8-9,
  a quarter at 9-10; 290 dogs per group; about one third of placebo dogs expected to die in 3 years.
ASSUMPTIONS (ours, not the paper's): a two-sided log-rank test at the 5% level, no early looks at
the data, no dropout, and the effect sizes listed in CASES.

This is a reconciliation, not an experiment: it explains why this project's low power figures and
the design paper's "greater than 80% power" do not contradict each other. Writes
triad_published_inputs.json.
"""
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy import integrate

A, B = 0.0066, 0.29          # published Gompertz level and slope (per year)
N_ARM = 290
Z = 1.959964                 # two-sided 5%

CASES = [                    # name, years of follow-up, effect
    ("no effect", 3, {}),
    ("20% risk cut for all 3 years", 3, {"r": 0.80}),
    ("30% risk cut for all 3 years", 3, {"r": 0.70}),
    ("38% risk cut for all 3 years", 3, {"r": 0.62}),
    ("22% risk cut, every dog followed until death", 60, {"r": 0.78}),
    ("20% risk cut during the dosing year only", 3, {"r_year1": 0.80}),
    ("20% slower rise of risk from enrolment", 3, {"s": 0.80}),
]


def entry_ages(n, rng):
    u = rng.random(n)
    return np.where(u < 0.5, 7, np.where(u < 0.75, 8, 9)) + rng.random(n)


def death_times(x0, rng, r=1.0, s=1.0, r_year1=None):
    """Years from enrolment to death. Hazard = A*exp(B*x0) * r * exp(B*s*t).
    r_year1: risk is multiplied by r_year1 in the first year only, then returns to normal."""
    e = rng.exponential(size=x0.size)
    h0 = A * np.exp(B * x0)
    if r_year1 is None:
        k = B * s
        return np.log1p(e * k / (h0 * r)) / k
    h_year1 = h0 * r_year1 * (np.exp(B) - 1) / B
    inside = np.log1p(e * B / (h0 * r_year1)) / B
    after = np.log(np.exp(B) + np.maximum(e - h_year1, 0) * B / h0) / B
    return np.where(e < h_year1, inside, after)


def logrank_z(t, died, treated):
    order = np.argsort(t)
    t, died, treated = t[order], died[order], treated[order]
    at_risk = t.size - np.arange(t.size)
    at_risk_treated = np.cumsum(treated[::-1])[::-1]
    share = at_risk_treated / at_risk
    observed, expected = (died * treated).sum(), (died * share).sum()
    return (observed - expected) / np.sqrt((died * share * (1 - share)).sum())


def run_case(follow, effect, n_sim, rng):
    hits, placebo_dead = 0, 0.0
    treated = np.repeat([0.0, 1.0], N_ARM)
    for _ in range(n_sim):
        x0 = entry_ages(2 * N_ARM, rng)
        t = np.concatenate([death_times(x0[:N_ARM], rng), death_times(x0[N_ARM:], rng, **effect)])
        died = (t < follow).astype(float)
        placebo_dead += died[:N_ARM].mean()
        hits += abs(logrank_z(np.minimum(t, follow), died, treated)) > Z
    return hits / n_sim, placebo_dead / n_sim


def mean_total_lifespan(r):
    """Average age at death for the published mix of entry ages, with risk multiplied by r for life."""
    out = []
    for x, weight in ((7.5, 0.5), (8.5, 0.25), (9.5, 0.25)):
        surv = lambda t: np.exp(-A * np.exp(B * x) * r * (np.exp(B * t) - 1) / B)
        out.append(weight * (x + integrate.quad(surv, 0, 40)[0]))
    return sum(out)


def main(n_sim=4000):
    rng = np.random.default_rng(20261111)
    base = mean_total_lifespan(1.0)
    rows = []
    for name, follow, effect in CASES:
        power, dead = run_case(follow, effect, n_sim, rng)
        row = {"case": name, "follow_up_years": follow, "power": round(power, 4),
               "placebo_share_dead_in_window": round(dead, 4)}
        if "r" in effect:
            row["mean_total_lifespan_gain_percent"] = round(100 * (mean_total_lifespan(effect["r"]) / base - 1), 2)
        rows.append(row)
        print(row)
    out = {"source": "Coleman et al. 2025, GeroScience 47:2851 (PMC12181551)", "gompertz_level": A,
           "gompertz_slope_per_year": B, "dogs_per_group": N_ARM, "simulated_trials_per_case": n_sim,
           "test": "two-sided log-rank, 5% level (ASSUMPTION)", "cases": rows}
    name = "triad_published_inputs.json" if n_sim == 4000 else "triad_published_inputs_quick.json"
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), name), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main(300 if "--quick" in sys.argv else 4000)
