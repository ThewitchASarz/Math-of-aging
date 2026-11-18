"""Experiment 009: what the two-period test reports when a benefit fades slowly or only part of the way.

The model: the rapamycin trial (TRIAD) as published: Gompertz level 0.0066 and slope 0.29 per year, entry
ages 7-8 (half), 8-9 (quarter), 9-10 (quarter), 12 months of dosing. Source: Coleman et al. 2025,
GeroScience 47:2851 (PMC12181551), as recorded in interventions/triad_published_inputs.py and
papers/literature-check.md. Not re-opened in this run.

Model for one dog (t = years since enrolment, x0 = entry age):
    hazard(t) = m(t) * a * exp(b * (x0 + t))
    m(t) = 1 for placebo dogs. For treated dogs:
        m(t) = m1                                    for t < 1   (the dosing year)
        m(t) = 1 - d * 2 ** (-(t - 1) / half_life)   for t >= 1  (half_life = None: m(t) = 1 - d for good)
ASSUMPTIONS: every effect size (m1, d), every half-life, and designs D2 (10 years) and D3 (5,000 dogs).
Design D1 (580 dogs, 3 years) is the published one.

Analyses of each simulated trial (both reused):
    two-period test     : experiment 006 (is the risk ratio after dosing different from 1?)
    any-difference test : two-sided log-rank test at the 5% level (interventions/triad_published_inputs.py)

Run:  python3 partial_fade.py            (about 2 minutes on 2 cores)
      python3 partial_fade.py --quick
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
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


two = _load("two_period", "experiments", "006_stop_and_watch", "two_period.py")
triad = _load("triad_published_inputs", "interventions", "triad_published_inputs.py")

SEED = 20261118
A, B = triad.A, triad.B                  # 0.0066 and 0.29 per year: published
TAU = 1.0                                # published: 12 months of dosing
Z = 1.959964
assert two.TAU == TAU

# name: (dogs, years followed, simulated trials). D1 is published; D2 and D3 are ASSUMPTIONS.
DESIGNS = {
    "D1: 580 dogs, 3 years (published design)": (580, 3.0, 1000),
    "D2: 580 dogs, 10 years": (580, 10.0, 1000),
    "D3: 5,000 dogs, 3 years": (5000, 3.0, 400),
}
# name: (m1, d, half_life in years or None). ALL ASSUMPTIONS.
EFFECTS = {
    "no effect": (1.0, 0.0, None),
    "fades at once": (0.5, 0.0, None),
    "half lasts": (0.5, 0.25, None),
    "a fifth lasts": (0.5, 0.10, None),
    "half of a smaller benefit lasts": (0.7, 0.15, None),
    "quick fade (halves every 6 months)": (0.5, 0.5, 0.5),
    "slow fade (halves every 2 years)": (0.5, 0.5, 2.0),
}


def multiplier(t, m1, d, half_life):
    """Treated dogs' risk multiplier at time t since enrolment."""
    t = np.asarray(t, float)
    lam = 0.0 if half_life is None else np.log(2.0) / half_life
    return np.where(t < TAU, m1, 1.0 - d * np.exp(-lam * np.maximum(t - TAU, 0.0)))


def cumulative(t, m1, d, half_life, b):
    """G(t) = integral from 0 to t of m(s) * exp(b * s) ds, in closed form."""
    t = np.asarray(t, float)
    lam = 0.0 if half_life is None else np.log(2.0) / half_life
    k = b - lam                                        # never zero for the half-lives used here
    first = m1 * np.expm1(b * np.minimum(t, TAU)) / b
    tt = np.maximum(t, TAU)
    later = (np.exp(b * tt) - np.exp(b * TAU)) / b - d * np.exp(lam * TAU) * (np.exp(k * tt) - np.exp(k * TAU)) / k
    return first + later


def simulate(n, follow, effect, a, b, rng):
    """Death times by inverting the cumulative risk. Arms alternate. Returns (x0, arm, time, died)."""
    m1, d, half_life = effect
    x0 = triad.entry_ages(n, rng)
    arm = (np.arange(n) % 2).astype(float)
    target = -np.log(rng.random(n)) / (a * np.exp(b * x0))        # value of G at death
    grid = np.linspace(0.0, follow + 0.01, int(round((follow + 0.01) / 0.0005)) + 1)
    t_placebo = np.log1p(target * b) / b
    t_treated = np.interp(target, cumulative(grid, m1, d, half_life, b), grid)   # past the grid: censored anyway
    t = np.where(arm == 1, t_treated, t_placebo)
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def one_case(args):
    design, effect, seed, nsim = args
    n, follow, _ = DESIGNS[design]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    fires = after = lasting = worse = 0
    on_hat, off_hat, dead_p, dead_t, deaths_after = [], [], [], [], []
    for _ in range(nsim):
        data = simulate(n, follow, EFFECTS[effect], A, B, rng)
        _, arm, t, died = data
        dead_p.append(died[arm == 0].mean())
        dead_t.append(died[arm == 1].mean())
        deaths_after.append(float((died * (t > TAU)).sum()))
        fires += abs(triad.logrank_z(t, died, arm)) > Z
        on = two.fit("on_only", data, A, B)
        both = two.fit("both", data, A, B, extra=[on[1][2], 0.0])
        flag = 2 * (both[0] - on[0]) > c1
        r_off = float(np.exp(both[1][3]))
        after += flag
        lasting += flag and r_off < 1.0
        worse += flag and r_off > 1.0
        on_hat.append(float(np.exp(both[1][2])))
        off_hat.append(r_off)
    m1, d, half_life = EFFECTS[effect]
    return dict(design=design, effect=effect, dosing_year_multiplier=m1, cut_just_after_dosing=d,
                half_life_years=half_life, true_multiplier_at_end_of_follow_up=float(multiplier(follow, m1, d, half_life)),
                n_dogs=n, follow_up_years=follow, n_trials=nsim,
                placebo_share_dead=float(np.mean(dead_p)), treated_share_dead=float(np.mean(dead_t)),
                mean_deaths_after_dosing=float(np.mean(deaths_after)),
                any_difference_detected=fires / nsim, after_dosing_difference_reported=after / nsim,
                reported_lasting_benefit=lasting / nsim, reported_worse_after_stopping=worse / nsim,
                median_ratio_in_dosing_year=float(np.median(on_hat)), median_ratio_after_dosing=float(np.median(off_hat)))


def self_check():
    """The closed-form cumulative risk against a numerical integral, and the simulator against exact survival."""
    worst = 0.0
    s = np.linspace(0.0, 10.0, 200001)
    for eff in EFFECTS.values():
        y = multiplier(s, *eff) * np.exp(B * s)
        num = np.concatenate([[0.0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(s))])
        worst = max(worst, float(np.max(np.abs(num - cumulative(s, *eff, B)) / (1e-9 + cumulative(s, *eff, B)))))
    rng = np.random.default_rng(1)
    eff = EFFECTS["slow fade (halves every 2 years)"]
    x0, arm, t, died = simulate(400_000, 10.0, eff, A, B, rng)
    gap = 0.0
    for lo in (7, 8, 9):
        sel = (arm == 1) & (x0 >= lo) & (x0 < lo + 1)
        xs = lo + (np.arange(200) + 0.5) / 200
        for tt in (1.0, 3.0, 6.0):
            exact = np.mean(np.exp(-A * np.exp(B * xs) * cumulative(tt, *eff, B)))
            gap = max(gap, abs(float((t[sel] > tt).mean()) - exact))
    return worst, gap


def main():
    quick = "--quick" in sys.argv
    worst, gap = self_check()
    assert worst < 1e-4 and gap < 0.006, (worst, gap)
    jobs = [(d, e, SEED + 1000 * i + j, 30 if quick else DESIGNS[d][2])
            for i, d in enumerate(DESIGNS) for j, e in enumerate(EFFECTS)]
    jobs.sort(key=lambda job: -DESIGNS[job[0]][0] * job[3])       # longest first, so the two workers balance
    with Pool(2) as pool:
        rows = pool.map(one_case, jobs, chunksize=1)
    order = {(d, e): k for k, (d, e) in enumerate((d, e) for d in DESIGNS for e in EFFECTS)}
    rows.sort(key=lambda r: order[(r["design"], r["effect"])])
    out = dict(
        note="Simulated trials with the rapamycin trial's published death rates, entry ages and one year of dosing. "
             "'reported_lasting_benefit' is experiment 006's two-period test firing with an after-dosing ratio "
             "below 1; 'any_difference_detected' is a two-sided 5% log-rank test. All effect sizes, fading speeds "
             "and designs D2 and D3 are assumptions.",
        triad_level_per_year=A, triad_slope_per_year=B, dosing_years=TAU, seed=SEED,
        self_check=dict(cumulative_risk_max_relative_error=worst, survival_max_gap_400k_dogs=gap), cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['design'][:2]} {r['effect'][:34]:34s} | dead P {r['placebo_share_dead']:.3f} T "
              f"{r['treated_share_dead']:.3f} | any {r['any_difference_detected']:.3f} | lasting "
              f"{r['reported_lasting_benefit']:.3f} worse {r['reported_worse_after_stopping']:.3f} | on "
              f"{r['median_ratio_in_dosing_year']:.2f} off {r['median_ratio_after_dosing']:.3f} | true at end "
              f"{r['true_multiplier_at_end_of_follow_up']:.3f}")


if __name__ == "__main__":
    main()
