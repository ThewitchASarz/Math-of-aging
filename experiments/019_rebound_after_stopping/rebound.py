"""Experiment 019: what a one-year-dosing, three-year trial reports when risk rebounds after dosing stops,
when part of the benefit lasts, and when the dogs that survive the dosing year are a different mix.

The trial: the rapamycin trial (TRIAD) as published: 580 dogs, Gompertz level 0.0066 and slope 0.29 per year,
entry ages 7-8 (half), 8-9 (quarter), 9-10 (quarter), 12 months of dosing, 36 months of follow-up.
Source: Coleman et al. 2025, GeroScience 47:2851 (PMC12181551), as recorded in
interventions/triad_published_inputs.py and papers/literature-check.md. Not re-opened in this run.

Model for one dog (t = years since enrolment, x0 = entry age, z = hidden personal risk factor):
    hazard(t) = z * m(t) * a * exp(b * (x0 + t))
    m(t) = 1 for placebo dogs; for treated dogs m1 for t < 1, m2 for 1 <= t < 1.5, m3 for t >= 1.5.
    z = 1 for every dog, or gamma-distributed with mean 1 and variance v ("hidden differences").
ASSUMPTIONS: every effect size (m1, m2, m3), the six-month length of the short rebound, and v = 0.30.
v = 0.10 is the tentative value fitted to the UK table in experiment 004.

Analyses of each simulated trial (all reused):
    any-difference test : two-sided log-rank test at the 5% level (interventions/triad_published_inputs.py)
    two-period test     : experiment 006 (is the risk ratio after dosing different from 1?)
    slope test          : experiment 001 (reported, not scored; it must not be used on this design)

Run:  python3 rebound.py            (about 1 minute on 2 cores)
      python3 rebound.py --quick
"""
import importlib.util
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from multiprocessing import Pool

import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "experiments", "001_trial_fingerprint", "fingerprint.py")
two = _load("two_period", "experiments", "006_stop_and_watch", "two_period.py")
triad = _load("triad_published_inputs", "interventions", "triad_published_inputs.py")

SEED = 20261119
A, B = triad.A, triad.B                  # 0.0066 and 0.29 per year: published
N_DOGS, FOLLOW, TAU = 580, 3.0, 1.0      # published design: 580 dogs, 36 months, 12 months of dosing
REBOUND_END = 1.5                        # ASSUMPTION: a short rebound lasts six months
Z = 1.959964
assert two.TAU == TAU

# name: (m1, m2, m3) = risk multiplier in the dosing year, months 12-18, months 18-36. ALL ASSUMPTIONS.
EFFECTS = {
    "no effect": (1.0, 1.0, 1.0),
    "fades": (0.5, 1.0, 1.0),
    "lasting rebound": (0.5, 1.5, 1.5),
    "short rebound": (0.5, 2.0, 1.0),
    "small lasting rebound": (0.8, 1.25, 1.25),
    "rebound with no benefit first": (1.0, 1.5, 1.5),
    "half the benefit lasts": (0.5, 0.75, 0.75),
    "half of a smaller benefit lasts": (0.7, 0.85, 0.85),
}
# (effect, variance of the hidden risk factor). 0.10: experiment 004 (tentative). 0.30: ASSUMPTION.
CASES = [(e, 0.0) for e in EFFECTS] + [("no effect", 0.10), ("fades", 0.10), ("no effect", 0.30), ("fades", 0.30)]


def simulate(x0, mult, v, a, b, rng, follow=FOLLOW):
    """Death times by inverse transform through three segments with multipliers mult = (m1, m2, m3).
    Arms alternate. v is the variance of a gamma hidden risk factor with mean 1 (0 = none)."""
    n = x0.size
    arm = (np.arange(n) % 2).astype(float)
    z = rng.gamma(1.0 / v, v, n) if v > 0 else np.ones(n)
    target = -np.log(rng.random(n)) / z               # baseline-scale cumulative hazard at death
    h0 = a * np.exp(b * x0)
    cuts = [0.0, TAU, REBOUND_END, np.inf]
    t = np.full(n, np.inf)
    done = np.zeros(n, bool)
    used = np.zeros(n)
    for k in range(3):
        m = np.where(arm == 1, mult[k], 1.0)
        lo, hi = cuts[k], cuts[k + 1]
        seg = m * h0 * (np.exp(b * hi) - np.exp(b * lo)) / b if np.isfinite(hi) else np.full(n, np.inf)
        here = ~done & (target - used < seg)
        t[here] = np.log(np.exp(b * lo) + (target[here] - used[here]) * b / (m[here] * h0[here])) / b
        done |= here
        used += np.where(np.isfinite(seg), seg, 0.0)
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def one_case(args):
    (effect, v), nsim, seed = args
    mult = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    fires = treated_more = after = after_above = after_below = slope = slope_faster = 0
    on_hat, off_hat, dead_p, dead_t, dead_p1 = [], [], [], [], []
    for _ in range(nsim):
        data = simulate(triad.entry_ages(N_DOGS, rng), mult, v, A, B, rng)
        _, arm, t, died = data
        dead_p.append(died[arm == 0].mean())
        dead_t.append(died[arm == 1].mean())
        dead_p1.append((died[arm == 0] * (t[arm == 0] < TAU)).mean())
        zlr = triad.logrank_z(t, died, arm)             # positive: more treated deaths than expected
        fires += abs(zlr) > Z
        treated_more += zlr > Z
        on = two.fit("on_only", data, A, B)
        both = two.fit("both", data, A, B, extra=[on[1][2], 0.0])
        flag = 2 * (both[0] - on[0]) > c1
        r_off = float(np.exp(both[1][3]))
        after += flag
        after_above += flag and r_off > 1.0
        after_below += flag and r_off < 1.0
        on_hat.append(float(np.exp(both[1][2])))
        off_hat.append(r_off)
        ll = fp.fit_all(data, A, B)
        sflag = 2 * (ll["full"][0] - ll["level"][0]) > c1
        slope += sflag
        slope_faster += sflag and np.exp(ll["full"][1][3]) > 1.0
    return dict(effect=effect, multipliers=list(mult), hidden_variance=v, n_trials=nsim,
                placebo_share_dead=float(np.mean(dead_p)), treated_share_dead=float(np.mean(dead_t)),
                placebo_share_dead_in_dosing_year=float(np.mean(dead_p1)),
                any_difference_detected=fires / nsim, any_difference_with_more_treated_deaths=treated_more / nsim,
                after_dosing_difference_reported=after / nsim, reported_worse_after_stopping=after_above / nsim,
                reported_lasting_benefit=after_below / nsim,
                median_ratio_in_dosing_year=float(np.median(on_hat)), median_ratio_after_dosing=float(np.median(off_hat)),
                slope_change_flagged=slope / nsim, flagged_as_faster_aging=slope_faster / nsim)


def survivor_ratio(m1, v, a, b, ages, weights, t_lo=TAU, t_hi=FOLLOW):
    """Direct calculation for a treatment that multiplies risk by m1 in the dosing year and does nothing after.
    With a gamma hidden factor of variance v, the death rate of the dogs still alive in a group is
    h(t) / (1 + v * H(t)), where H is that group's cumulative risk. Returns the treated-to-placebo ratio of
    death rates among survivors: just after dosing stops, at the end, and averaged over the after-dosing
    period weighting each moment by the placebo deaths expected then. These three compare dogs of the same
    entry age. 'pooled_deaths_per_dog_year' compares the whole groups."""
    def pieces(x0, t):
        h0 = a * np.exp(b * x0)
        Hc = h0 * np.expm1(b * t) / b
        Ht = Hc - (1.0 - m1) * h0 * np.expm1(b * TAU) / b
        ratio = (1.0 + v * Hc) / (1.0 + v * Ht)
        dens = h0 * np.exp(b * t) * (1.0 + v * Hc) ** (-1.0 / v - 1.0)      # placebo death density
        return ratio, dens
    num = sum(w * quad(lambda t: np.prod(pieces(x, t)), t_lo, t_hi)[0] for x, w in zip(ages, weights))
    den = sum(w * quad(lambda t: pieces(x, t)[1], t_lo, t_hi)[0] for x, w in zip(ages, weights))
    at = lambda t: float(sum(w * pieces(x, t)[0] for x, w in zip(ages, weights)))

    # Pooled version: deaths per dog-year after dosing in each whole group, entry ages combined. This also
    # counts the shift in the age mix of the survivors, and it is what big_sample_ratio measures.
    def rate(m):
        deaths = years = 0.0
        for x, w in zip(ages, weights):
            h0 = a * np.exp(b * x)
            surv = lambda t: (1.0 + v * (h0 * np.expm1(b * t) / b - (1.0 - m) * h0 * np.expm1(b * TAU) / b)) ** (-1.0 / v)
            deaths += w * (surv(t_lo) - surv(t_hi))
            years += w * quad(surv, t_lo, t_hi)[0]
        return deaths / years
    return dict(just_after_dosing=at(t_lo), at_end=at(t_hi), average_after_dosing=float(num / den),
                pooled_deaths_per_dog_year=float(rate(m1) / rate(1.0)))


def big_sample_ratio(v, rng, n=2_000_000):
    """Check of survivor_ratio by brute force: deaths per dog-year after dosing, treated over placebo."""
    _, arm, t, died = simulate(triad.entry_ages(n, rng), EFFECTS["fades"], v, A, B, rng)
    late = t > TAU
    rate = lambda g: (died[late & (arm == g)]).sum() / (t[late & (arm == g)] - TAU).sum()
    return float(rate(1) / rate(0))


def main():
    quick = "--quick" in sys.argv
    nsim = 30 if quick else 1000
    jobs = [(c, nsim, SEED + 100 * i) for i, c in enumerate(CASES)]
    with Pool(2) as pool:
        rows = pool.map(one_case, jobs, chunksize=1)
    uk_a, uk_b = fp.dog_parameters()
    grid = 7.0 + (np.arange(300) + 0.5) / 100.0             # published entry-age mix, as weights on cell midpoints
    w = np.where(grid < 8.0, 0.5, 0.25) / 100.0
    old = 10.0 + (np.arange(400) + 0.5) / 100.0
    side = {"rapamycin design, fades (0.5 then 1)": {
                f"hidden variance {v:g}": survivor_ratio(0.5, v, A, B, grid, w) for v in (0.10, 0.30)},
            "older dogs (entry 10 to 14, UK rates; ASSUMPTION), fades (0.5 then 1)": {
                f"hidden variance {v:g}": survivor_ratio(0.5, v, uk_a, uk_b, old, np.full(old.size, 1 / old.size))
                for v in (0.10, 0.30)}}
    if not quick:
        rng = np.random.default_rng(SEED + 7)
        side["brute-force check, 2 million dogs, rapamycin design, hidden variance 0.3"] = big_sample_ratio(0.30, rng)
    out = dict(
        note="Simulated trials with the rapamycin trial's published design. 'any_difference_detected' is a "
             "two-sided 5% log-rank test; 'after_dosing_difference_reported' is experiment 006's two-period test; "
             "'slope_change_flagged' is experiment 001's test. All effect sizes are assumptions.",
        triad_level_per_year=A, triad_slope_per_year=B, n_dogs=N_DOGS, follow_up_years=FOLLOW, dosing_years=TAU,
        short_rebound_ends_at_years=REBOUND_END, seed=SEED, trials_per_case=nsim,
        survivor_risk_ratio_after_dosing=side, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['effect'][:31]:31s} v={r['hidden_variance']:.2f} | dead P {r['placebo_share_dead']:.3f} "
              f"T {r['treated_share_dead']:.3f} | any {r['any_difference_detected']:.3f} "
              f"(treated more {r['any_difference_with_more_treated_deaths']:.3f}) | after "
              f"{r['after_dosing_difference_reported']:.3f} (worse {r['reported_worse_after_stopping']:.3f}, lasting "
              f"{r['reported_lasting_benefit']:.3f}) | on {r['median_ratio_in_dosing_year']:.2f} "
              f"off {r['median_ratio_after_dosing']:.2f} | slope {r['slope_change_flagged']:.3f} "
              f"(faster {r['flagged_as_faster_aging']:.3f})")
    print(json.dumps(side, indent=1))


if __name__ == "__main__":
    main()
