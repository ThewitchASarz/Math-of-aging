"""Experiment 015: the 10-year trial of experiment 011 with a delayed effect and owners who do not treat.

Trial: the "realistic" scenario of experiment 011, unchanged (size mix, hidden risk factor, 5 percent of dogs
lost to follow-up a year, entry ages 7 to 12, 10 years, half the dogs assigned to treatment).

Model for one dog of size class c, entering at age x0, t years after enrolment:
    untreated:  hazard(t) = z * a * exp(b_c * x0) * exp(b_c * t)
    treated:    hazard(t) = z * a * exp(b_c * x0) * exp(b_c * g(t) + L(t))
  g(t)  "aging clock". Slowing d that builds up in a straight line over D years:
            g(t) = t - d * t^2 / (2 D)        for t <= D
            g(t) = (1 - d) * t + d * D / 2    for t >  D        (D = 0: g(t) = (1 - d) * t)
  L(t)  log of the risk ratio. A cut to r that builds up in a straight line (on the log scale) over D years:
            L(t) = ln(r) * min(t / D, 1)                         (D = 0: L(t) = ln r)
  z     hidden personal risk factor, exactly as in experiment 011 (gamma, variance V from experiment 004).

ASSUMPTIONS (no published figure was opened for any of these):
  - the straight-line build-up and its length D (0, 1 or 2 years);
  - "never given": a share p (0.1, 0.2, 0.3) of dogs assigned to treatment are never treated;
  - "stopped": each treated owner stops for good at a constant rate of 10 percent a year; after stopping at
    time u the aging clock runs at full speed again from where it stood: g(t) = (1 - d) * u + (t - u);
  - everything assumed in experiment 011 (dropout 5 percent a year unrelated to health, entry ages, size mix).

Analysis: by assigned group. The size-aware slope test of experiment 011 (likelihood ratio, 5 percent level).
"From year 2" analysis: the same test on dogs still followed 2 years after enrolment, each treated as entering
2 years older, with time counted from that point.

Run:  python3 delayed_partial.py            (about 5 minutes on 2 cores)
      python3 delayed_partial.py --quick
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


rt = _load("realistic", "011_realistic_trial", "realistic_trial.py")   # slope_test, frailty_variance, size code
br = rt.br

SEED = 20261118
AGE0, ENTRY, FOLLOW = rt.AGE0, rt.ENTRY, rt.FOLLOW      # 2.0, (7, 12), 10.0  as in experiment 011
DROP = 0.05                  # ASSUMPTION, experiment 011's "realistic" scenario
LANDMARK = 2.0               # the "from year 2" analysis
GRID = np.linspace(0.0, FOLLOW, 10001)

# name: (slowing d, risk ratio r, build-up years D, share never given, stopping rate per year)   ASSUMPTIONS
WORLDS = {
    "no effect": (0.0, 1.0, 0.0, 0.0, 0.0),
    "slowed 20%, at once": (0.2, 1.0, 0.0, 0.0, 0.0),
    "slowed 20%, builds over 1 year": (0.2, 1.0, 1.0, 0.0, 0.0),
    "slowed 20%, builds over 2 years": (0.2, 1.0, 2.0, 0.0, 0.0),
    "slowed 20%, 10% never given": (0.2, 1.0, 0.0, 0.1, 0.0),
    "slowed 20%, 20% never given": (0.2, 1.0, 0.0, 0.2, 0.0),
    "slowed 20%, 30% never given": (0.2, 1.0, 0.0, 0.3, 0.0),
    "slowed 20%, 10% stop each year": (0.2, 1.0, 0.0, 0.0, 0.1),
    "slowed 20%, builds over 1 year, 20% never given": (0.2, 1.0, 1.0, 0.2, 0.0),
    "risk cut 50%, at once": (0.0, 0.5, 0.0, 0.0, 0.0),
    "risk cut 50%, builds over 1 year": (0.0, 0.5, 1.0, 0.0, 0.0),
    "risk cut 50%, builds over 2 years": (0.0, 0.5, 2.0, 0.0, 0.0),
    "risk cut 50%, 20% never given": (0.0, 0.5, 0.0, 0.2, 0.0),
    "risk cut 50%, 30% never given": (0.0, 0.5, 0.0, 0.3, 0.0),
}
COMBINED = "slowed 20%, builds over 1 year, 20% never given"
WITH_LANDMARK = ("no effect", "slowed 20%, at once", "risk cut 50%, builds over 1 year",
                 "risk cut 50%, builds over 2 years")
# (world, dogs)
CASES = [(w, 1300) for w in WORLDS] + [(COMBINED, n) for n in (1600, 2000, 2600)]


def clock_and_level(t, d, r, D):
    """g(t) and L(t) of the header, for an array of times."""
    if D > 0:
        g = np.where(t <= D, t - d * t ** 2 / (2 * D), (1 - d) * t + d * D / 2)
        L = np.log(r) * np.minimum(t / D, 1.0)
    else:
        g = (1 - d) * t
        L = np.full_like(t, np.log(r))
    return g, L


def treated_time(target, b, cls, b_class, d, r, D):
    """Death time of a treated dog: solve  integral_0^T exp(b g(u) + L(u)) du = target  on a fine grid.
    Returns infinity when the answer lies beyond the end of follow-up."""
    g, L = clock_and_level(GRID, d, r, D)
    out = np.full(target.shape, np.inf)
    for c in range(len(b_class)):
        f = np.exp(b_class[c] * g + L)
        cum = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(GRID))])
        m = (cls == c) & (target < cum[-1])
        out[m] = np.interp(target[m], cum, GRID)
    return out


def simulate(n, world, a, b_class, shares, v, rng):
    d, r, D, never, stop = WORLDS[world]
    cls = rng.choice(len(shares), size=n, p=shares)
    x0 = rng.uniform(ENTRY[0], ENTRY[1], n)
    arm = (np.arange(n) % 2).astype(float)
    b = b_class[cls]
    H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
    z = rng.gamma(1.0 / v, v / (1.0 + v * H0))               # hidden risk factor among dogs alive at entry
    h0 = z * a * np.exp(b * x0)
    e = -np.log(rng.random(n))                                # the cumulative risk at which each dog dies
    t = np.log1p(e * b / h0) / b                              # untreated death time
    given = (arm == 1) & (rng.random(n) >= never)
    if stop > 0:                                              # immediate slowing until the owner stops at u
        u = rng.exponential(1.0 / stop, n)
        bs = b * (1 - d)
        h1 = h0 * np.expm1(bs * u) / bs                       # cumulative risk up to the stopping time
        before = np.log1p(e * bs / h0) / bs
        after = u + np.log1p(np.maximum(e - h1, 0.0) * b / (h0 * np.exp(bs * u))) / b
        t = np.where(given, np.where(e < h1, before, after), t)
    elif d > 0 or r < 1:
        t = np.where(given, treated_time(e / h0, b, cls, b_class, d, r, D), t)
    leave = rng.exponential(1.0 / DROP, n)
    end = np.minimum(leave, FOLLOW)
    died = (t <= end).astype(float)
    return x0, arm, np.minimum(t, end), died, cls


def from_landmark(data, cls):
    x0, arm, t, died = data
    keep = t > LANDMARK
    return (x0[keep] + LANDMARK, arm[keep], t[keep] - LANDMARK, died[keep]), cls[keep]


def one_case(args):
    world, n, nsim, a, b_class, shares, v, seed = args
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    b0 = float(np.average(b_class, weights=shares))
    lm = world in WITH_LANDMARK and n == 1300
    flag = flag_lm = 0
    s_hat, s_lm, seen, kept = [], [], [], []
    for _ in range(nsim):
        *data, cls5 = simulate(n, world, a, b_class, shares, v, rng)
        cls = br.STRATUM[cls5]
        stat, sh = rt.slope_test(data, cls, a, b0)
        flag += stat > c1
        s_hat.append(sh)
        seen.append(float(data[3].mean()))
        if lm:
            d2, c2 = from_landmark(data, cls)
            stat2, sh2 = rt.slope_test(d2, c2, a, b0)
            flag_lm += stat2 > c1
            s_lm.append(sh2)
            kept.append(len(c2) / n)
    row = dict(world=world, dogs=n, n_trials=nsim, slope_change_flagged=flag / nsim,
               median_slope_estimate=float(np.median(s_hat)), share_of_deaths_seen=float(np.mean(seen)))
    if lm:
        row.update(from_year_2_slope_change_flagged=flag_lm / nsim,
                   from_year_2_median_slope_estimate=float(np.median(s_lm)),
                   from_year_2_share_of_dogs_still_followed=float(np.mean(kept)))
    return row


def self_checks(a, b_class, shares, v):
    """The grid solver against the exact formula, with no build-up."""
    rng = np.random.default_rng(5)
    n = 4000
    cls = rng.choice(len(shares), size=n, p=shares)
    b = b_class[cls]
    target = rng.exponential(30.0, n)
    worst = 0.0
    for d, r in ((0.2, 1.0), (0.0, 0.5)):
        bs = b * (1 - d)
        exact = np.log1p(target * bs / r) / bs
        got = treated_time(target, b, cls, b_class, d, r, 0.0)
        m = exact < FOLLOW - 1e-6
        assert np.all(np.isinf(got[exact > FOLLOW + 1e-6]))
        worst = max(worst, float(np.max(np.abs(got[m] - exact[m]))))
    return worst


def main():
    quick = "--quick" in sys.argv
    nsim = 20 if quick else 800
    a, b_class, shares, names = br.size_classes()
    v = rt.frailty_variance()
    grid_err = self_checks(a, b_class, shares, v)
    assert grid_err < 1e-5, grid_err
    cases = [(w, n, nsim, a, b_class, shares, v, SEED + 1000 * i) for i, (w, n) in enumerate(CASES)]
    if quick:
        cases = [c for c in cases if c[0] in ("slowed 20%, at once", "risk cut 50%, builds over 2 years",
                                              "slowed 20%, 10% stop each year", COMBINED) and c[1] == 1300]
    cost = lambda c: c[1] * (2 if c[0] in WITH_LANDMARK else 1)
    order = sorted(range(len(cases)), key=lambda i: -cost(cases[i]))    # big cases first, to balance the cores
    with Pool(2) as pool:
        done = pool.map(one_case, [cases[i] for i in order], chunksize=1)
    rows = [None] * len(cases)
    for i, row in zip(order, done):
        rows[i] = row
    ok = [r["dogs"] for r in rows if r["world"] == COMBINED and r["slope_change_flagged"] >= 0.8]
    out = dict(
        note="Simulated 10-year trials in the 'realistic' scenario of experiment 011 (size mix, hidden risk factor, "
             "5 percent of dogs lost a year, entry ages 7 to 12). Dogs are analysed in the group they were assigned "
             "to. 'slope_change_flagged' is the share of trials in which the size-aware slope test reports a change "
             "in the rate of aging; 'from_year_2' is the same test on dogs still followed 2 years after enrolment. "
             "ASSUMPTIONS: straight-line build-up over 1 or 2 years; shares of dogs never given the treatment; "
             "owners stopping at 10 percent a year with the gain kept; and all assumptions of experiment 011.",
        seed=SEED, trials_per_case=nsim, frailty_variance_from_experiment_004=v, dropout_per_year=DROP,
        size_classes=names, rate_of_aging_per_year=[float(x) for x in b_class],
        population_share=[float(x) for x in shares], common_level_a=a,
        grid_solver_max_error_years=grid_err,
        smallest_grid_size_reaching_80_percent_combined_case=min(ok) if ok else None, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        extra = (f" | from year 2: flagged {r['from_year_2_slope_change_flagged']:.3f} "
                 f"s-hat {r['from_year_2_median_slope_estimate']:.2f}" if "from_year_2_slope_change_flagged" in r else "")
        print(f"{r['world']:48s} | {r['dogs']:5d} | flagged {r['slope_change_flagged']:.3f} | "
              f"s-hat {r['median_slope_estimate']:.3f} | deaths seen {r['share_of_deaths_seen']:.2f}{extra}")
    print("smallest for the combined case:", out["smallest_grid_size_reaching_80_percent_combined_case"])


if __name__ == "__main__":
    main()
