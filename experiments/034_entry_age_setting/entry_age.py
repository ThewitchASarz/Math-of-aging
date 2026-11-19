"""Experiment 034: does giving entry age its own setting repair the score-adjusted slope test?

Model that makes the data (experiment 032's, unchanged; simulate() is imported from it):
    hazard(t) = z * a * exp(b * x0) * r * exp(b * s * t)
  z     hidden personal risk factor: gamma, mean 1 and variance V = 0.3 among dogs alive at age 2 (ASSUMPTION).
  r, s  treatment: r < 1 cuts risk, s < 1 slows aging; both 1 in the control group.
  a, b  UK fit, dog_lifespan/results.json (read through experiment 001).
Trial: 1,300 dogs, half treated, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

Frailty scores (ASSUMPTION):
  rank score   w = RHO * q + sqrt(1 - RHO^2) * noise, q = bell-curve number with the dog's rank in z among
               dogs of its entry age (experiment 032's score)
  risk score   the same with q replaced by ln z, put on a standard scale across the dogs in the trial

True effects (sizes are ASSUMPTIONS): no effect; "sturdy helped more": treated dogs with z below the median
for their entry age get r = 0.3, the others r = 0.7; aging slowed: s = 0.8 for all.

Analyses. Each fits   ln hazard = ln a + [entry-age part] + X . beta + b * s^arm * t   with no hidden
differences, and tests s = 1 with a 5% likelihood-ratio test (1 degree of freedom).
  old          entry-age part = b * x0 (one rate for entry age and time);  X = [w, arm, arm * w]
  entry-age    entry-age part = c * (x0 - 9.5), c its own setting;         X = [x0 - 9.5, w, arm, arm * w]
  five groups  as entry-age, but X = [x0 - 9.5, w, arm * (score in fifth 1), ..., arm * (score in fifth 5)]
Log-likelihood: sum of died * (A + b s t) - exp(A) * (exp(b s t) - 1) / (b s),  A = ln hazard at t = 0.

Run:  python3 entry_age.py            (about 5 minutes on 2 cores)
      python3 entry_age.py --quick
"""
import importlib.util
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, norm
from scipy.stats import gamma as gamma_dist

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "001_trial_fingerprint", "fingerprint.py")
e32 = _load("frailty_score", "032_frailty_score", "frailty_score.py")

SEED = 20261134
V, AGE0 = e32.V, e32.AGE0                   # spread 0.3: ASSUMPTION
MID = 9.5                                   # middle of the entry ages, used only to centre x0
# (kind of score, quality)   ASSUMPTION
SCORES_ALL = (("rank", 0.5), ("rank", 0.7), ("rank", 0.9), ("rank", 1.0), ("risk", 0.7), ("risk", 1.0))
SCORES_FEW = (("rank", 0.7), ("rank", 1.0), ("risk", 0.7), ("risk", 1.0))
# name: (r for sturdier half, r for frailer half, s, scores analysed)   all sizes: ASSUMPTION
CASES = {
    "no effect": (1.0, 1.0, 1.0, SCORES_ALL),
    "sturdy helped more (70% and 30% cuts)": (0.3, 0.7, 1.0, SCORES_ALL),
    "aging slowed 20% in every dog": (1.0, 1.0, 0.8, SCORES_FEW),
}
ANALYSES = ("old", "entry-age", "five groups")


def negll_and_grad(theta, x0, t, died, arm, X, tied):
    """theta = [ln a, ln b, beta (one per column of X), ln s]. Minus log-likelihood and its gradient."""
    lna, lnb, lns = theta[0], theta[1], theta[-1]
    b = np.exp(lnb)
    s = np.exp(lns * arm)
    A = lna + X @ theta[2:-1]
    if tied:
        A = A + b * x0
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    H = eA * (E - 1.0) / bs
    ll = np.sum(died * (A + bs * t) - H)
    res = died - H
    dH_b = eA * t * E / b - H / b            # dH/db apart from the entry-age part
    g_b = died * s * t - dH_b
    if tied:
        g_b = g_b + res * x0
    grad = np.concatenate(([np.sum(res)], [b * np.sum(g_b)], X.T @ res,
                           [np.sum(arm * (died * bs * t - eA * t * E + H))]))
    return -ll, -grad


def fit(data, X, tied, start, free_s, a, b):
    x0, arm, t, died = data
    k = X.shape[1]
    bounds = ([(np.log(a) - 6, np.log(a) + 8), (np.log(b) - 1, np.log(b) + 1)] + [(-4.0, 4.0)] * k
              + [(np.log(0.05), np.log(3.0)) if free_s else (start[-1], start[-1])])
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, np.array(start, float), args=(x0, t, died, arm, X, tied), jac=True,
                       method="L-BFGS-B", bounds=bounds, options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-7})
    return -res.fun, res.x


def slope_test(data, X, tied, a, b, start_beta):
    """Level-only fit, then the fit with a free slope change from two starts. Returns statistic, s-hat, level fit."""
    lna0 = np.log(a) + (0.0 if tied else b * MID)
    lvl = fit(data, X, tied, [lna0, np.log(b), *start_beta, 0.0], False, a, b)
    full = max((fit(data, X, tied, lvl[1], True, a, b),
                fit(data, X, tied, [*lvl[1][:-1], np.log(0.8)], True, a, b)), key=lambda f: f[0])
    return 2 * (full[0] - lvl[0]), float(np.exp(full[1][-1])), lvl


def analyse(data, w, a, b):
    """The three slope tests and two tests of 'the benefit differs by score'."""
    x0, arm, t, died = data
    xc = x0 - MID
    out = {}
    out["old"] = slope_test(data, np.column_stack([w, arm, arm * w]), True, a, b, [0.0, 0.0, 0.0])[:2]
    X_line = np.column_stack([xc, w, arm, arm * w])
    stat, s_hat, lvl_line = slope_test(data, X_line, False, a, b, [b, 0.0, 0.0, 0.0])
    out["entry-age"] = (stat, s_hat)
    fifth = np.searchsorted(np.quantile(w, [0.2, 0.4, 0.6, 0.8]), w)
    X_five = np.column_stack([xc, w] + [arm * (fifth == k) for k in range(5)])
    stat, s_hat, lvl_five = slope_test(data, X_five, False, a, b, [b, 0.0] + [0.0] * 5)
    out["five groups"] = (stat, s_hat)
    # does the benefit differ by score?  Same model with one risk cut for every treated dog as the base.
    same = fit(data, np.column_stack([xc, w, arm]), False, [np.log(a) + b * MID, np.log(b), b, 0.0, 0.0, 0.0],
               False, a, b)
    return out, 2 * (lvl_line[0] - same[0]), 2 * (lvl_five[0] - same[0])


def make_score(kind, rho, data, q, noise, a, b):
    if kind == "rank":
        base = q
    else:
        x0 = data[0]
        H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
        lnz = np.log(gamma_dist.ppf(norm.cdf(q), 1.0 / V) * V / (1.0 + V * H0))
        base = (lnz - lnz.mean()) / lnz.std()
    return rho * base + np.sqrt(1.0 - rho * rho) * noise


def check_gradient(a, b):
    rng = np.random.default_rng(5)
    data, q, noise = e32.simulate(0.4, 0.8, 0.9, a, b, rng)
    x0, arm, t, died = data
    w = 0.7 * q + np.sqrt(1 - 0.49) * noise
    worst = 0.0
    for tied, X in ((True, np.column_stack([w, arm, arm * w])),
                    (False, np.column_stack([x0 - MID, w, arm, arm * w]))):
        theta = np.concatenate(([np.log(a) + (0.1 if tied else b * MID), np.log(b) + 0.05],
                                np.linspace(0.3, -0.2, X.shape[1]), [-0.1]))
        _, g = negll_and_grad(theta, x0, t, died, arm, X, tied)
        for i in range(theta.size):
            d = np.zeros(theta.size)
            d[i] = 1e-6
            num = (negll_and_grad(theta + d, x0, t, died, arm, X, tied)[0]
                   - negll_and_grad(theta - d, x0, t, died, arm, X, tied)[0]) / 2e-6
            worst = max(worst, abs(g[i] - num) / (abs(num) + 1e-9))
    return float(worst)


def check_against_experiment_032(a, b):
    """The 'old' likelihood here must equal experiment 032's with no hidden differences."""
    rng = np.random.default_rng(6)
    data, q, noise = e32.simulate(0.5, 0.5, 0.9, a, b, rng)
    x0, arm, t, died = data
    w = 0.7 * q + np.sqrt(1 - 0.49) * noise
    th = [np.log(a) + 0.1, np.log(b) - 0.05, 0.3, -0.3, 0.15, -0.1]       # ln a, ln b, g, ln r, d, ln s
    mine = negll_and_grad(np.array(th), x0, t, died, arm, np.column_stack([w, arm, arm * w]), True)[0]
    theirs = e32.negll_and_grad(np.array([th[0], th[1], np.log(e32.V_TINY), th[2], th[3], th[4], th[5]]),
                                *data, w)[0]
    return float(abs(mine - theirs) / abs(theirs))


def one_chunk(args):
    name, nsim, a, b, seed = args
    r_sturdy, r_frail, s, scores = CASES[name]
    rng = np.random.default_rng(seed)
    c1, c4 = chi2.ppf(0.95, 1), chi2.ppf(0.95, 4)
    out = {(sc, an): dict(flag=0, slower=0, s=[]) for sc in scores for an in ANALYSES}
    differs = {sc: [0, 0] for sc in scores}
    for _ in range(nsim):
        data, q, noise = e32.simulate(r_sturdy, r_frail, s, a, b, rng)
        for sc in scores:
            w = make_score(sc[0], sc[1], data, q, noise, a, b)
            tests, d_line, d_five = analyse(data, w, a, b)
            differs[sc][0] += int(d_line > c1)
            differs[sc][1] += int(d_five > c4)
            for an, (stat, s_hat) in tests.items():
                o = out[(sc, an)]
                o["flag"] += int(stat > c1)
                o["slower"] += int(stat > c1 and s_hat < 1.0)
                o["s"].append(s_hat)
    return name, nsim, out, differs


def main():
    quick = "--quick" in sys.argv
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    same = check_against_experiment_032(a, b)
    assert same < 1e-6, same
    chunks, per = (1, 6) if quick else (4, 150)
    jobs = [(name, per, a, b, SEED + 1000 * i + k) for i, name in enumerate(CASES) for k in range(chunks)]
    with Pool(2) as pool:
        parts = pool.map(one_chunk, jobs, chunksize=1)
    rows = []
    for name, (r_sturdy, r_frail, s, scores) in CASES.items():
        mine = [p for p in parts if p[0] == name]
        n = sum(p[1] for p in mine)
        for sc in scores:
            for an in ANALYSES:
                flag = sum(p[2][(sc, an)]["flag"] for p in mine)
                slower = sum(p[2][(sc, an)]["slower"] for p in mine)
                row = dict(case=name, score_kind=sc[0], score_quality=sc[1], analysis=an, n_trials=n,
                           risk_multiplier_sturdier_half=r_sturdy, risk_multiplier_frailer_half=r_frail,
                           rate_multiplier=s, slope_change_flagged=flag / n, flagged_as_slower_aging=slower / n,
                           flagged_as_faster_aging=(flag - slower) / n,
                           median_slope_estimate=float(np.median([x for p in mine for x in p[2][(sc, an)]["s"]])))
                if an != "old":
                    row["benefit_differs_by_score_found"] = sum(p[3][sc][an == "five groups"] for p in mine) / n
                rows.append(row)
    out = dict(note="Simulated 10-year trials of 1,300 dogs, half treated from enrolment, analysed with a frailty "
                    "score measured at entry. 'slope_change_flagged' is the share of trials in which a 5% "
                    "likelihood-ratio test rejects 'the slope did not change'. 'old' is experiment 032's plain "
                    "test with score terms; 'entry-age' gives entry age its own setting; 'five groups' also "
                    "replaces the straight-line 'benefit differs by score' term with a risk cut per fifth of "
                    "the score. ASSUMPTIONS: spread 0.3, entry ages 7 to 12, the form and quality of the scores, "
                    "and every size of benefit.",
               seed=SEED, uk_a_per_year=a, uk_b_per_year=b, spread=V, n_dogs=e32.N,
               entry_ages=[e32.LO, e32.HI], years=e32.FOLLOW, gradient_check_max_relative_error=grad_err,
               difference_from_experiment_032_likelihood=same, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case'][:30]:30s} | {r['score_kind']} {r['score_quality']:.1f} | {r['analysis']:11s}"
              f" | flagged {r['slope_change_flagged']:.3f} (slower {r['flagged_as_slower_aging']:.3f},"
              f" faster {r['flagged_as_faster_aging']:.3f}) | s-hat {r['median_slope_estimate']:.3f}"
              + (f" | differs found {r['benefit_differs_by_score_found']:.2f}"
                 if "benefit_differs_by_score_found" in r else ""))


if __name__ == "__main__":
    main()
