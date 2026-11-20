"""Experiment 035: a score-adjusted slope test with a setting for left-over hidden differences.

Model that makes the data (experiment 032's, unchanged; simulate() is imported from it):
    hazard(t) = z * a * exp(b * x0) * r * exp(b * s * t)
  z     hidden personal risk factor: gamma, mean 1 and variance V = 0.3 among dogs alive at age 2 (ASSUMPTION).
  r, s  treatment: r < 1 cuts risk, s < 1 slows aging; both 1 in the control group.
  a, b  UK fit, dog_lifespan/results.json (read through experiment 001).
Trial: 1,300 dogs, half treated, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

Frailty scores (ASSUMPTION; experiment 034's make_score, imported): rank score and risk score, quality RHO.
"none" means the analysis has no score terms.

True effects (sizes are ASSUMPTIONS): no effect; "sturdy helped more": treated dogs with z below the median
for their entry age get r = 0.3, the others r = 0.7; r = 0.5 for all; aging slowed: s = 0.8 for all.

Analyses. Both test s = 1 with a 5% likelihood-ratio test (1 degree of freedom).
  entry-age              experiment 034's (imported):  ln hazard = ln a + X . beta + b * s^arm * t,
                         X = [x0 - 9.5, w, arm, arm * w]   (X = [x0 - 9.5, arm] with no score)
  entry-age with spread  the same with one more setting v. Each dog's hazard is multiplied by a left-over hidden
                         factor, gamma with mean 1 and variance v among dogs at entry. With
                             A = ln a + X . beta,   L(t) = exp(A) * (exp(b s t) - 1) / (b s)
                         the chance of surviving to t is (1 + v L(t)) ** (-1/v) and
                             ln lik = died * (A + b s t - ln(1 + v L)) - ln(1 + v L) / v
                         No term for deaths before entry: the entry-age setting stands in for it.

Run:  python3 score_and_spread.py            (about 5 minutes on 2 cores)
      python3 score_and_spread.py --quick
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
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "001_trial_fingerprint", "fingerprint.py")
e32 = _load("frailty_score", "032_frailty_score", "frailty_score.py")
e34 = _load("entry_age", "034_entry_age_setting", "entry_age.py")

SEED = 20261135
MID = e34.MID
V_LO, V_HI = 1e-4, 5.0                      # limits on the fitted left-over spread
NONE = ("none", 0.0)
# (kind of score, quality)   ASSUMPTION
S_ALL = (NONE, ("rank", 0.5), ("rank", 0.7), ("rank", 0.9), ("rank", 1.0), ("risk", 0.7), ("risk", 1.0))
# name: (r for sturdier half, r for frailer half, s, scores analysed)   all sizes: ASSUMPTION
CASES = {
    "no effect": (1.0, 1.0, 1.0, (NONE, ("rank", 0.7), ("rank", 0.9), ("rank", 1.0), ("risk", 0.7), ("risk", 1.0))),
    "sturdy helped more (70% and 30% cuts)": (0.3, 0.7, 1.0, S_ALL),
    "risk cut 50% in every dog": (0.5, 0.5, 1.0, (NONE, ("rank", 0.7), ("risk", 0.7))),
    "aging slowed 20% in every dog": (1.0, 1.0, 0.8, (NONE, ("rank", 0.7), ("rank", 1.0), ("risk", 0.7))),
}
ANALYSES = ("entry-age", "entry-age with spread")


def negll_and_grad(theta, t, died, arm, X):
    """theta = [ln a, ln b, ln v, beta (one per column of X), ln s]. Minus log-likelihood and its gradient."""
    lna, lnb, lnv, lns = theta[0], theta[1], theta[2], theta[-1]
    b, v = np.exp(lnb), np.exp(lnv)
    s = np.exp(lns * arm)
    A = lna + X @ theta[3:-1]
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    L = eA * (E - 1.0) / bs
    D = 1.0 + v * L
    lnD = np.log1p(v * L)
    ll = np.sum(died * (A + bs * t - lnD) - lnD / v)
    c = (1.0 + v * died) / D
    res = died - c * L
    dL = eA * t * E - L                      # dL / d ln b, and dL / d ln s for treated dogs
    slope = died * bs * t - c * dL
    grad = np.concatenate(([np.sum(res)], [np.sum(slope)],
                           [np.sum(-died * v * L / D + lnD / v - L / D)],
                           X.T @ res, [np.sum(arm * slope)]))
    return -ll, -grad


def fit(data, X, start, free_s, a, b):
    x0, arm, t, died = data
    k = X.shape[1]
    bounds = ([(np.log(a) - 6, np.log(a) + 8), (np.log(b) - 1, np.log(b) + 1), (np.log(V_LO), np.log(V_HI))]
              + [(-4.0, 4.0)] * k + [(np.log(0.05), np.log(3.0)) if free_s else (start[-1], start[-1])])
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, np.array(start, float), args=(t, died, arm, X), jac=True,
                       method="L-BFGS-B", bounds=bounds, options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-7})
    return -res.fun, res.x


def best(fits):
    return max(fits, key=lambda f: f[0])


def with_v(p, v0):
    """Settings of a no-spread fit [ln a, ln b, beta, ln s] with a starting spread put in."""
    return [p[0], p[1], np.log(v0), *p[2:]]


def analyse(data, w, a, b):
    """Both slope tests. Returns {analysis: (statistic, s-hat)} and the fitted spread without a slope change."""
    x0, arm, t, died = data
    xc = x0 - MID
    if w is None:
        X, start_beta = np.column_stack([xc, arm]), [b, 0.0]
    else:
        X, start_beta = np.column_stack([xc, w, arm, arm * w]), [b, 0.0, 0.0, 0.0]
    stat, s_hat, lvl = e34.slope_test(data, X, False, a, b, start_beta)
    out = {"entry-age": (stat, s_hat)}
    n_lvl = best([fit(data, X, with_v(lvl[1], v0), False, a, b) for v0 in (0.02, 0.5)])
    q = n_lvl[1]
    n_full = best([fit(data, X, q, True, a, b),
                   fit(data, X, [*q[:-1], np.log(0.8)], True, a, b),
                   fit(data, X, [*q[:-1], np.log(1.2)], True, a, b)])
    out["entry-age with spread"] = (2 * (n_full[0] - n_lvl[0]), float(np.exp(n_full[1][-1])))
    return out, float(np.exp(q[2]))


def check_gradient(a, b):
    rng = np.random.default_rng(5)
    data, q, noise = e32.simulate(0.4, 0.8, 0.9, a, b, rng)
    x0, arm, t, died = data
    w = 0.7 * q + np.sqrt(1 - 0.49) * noise
    X = np.column_stack([x0 - MID, w, arm, arm * w])
    theta = np.concatenate(([np.log(a) + b * MID, np.log(b) + 0.05, np.log(0.2)], [0.3, 0.2, -0.3, 0.1], [-0.1]))
    _, g = negll_and_grad(theta, t, died, arm, X)
    worst = 0.0
    for i in range(theta.size):
        d = np.zeros(theta.size)
        d[i] = 1e-6
        num = (negll_and_grad(theta + d, t, died, arm, X)[0] - negll_and_grad(theta - d, t, died, arm, X)[0]) / 2e-6
        worst = max(worst, abs(g[i] - num) / (abs(num) + 1e-9))
    return float(worst)


def check_against_experiment_034(a, b):
    """With the spread at its lower limit this likelihood must be close to experiment 034's entry-age one."""
    rng = np.random.default_rng(6)
    data, q, noise = e32.simulate(0.5, 0.5, 0.9, a, b, rng)
    x0, arm, t, died = data
    w = 0.7 * q + np.sqrt(1 - 0.49) * noise
    X = np.column_stack([x0 - MID, w, arm, arm * w])
    th = [np.log(a) + b * MID, np.log(b) - 0.05, 0.3, 0.2, -0.3, 0.15, -0.1]
    theirs = e34.negll_and_grad(np.array(th), x0, t, died, arm, X, False)[0]
    mine = negll_and_grad(np.array(with_v(th, 1e-9)), t, died, arm, X)[0]
    return float(abs(mine - theirs) / abs(theirs))


def one_chunk(args):
    name, nsim, a, b, seed = args
    r_sturdy, r_frail, s, scores = CASES[name]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    out = {(sc, an): dict(flag=0, slower=0, s=[]) for sc in scores for an in ANALYSES}
    spread = {sc: [] for sc in scores}
    for _ in range(nsim):
        data, q, noise = e32.simulate(r_sturdy, r_frail, s, a, b, rng)
        for sc in scores:
            w = None if sc == NONE else e34.make_score(sc[0], sc[1], data, q, noise, a, b)
            tests, v_hat = analyse(data, w, a, b)
            spread[sc].append(v_hat)
            for an, (stat, s_hat) in tests.items():
                o = out[(sc, an)]
                o["flag"] += int(stat > c1)
                o["slower"] += int(stat > c1 and s_hat < 1.0)
                o["s"].append(s_hat)
    return name, nsim, out, spread


def main():
    quick = "--quick" in sys.argv
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    same = check_against_experiment_034(a, b)
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
            v_all = np.array([x for p in mine for x in p[3][sc]])
            for an in ANALYSES:
                flag = sum(p[2][(sc, an)]["flag"] for p in mine)
                slower = sum(p[2][(sc, an)]["slower"] for p in mine)
                row = dict(case=name, score_kind=sc[0], score_quality=sc[1], analysis=an, n_trials=n,
                           risk_multiplier_sturdier_half=r_sturdy, risk_multiplier_frailer_half=r_frail,
                           rate_multiplier=s, slope_change_flagged=flag / n, flagged_as_slower_aging=slower / n,
                           flagged_as_faster_aging=(flag - slower) / n,
                           median_slope_estimate=float(np.median([x for p in mine for x in p[2][(sc, an)]["s"]])))
                if an == "entry-age with spread":
                    row["median_fitted_spread"] = float(np.median(v_all))
                    row["share_fitted_spread_at_lower_limit"] = float(np.mean(v_all < 2 * V_LO))
                rows.append(row)
    out = dict(note="Simulated 10-year trials of 1,300 dogs, half treated from enrolment, analysed with a frailty "
                    "score measured at entry ('none': no score terms). 'slope_change_flagged' is the share of "
                    "trials in which a 5% likelihood-ratio test rejects 'the slope did not change'. 'entry-age' "
                    "is experiment 034's test; 'entry-age with spread' adds one setting for hidden differences "
                    "left over after the score (gamma, mean 1, among dogs at entry). 'median_fitted_spread' is "
                    "that setting in the fit without a slope change. ASSUMPTIONS: spread 0.3, entry ages 7 to 12, "
                    "the form and quality of the scores, and every size of benefit.",
               seed=SEED, uk_a_per_year=a, uk_b_per_year=b, spread=e32.V, n_dogs=e32.N,
               entry_ages=[e32.LO, e32.HI], years=e32.FOLLOW, gradient_check_max_relative_error=grad_err,
               difference_from_experiment_034_likelihood=same, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case'][:30]:30s} | {r['score_kind']} {r['score_quality']:.1f} | {r['analysis']:21s}"
              f" | flagged {r['slope_change_flagged']:.3f} (slower {r['flagged_as_slower_aging']:.3f},"
              f" faster {r['flagged_as_faster_aging']:.3f}) | s-hat {r['median_slope_estimate']:.3f}"
              + (f" | spread {r['median_fitted_spread']:.2f}" if "median_fitted_spread" in r else ""))


if __name__ == "__main__":
    main()
