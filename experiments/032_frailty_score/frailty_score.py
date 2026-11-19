"""Experiment 032: can a frailty score measured at entry stop a benefit that differs between dogs from being
read as slowed aging?

Model for one dog entering at age x0, t years after enrolment (experiment 017, Part 3):
    hazard(t) = z * a * exp(b * x0) * r * exp(b * s * t)
  z     hidden personal risk factor: gamma with mean 1 and variance V among dogs alive at age 2, so among dogs
        alive at entry age x0 it is gamma with shape 1/V and scale V / (1 + V * H0),
        H0 = (a / b) * (exp(b * x0) - exp(2 * b)).   V = 0.3 is an ASSUMPTION (no measurement exists).
  r, s  treatment: r < 1 cuts risk, s < 1 slows aging; both 1 in the control group.
  a, b  UK fit, dog_lifespan/results.json (read through experiment 001).

Trial (experiment 012's design): 1,300 dogs, half treated, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

Frailty score (ASSUMPTION): q = standard normal number with the same rank as the dog's z among dogs of its
entry age; score w = RHO * q + sqrt(1 - RHO^2) * noise. RHO is the quality of the score.

True effects: "sturdy helped more": treated dogs with z below the median for their entry age get r = 0.3, the
others r = 0.7 (ASSUMPTION, experiment 017); "frail helped more" swaps them; no effect; r = 0.5 for all;
s = 0.8 for all.

Analyses, each a 5% likelihood-ratio test of "no slope change" (s = 1) against "slope change":
  plain, new            experiment 017's two tests (imported, unchanged): no score
  plain+score, new+score  the same with two more settings: ln hazard gets g * w for every dog and d * w for
                        treated dogs ("benefit differs by score").
  Likelihood with score: A = ln a + b*x0 + g*w + (ln r + d*w)*arm,  L = exp(A) (exp(b*s*t) - 1) / (b*s),
        H0w = exp(g*w) * H0,   ln lik = died * (A + b*s*t - ln(1 + v*H0w + v*L)) - ln(1 + v*L/(1 + v*H0w)) / v
  "plain+score" fixes v at 1e-8 (no hidden differences); "new+score" estimates v.
Also: a 5% test of d = 0 in the plain+score model without a slope change ("benefit differs by score" found).

Run:  python3 frailty_score.py            (about 5 minutes on 2 cores)
      python3 frailty_score.py --quick
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
h17 = _load("harder", "017_hidden_differences_harder", "harder.py")

SEED = 20261119
AGE0 = 2.0
N, LO, HI, FOLLOW = 1300, 7.0, 12.0, 10.0   # entry ages: ASSUMPTION
V = 0.3                                     # ASSUMPTION: spread of hidden risk among dogs alive at age 2
V_TINY = 1e-8                               # "no hidden differences" in the plain+score analysis
V_LO, V_HI = h17.V_LO, h17.V_HI
ALL_RHO = (0.3, 0.5, 0.7, 0.9, 1.0)         # ASSUMPTION: quality of the score
# name: (r for sturdier half, r for frailer half, s, score qualities analysed)   all sizes: ASSUMPTION
CASES = {
    "sturdy helped more (70% and 30% cuts)": (0.3, 0.7, 1.0, ALL_RHO),
    "frail helped more (30% and 70% cuts)": (0.7, 0.3, 1.0, (0.7,)),
    "no effect": (1.0, 1.0, 1.0, (0.7, 0.9, 1.0)),          # 0.9 and 1.0 added after the first run, see RESULT.md
    "risk cut 50% in every dog": (0.5, 0.5, 1.0, (0.7,)),
    "aging slowed 20% in every dog": (1.0, 1.0, 0.8, (0.7, 1.0)),   # 1.0 added after the first run
}


def simulate(r_sturdy, r_frail, s, a, b, rng):
    x0 = rng.uniform(LO, HI, N)
    arm = (np.arange(N) % 2).astype(float)
    H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
    u = rng.gamma(1.0 / V, 1.0, N)                    # the dog's standing among dogs of its entry age
    z = u * V / (1.0 + V * H0)
    p = np.clip(gamma_dist.cdf(u, 1.0 / V), 1e-12, 1 - 1e-12)
    q = norm.ppf(p)                                   # the same rank as a standard normal number
    frail = p > 0.5
    r = np.where(arm == 1, np.where(frail, r_frail, r_sturdy), 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(N)) * bs / (z * a * np.exp(b * x0) * r)) / bs
    died = (t <= FOLLOW).astype(float)
    noise = rng.standard_normal(N)
    return (x0, arm, np.minimum(t, FOLLOW), died), q, noise


def negll_and_grad(theta, x0, arm, t, died, w):
    """theta = [ln a, ln b, ln v, g, ln r, d, ln s]. Returns minus log-likelihood and its gradient."""
    lna, lnb, lnv, g, lnr, d, lns = theta
    a, b, v = np.exp(lna), np.exp(lnb), np.exp(lnv)
    s = np.exp(lns * arm)
    A = lna + b * x0 + g * w + (lnr + d * w) * arm
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    L = eA * (E - 1.0) / bs
    egw = np.exp(g * w)
    ebx, e2b = np.exp(b * x0), np.exp(b * AGE0)
    H0 = egw * (a / b) * (ebx - e2b)
    D0 = 1.0 + v * H0
    D = D0 + v * L
    lnratio = np.log1p(v * L / D0)
    ll = np.sum(died * (A + bs * t - np.log(D)) - lnratio / v)
    c = (1.0 + v * died) / D
    dL_b = b * x0 * L + eA * t * E - L
    dH0_b = -H0 + egw * a * (x0 * ebx - AGE0 * e2b)
    base = died - c * (H0 + L) + H0 / D0
    trt = arm * (died - c * L)
    grad = np.array([np.sum(base),
                     np.sum(died * b * (x0 + s * t) - c * (dH0_b + dL_b) + dH0_b / D0),
                     np.sum(-died * v * (H0 + L) / D + lnratio / v - (H0 + L) / D + H0 / D0),
                     np.sum(w * base),
                     np.sum(trt),
                     np.sum(w * trt),
                     np.sum(arm * (died * bs * t - c * (eA * t * E - L)))])
    return -ll, -grad


def fit(data, w, a, b, start, free_v, free_d, free_s):
    """Maximum likelihood. A setting that is not free is held at its start value (equal bounds)."""
    bounds = [(np.log(a) - 4, np.log(a) + 4), (np.log(b) - 1, np.log(b) + 1),
              (np.log(V_LO), np.log(V_HI)) if free_v else (start[2], start[2]),
              (-3.0, 3.0), (np.log(0.05), np.log(5.0)),
              (-3.0, 3.0) if free_d else (start[5], start[5]),
              (np.log(0.05), np.log(3.0)) if free_s else (start[6], start[6])]
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, np.array(start, float), args=(*data, w), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-7})
    return -res.fun, res.x


def best(fits):
    return max(fits, key=lambda f: f[0])


def score_tests(data, w, a, b):
    """Slope tests with the score terms. Returns {analysis: (statistic, s-hat)}, statistic for d = 0, d-hat, g-hat."""
    tiny = np.log(V_TINY)
    lna, lnb = np.log(a), np.log(b)
    p_nod = fit(data, w, a, b, [lna, lnb, tiny, 0.0, 0.0, 0.0, 0.0], False, False, False)
    p_lvl = fit(data, w, a, b, p_nod[1], False, True, False)
    q = p_lvl[1]
    p_full = best([fit(data, w, a, b, q, False, True, True),
                   fit(data, w, a, b, [*q[:6], np.log(0.8)], False, True, True)])
    n_lvl = best([fit(data, w, a, b, [q[0], q[1], np.log(v0), *q[3:]], True, True, False) for v0 in (0.02, 0.5)])
    f = p_full[1]
    n_full = best([fit(data, w, a, b, n_lvl[1], True, True, True),
                   fit(data, w, a, b, [f[0], f[1], n_lvl[1][2], *f[3:]], True, True, True)])
    return ({"plain+score": (2 * (p_full[0] - p_lvl[0]), float(np.exp(f[6]))),
             "new+score": (2 * (n_full[0] - n_lvl[0]), float(np.exp(n_full[1][6])))},
            2 * (p_lvl[0] - p_nod[0]), float(q[5]), float(q[3]))


def check_gradient(a, b):
    rng = np.random.default_rng(5)
    data, q, noise = simulate(0.4, 0.8, 0.9, a, b, rng)
    w = 0.7 * q + np.sqrt(1 - 0.49) * noise
    theta = np.array([np.log(a) + 0.1, np.log(b) + 0.05, np.log(0.2), 0.3, -0.2, 0.15, -0.1])
    _, g = negll_and_grad(theta, *data, w)
    worst = 0.0
    for i in range(theta.size):
        d = np.zeros(theta.size)
        d[i] = 1e-6
        num = (negll_and_grad(theta + d, *data, w)[0] - negll_and_grad(theta - d, *data, w)[0]) / 2e-6
        worst = max(worst, abs(g[i] - num) / (abs(num) + 1e-9))
    return float(worst)


def check_against_experiment_001(a, b):
    """With no score terms and v tiny, this likelihood must equal experiment 001's."""
    rng = np.random.default_rng(6)
    data, q, noise = simulate(0.5, 0.5, 0.9, a, b, rng)
    theta = [np.log(a) + 0.1, np.log(b) - 0.05, np.log(V_TINY), 0.0, -0.3, 0.0, -0.1]
    mine = negll_and_grad(np.array(theta), *data, q)[0]
    theirs = fp.negll_and_grad(np.array([theta[0], theta[1], theta[4], theta[6]]), *data, True, True)[0]
    return float(abs(mine - theirs) / abs(theirs))


def one_chunk(args):
    name, nsim, a, b, seed = args
    r_sturdy, r_frail, s, rhos = CASES[name]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    keys = [("plain", None), ("new", None)] + [(k, rho) for rho in rhos for k in ("plain+score", "new+score")]
    out = {k: dict(flag=0, slower=0, s=[]) for k in keys}
    inter = {rho: dict(found=0, d=[], g=[]) for rho in rhos}
    dead = 0.0
    for _ in range(nsim):
        data, q, noise = simulate(r_sturdy, r_frail, s, a, b, rng)
        dead += data[3].mean()
        res = {(k, None): v[:2] for k, v in h17.tests_one_size(data, a, b).items()}
        for rho in rhos:
            w = rho * q + np.sqrt(1.0 - rho * rho) * noise
            tests, stat_d, d_hat, g_hat = score_tests(data, w, a, b)
            for k, v in tests.items():
                res[(k, rho)] = v
            inter[rho]["found"] += int(stat_d > c1)
            inter[rho]["d"].append(d_hat)
            inter[rho]["g"].append(g_hat)
        for k, (stat, s_hat) in res.items():
            out[k]["flag"] += int(stat > c1)
            out[k]["slower"] += int(stat > c1 and s_hat < 1.0)
            out[k]["s"].append(s_hat)
    return name, nsim, dead, out, inter


def one_extra(args):
    """Added after the first run (see RESULT.md): the score is the dog's hidden risk itself, ln z, read
    without error and put on a standard scale, in place of its rank among dogs of the same entry age."""
    name, nsim, a, b, seed = args
    r_sturdy, r_frail, s, _ = CASES[name]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    flag = slower = found = 0
    s_all = []
    for _ in range(nsim):
        data, q, _noise = simulate(r_sturdy, r_frail, s, a, b, rng)
        x0 = data[0]
        H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
        lnz = np.log(gamma_dist.ppf(norm.cdf(q), 1.0 / V) * V / (1.0 + V * H0))
        w = (lnz - lnz.mean()) / lnz.std()
        tests, stat_d, _d, _g = score_tests(data, w, a, b)
        stat, s_hat = tests["plain+score"]
        flag += int(stat > c1)
        slower += int(stat > c1 and s_hat < 1.0)
        found += int(stat_d > c1)
        s_all.append(s_hat)
    return dict(case=name, analysis="plain+score", score="hidden risk itself, no error", n_trials=nsim,
                slope_change_flagged=flag / nsim, flagged_as_slower_aging=slower / nsim,
                benefit_differs_by_score_found=found / nsim, median_slope_estimate=float(np.median(s_all)))


EXTRA = ("sturdy helped more (70% and 30% cuts)", "no effect", "aging slowed 20% in every dog")


def main():
    quick = "--quick" in sys.argv
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    same = check_against_experiment_001(a, b)
    assert same < 1e-6, same
    # (chunks, trials per chunk): the five-score case is split finer so two cores stay busy
    plan = {name: ((1, 8) if quick else ((4, 150) if len(c[3]) > 1 else (2, 300))) for name, c in CASES.items()}
    jobs = [(name, plan[name][1], a, b, SEED + 1000 * i + k)
            for i, name in enumerate(CASES) for k in range(plan[name][0])]
    with Pool(2) as pool:
        parts = pool.map(one_chunk, jobs, chunksize=1)
        extra = pool.map(one_extra, [(name, 8 if quick else 400, a, b, SEED + 90000 + i)
                                     for i, name in enumerate(EXTRA)], chunksize=1)
    rows = []
    for name, (r_sturdy, r_frail, s, rhos) in CASES.items():
        mine = [p for p in parts if p[0] == name]
        n = sum(p[1] for p in mine)
        for key in mine[0][3]:
            analysis, rho = key
            s_all = [x for p in mine for x in p[3][key]["s"]]
            row = dict(case=name, analysis=analysis, score_quality=rho, n_trials=n,
                       risk_multiplier_sturdier_half=r_sturdy, risk_multiplier_frailer_half=r_frail,
                       rate_multiplier=s, share_of_dogs_dead=sum(p[2] for p in mine) / n,
                       slope_change_flagged=sum(p[3][key]["flag"] for p in mine) / n,
                       flagged_as_slower_aging=sum(p[3][key]["slower"] for p in mine) / n,
                       median_slope_estimate=float(np.median(s_all)))
            if analysis == "plain+score":
                row["benefit_differs_by_score_found"] = sum(p[4][rho]["found"] for p in mine) / n
                row["median_benefit_by_score_estimate"] = float(np.median([x for p in mine for x in p[4][rho]["d"]]))
                row["median_score_effect_estimate"] = float(np.median([x for p in mine for x in p[4][rho]["g"]]))
            rows.append(row)
    out = dict(note="Simulated 10-year trials of 1,300 dogs, half treated from enrolment. 'slope_change_flagged' is "
                    "the share of trials in which a 5% likelihood-ratio test rejects 'the slope did not change'; "
                    "'flagged_as_slower_aging' is the share flagged with an estimated slope below the untreated "
                    "one. 'plain' ignores hidden differences; 'new' adds one setting for their spread; '+score' "
                    "adds the score's effect on risk and a 'benefit differs by score' term. ASSUMPTIONS: spread "
                    "0.3, entry ages 7 to 12, the form and quality of the score, and every size of benefit.",
               seed=SEED, uk_a_per_year=a, uk_b_per_year=b, spread=V, n_dogs=N, entry_ages=[LO, HI],
               years=FOLLOW, gradient_check_max_relative_error=grad_err,
               difference_from_experiment_001_likelihood=same, cases=rows,
               extra_score_is_hidden_risk_itself=extra)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case'][:38]:38s} | {r['analysis']:11s} | rho {str(r['score_quality']):4s}"
              f" | flagged {r['slope_change_flagged']:.3f} (slower {r['flagged_as_slower_aging']:.3f})"
              f" | s-hat {r['median_slope_estimate']:.3f}"
              + (f" | d found {r['benefit_differs_by_score_found']:.2f} d-hat {r['median_benefit_by_score_estimate']:+.2f}"
                 f" g-hat {r['median_score_effect_estimate']:+.2f}" if "benefit_differs_by_score_found" in r else ""))


    for r in extra:
        print(f"extra, score = hidden risk itself | {r['case'][:38]:38s} | flagged {r['slope_change_flagged']:.3f}"
              f" (slower {r['flagged_as_slower_aging']:.3f}) | s-hat {r['median_slope_estimate']:.3f}"
              f" | d found {r['benefit_differs_by_score_found']:.2f}")


if __name__ == "__main__":
    main()
