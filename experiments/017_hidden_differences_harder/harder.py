"""Experiment 017: the hidden-differences slope test of experiment 012 in harder conditions.

Model for one dog of size class c, entering at age x0, t years after enrolment:
    hazard(t) = z * a * exp(b_c * x0) * r * exp(b_c * s * t)
  z     hidden personal risk factor: gamma with mean 1 and variance V among dogs alive at age 2, so among dogs
        alive at entry age x0 it is gamma with shape 1/V and scale V / (1 + V * H0), where
        H0 = (a / b_c) * (exp(b_c * x0) - exp(2 * b_c)).   V = 0.3 is an ASSUMPTION (no measurement exists).
  r, s  treatment: r < 1 cuts risk, s < 1 slows aging; both 1 in the control group.

Part 1  4-year trial: entry ages 10 to 14 (ASSUMPTION), one size, no dropout, UK fit (dog_lifespan/results.json).
Part 2  10-year trial: entry ages 7 to 12 (ASSUMPTION), five size classes with their own rate of aging
        (dog_lifespan/data/size_class_gompertz_calibration.csv, scenario A, read through experiment 007),
        size mix = US population shares (ASSUMPTION), dropout at a constant 0.05 a year unrelated to health
        (ASSUMPTION), the same V inside every size class (ASSUMPTION).
Part 3  A treatment that helps frail dogs more. No dog's rate of aging changes. Treated dogs whose z is above
        the median for dogs of their entry age get r = 0.3, the others r = 0.7 (ASSUMPTION); a mirror case
        swaps the two. 10-year design of experiment 012, and the 4-year design of Part 1.

Two analyses of every simulated trial, each a 5% likelihood-ratio test of "risk cut only" against
"risk cut and slope change":
  plain  experiment 001's likelihood (experiment 007's, with four size groups, in Part 2)
  new    the same with one extra setting v for a gamma-shaped hidden risk factor (experiment 012's likelihood;
         in Part 2 the version below, with a level and a rate for each of four size groups).

Run:  python3 harder.py            (about 5 minutes on 2 cores)
      python3 harder.py --quick
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
from scipy.stats import gamma as gamma_dist

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "001_trial_fingerprint", "fingerprint.py")
br = _load("breeds", "007_breed_differences", "breeds.py")
hd = _load("hidden_differences", "012_hidden_differences_test", "hidden_differences.py")

SEED = 20261118
AGE0 = 2.0
N = 1300
V = 0.3                                   # ASSUMPTION: spread of hidden risk among dogs alive at age 2
DROP = 0.05                               # ASSUMPTION: dropout rate per year in Part 2
K = 4                                     # size groups in the analysis (giant analysed with large)
V_LO, V_HI = hd.V_LO, hd.V_HI
SHORT = (10.0, 14.0, 4.0)                 # entry ages: ASSUMPTION
LONG = (7.0, 12.0, 10.0)                  # entry ages: ASSUMPTION
# name: (part, design, spread, sizes and dropout?, r for sturdier half, r for frailer half, s)
CASES = {
    "4 years | no effect": (1, SHORT, V, False, 1.0, 1.0, 1.0),
    "4 years | risk cut 50%": (1, SHORT, V, False, 0.5, 0.5, 1.0),
    "4 years | aging slowed 20%": (1, SHORT, V, False, 1.0, 1.0, 0.8),
    "4 years, no hidden differences | aging slowed 20%": (1, SHORT, 0.0, False, 1.0, 1.0, 0.8),
    "10 years, sizes, dropout | no effect": (2, LONG, V, True, 1.0, 1.0, 1.0),
    "10 years, sizes, dropout | risk cut 50%": (2, LONG, V, True, 0.5, 0.5, 1.0),
    "10 years, sizes, dropout | aging slowed 20%": (2, LONG, V, True, 1.0, 1.0, 0.8),
    "10 years | frail helped more (70% and 30% cuts)": (3, LONG, V, False, 0.7, 0.3, 1.0),
    "10 years | sturdy helped more (30% and 70% cuts)": (3, LONG, V, False, 0.3, 0.7, 1.0),
    "4 years | frail helped more (70% and 30% cuts)": (3, SHORT, V, False, 0.7, 0.3, 1.0),
}


def simulate(design, v, mixed, r_sturdy, r_frail, s, a, b_one, b_class, shares, rng):
    lo, hi, follow = design
    x0 = rng.uniform(lo, hi, N)
    arm = (np.arange(N) % 2).astype(float)
    if mixed:
        cls = rng.choice(len(shares), size=N, p=shares)
        b = b_class[cls]
    else:
        cls = np.zeros(N, int)
        b = np.full(N, b_one)
    if v > 0:
        H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
        u = rng.gamma(1.0 / v, 1.0, N)                         # rank of the dog among dogs of its entry age
        z = u * v / (1.0 + v * H0)
        frail = u > gamma_dist.ppf(0.5, 1.0 / v)
    else:
        z = np.ones(N)
        frail = np.zeros(N, bool)
    r = np.where(arm == 1, np.where(frail, r_frail, r_sturdy), 1.0)
    h0 = z * a * np.exp(b * x0) * r
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(N)) * bs / h0) / bs
    leave = rng.exponential(1.0 / DROP, N) if mixed else np.full(N, np.inf)
    end = np.minimum(leave, follow)
    died = (t <= end).astype(float)
    return (x0, arm, np.minimum(t, end), died), cls


def negll_and_grad(theta, x0, arm, t, died, cls, has_s):
    """New analysis with size groups. theta = [ln a (K), ln b (K), ln v, ln r, (ln s)]."""
    lna = theta[:K][cls]
    b = np.exp(theta[K:2 * K])[cls]
    v = np.exp(theta[2 * K])
    lnr = theta[2 * K + 1]
    lns = theta[2 * K + 2] if has_s else 0.0
    a = np.exp(lna)
    s = np.exp(lns * arm)
    A = lna + b * x0 + lnr * arm
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    L = eA * (E - 1.0) / bs
    ebx, e2b = np.exp(b * x0), np.exp(b * AGE0)
    H0 = (a / b) * (ebx - e2b)
    D0 = 1.0 + v * H0
    D = D0 + v * L
    lnratio = np.log1p(v * L / D0)
    ll = np.sum(died * (A + bs * t - np.log(D)) - lnratio / v)
    c = (1.0 + v * died) / D
    dL_b = b * x0 * L + eA * t * E - L
    dH0_b = -H0 + a * (x0 * ebx - AGE0 * e2b)
    ga = died - c * (H0 + L) + H0 / D0
    gb = died * b * (x0 + s * t) - c * (dH0_b + dL_b) + dH0_b / D0
    grad = [np.bincount(cls, weights=ga, minlength=K), np.bincount(cls, weights=gb, minlength=K),
            [np.sum(-died * v * (H0 + L) / D + lnratio / v - (H0 + L) / D + H0 / D0)],
            [np.sum(arm * (died - c * L))]]
    if has_s:
        grad.append([np.sum(arm * (died * bs * t - c * (eA * t * E - L)))])
    return -ll, -np.concatenate(grad)


def fit_strat(has_s, data, cls, a0, b0, start):
    bounds = ([(np.log(a0) - 4, np.log(a0) + 4)] * K + [(np.log(b0) - 1, np.log(b0) + 1)] * K
              + [(np.log(V_LO), np.log(V_HI)), (np.log(0.05), np.log(5.0))]
              + ([(np.log(0.05), np.log(3.0))] if has_s else []))
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, start, args=(*data, cls, has_s), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-7})
    return -res.fun, res.x


def best(fits):
    return max(fits, key=lambda f: f[0])


def tests_one_size(data, a, b):
    """Plain and new slope tests with one baseline. Returns {analysis: (statistic, s-hat, v-hat)}."""
    p_lvl = fp.fit("level", data, a, b)
    p_full = best([fp.fit("full", data, a, b, start_extra=[p_lvl[1][2], 0.0]),
                   fp.fit("full", data, a, b, start_extra=[0.0, np.log(0.8)])])
    q, f = p_lvl[1], p_full[1]
    n_lvl = best([hd.fit("level", data, a, b, [q[0], q[1], np.log(v0), q[2]]) for v0 in (0.02, 0.5)])
    m = n_lvl[1]
    n_full = best([hd.fit("full", data, a, b, [m[0], m[1], m[2], m[3], 0.0]),
                   hd.fit("full", data, a, b, [f[0], f[1], np.log(0.02), f[2], f[3]]),
                   hd.fit("full", data, a, b, [f[0], f[1], np.log(0.5), f[2], f[3]])])
    return {"plain": (2 * (p_full[0] - p_lvl[0]), float(np.exp(f[3])), None),
            "new": (2 * (n_full[0] - n_lvl[0]), float(np.exp(n_full[1][4])), float(np.exp(m[2])))}


def tests_sizes(data, cls, a, b0):
    """The same two tests with a level and a rate for each of four size groups."""
    p_lvl = br.fit("level", data, cls, K, a, b0)
    p_full = best([br.fit("full", data, cls, K, a, b0, extra=[p_lvl[1][2 * K], 0.0]),
                   br.fit("full", data, cls, K, a, b0, extra=[0.0, np.log(0.8)])])
    q, f = p_lvl[1], p_full[1]
    n_lvl = best([fit_strat(False, data, cls, a, b0, [*q[:2 * K], np.log(v0), q[2 * K]]) for v0 in (0.02, 0.5)])
    m = n_lvl[1]
    n_full = best([fit_strat(True, data, cls, a, b0, [*m, 0.0]),
                   fit_strat(True, data, cls, a, b0, [*f[:2 * K], np.log(0.02), f[2 * K], f[2 * K + 1]]),
                   fit_strat(True, data, cls, a, b0, [*f[:2 * K], np.log(0.5), f[2 * K], f[2 * K + 1]])])
    return {"plain": (2 * (p_full[0] - p_lvl[0]), float(np.exp(f[-1])), None),
            "new": (2 * (n_full[0] - n_lvl[0]), float(np.exp(n_full[1][-1])), float(np.exp(m[2 * K])))}


def check_gradient(a, b_class, shares):
    rng = np.random.default_rng(5)
    data, cls5 = simulate(LONG, V, True, 0.8, 0.8, 0.9, a, 0.3, b_class, shares, rng)
    cls = br.STRATUM[cls5]
    theta = np.concatenate([np.full(K, np.log(a) + 0.1), np.full(K, np.log(0.35)), [np.log(0.2), -0.2, -0.1]])
    _, g = negll_and_grad(theta, *data, cls, True)
    worst = 0.0
    for i in range(theta.size):
        d = np.zeros(theta.size)
        d[i] = 1e-6
        num = (negll_and_grad(theta + d, *data, cls, True)[0] - negll_and_grad(theta - d, *data, cls, True)[0]) / 2e-6
        worst = max(worst, abs(g[i] - num) / (abs(num) + 1e-9))
    return float(worst)


def one_chunk(args):
    name, nsim, a_uk, b_uk, a_sz, b_class, shares, seed = args
    part, design, v, mixed, r_sturdy, r_frail, s = CASES[name]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    a = a_sz if mixed else a_uk
    b0 = float(np.average(b_class, weights=shares))
    out = {k: dict(flag=0, slower=0, s=[], v=[]) for k in ("plain", "new")}
    dead = left = 0.0
    for _ in range(nsim):
        data, cls5 = simulate(design, v, mixed, r_sturdy, r_frail, s, a, b_uk, b_class, shares, rng)
        dead += data[3].mean()
        left += np.mean((data[3] == 0) & (data[2] < design[2]))
        res = tests_sizes(data, br.STRATUM[cls5], a, b0) if mixed else tests_one_size(data, a, b_uk)
        for label, (stat, s_hat, v_hat) in res.items():
            out[label]["flag"] += int(stat > c1)
            out[label]["slower"] += int(stat > c1 and s_hat < 1.0)
            out[label]["s"].append(s_hat)
            if v_hat is not None:
                out[label]["v"].append(v_hat)
    return name, nsim, dead, left, out


def main():
    quick = "--quick" in sys.argv
    chunks, per = (1, 10) if quick else (2, 400)
    a_uk, b_uk = fp.dog_parameters()
    a_sz, b_class, shares, names = br.size_classes()
    grad_err = check_gradient(a_sz, b_class, shares)
    assert grad_err < 1e-5, grad_err
    jobs = [(name, per, a_uk, b_uk, a_sz, b_class, shares, SEED + 1000 * i + k)
            for i, name in enumerate(CASES) for k in range(chunks)]
    if quick:
        jobs = [j for j in jobs if j[0] in ("4 years | risk cut 50%", "10 years, sizes, dropout | risk cut 50%",
                                            "10 years | frail helped more (70% and 30% cuts)")]
    order = sorted(range(len(jobs)), key=lambda i: -CASES[jobs[i][0]][0] % 3)      # scheduling only: later parts first
    with Pool(2) as pool:
        parts = pool.map(one_chunk, [jobs[i] for i in order], chunksize=1)
    rows = []
    for name in CASES:
        mine = [p for p in parts if p[0] == name]
        if not mine:
            continue
        n = sum(p[1] for p in mine)
        part, design, v, mixed, r_sturdy, r_frail, s = CASES[name]
        for label in ("plain", "new"):
            s_all = [x for p in mine for x in p[4][label]["s"]]
            row = dict(case=name, part=part, analysis=label, n_trials=n, years=design[2],
                       entry_ages=list(design[:2]), spread=v, sizes_and_dropout=mixed,
                       risk_multiplier_sturdier_half=r_sturdy, risk_multiplier_frailer_half=r_frail,
                       rate_multiplier=s,
                       share_of_dogs_dead=sum(p[2] for p in mine) / len(mine) / mine[0][1],
                       share_of_dogs_dropped_out=sum(p[3] for p in mine) / len(mine) / mine[0][1],
                       slope_change_flagged=sum(p[4][label]["flag"] for p in mine) / n,
                       flagged_as_slower_aging=sum(p[4][label]["slower"] for p in mine) / n,
                       median_slope_estimate=float(np.median(s_all)))
            if label == "new":
                v_all = np.array([x for p in mine for x in p[4]["new"]["v"]])
                row["median_spread_estimate"] = float(np.median(v_all))
                row["spread_estimate_at_lower_bound"] = float(np.mean(v_all < 2 * V_LO))
            rows.append(row)
    out = dict(note="Simulated trials of 1,300 dogs, half treated from enrolment. 'slope_change_flagged' is the "
                    "share of trials in which a 5% likelihood-ratio test rejects 'the slope did not change'. "
                    "'plain' ignores hidden differences; 'new' adds one setting for the spread of a gamma-shaped "
                    "hidden risk factor. ASSUMPTIONS: spread 0.3, entry ages, dropout 0.05 a year, size mix, "
                    "and the 70% and 30% risk cuts of Part 3.",
               seed=SEED, uk_a_per_year=a_uk, uk_b_per_year=b_uk, size_class_a_per_year=a_sz,
               size_class_b_per_year=[float(x) for x in b_class], size_class_shares=[float(x) for x in shares],
               size_classes=names, dropout_per_year=DROP, trials_per_case=chunks * per,
               gradient_check_max_relative_error=grad_err, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case'][:50]:50s} | {r['analysis']:5s} | dead {r['share_of_dogs_dead']:.2f}"
              f" | flagged {r['slope_change_flagged']:.3f} (slower {r['flagged_as_slower_aging']:.3f})"
              f" | s-hat {r['median_slope_estimate']:.3f}"
              + (f" | v-hat {r['median_spread_estimate']:.2f} (at limit {r['spread_estimate_at_lower_bound']:.2f})"
                 if "median_spread_estimate" in r else ""))


if __name__ == "__main__":
    main()
