"""Experiment 007: are differences between breeds and sizes large enough to break the slope test?

Part 1  Spread of risk between breeds (UK) and size classes (US), from published life expectancies.
        Under the doubling law, a life expectancy that is D years shorter corresponds to a level of risk that
        is exp(b * D) times higher. The spread is reported as the variance of those multipliers divided by the
        square of their mean, comparable with the frailty variance of experiment 004.
Part 2  Trials in which each dog belongs to a size class with its own rate of aging (this repository's
        size-class calibration, "scenario A": larger dogs age faster). Class shares follow the US population.
        Each simulated trial is analysed twice: ignoring size, and with a separate baseline (a and b) for each
        of four strata (toy, small, medium, large-or-giant).

Run:  python3 breeds.py            (about 6 minutes on 2 cores)
      python3 breeds.py --quick
"""
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
ROOT = os.path.join(HERE, "..", "..")
SEED = 20261114
MODELS = {"null": (False, False), "level": (True, False), "rate": (False, True), "full": (True, True)}
DESIGNS = {"1,300 dogs, 4 years": (1300, 10.0, 14.0, 4.0), "1,300 dogs, 10 years": (1300, 7.0, 12.0, 10.0)}
EFFECTS = {"no effect": (1.0, 1.0), "risk cut 30%": (0.7, 1.0), "risk cut 50%": (0.5, 1.0),
           "aging slowed 20%": (1.0, 0.8)}
STRATUM = np.array([0, 1, 2, 3, 3])               # giant dogs are analysed together with large dogs


def read_csv(name):
    path = os.path.join(ROOT, "dog_lifespan", "data", name)
    rows = [l.strip().split(",") for l in open(path) if not l.startswith("#")]
    return rows[0], rows[1:]


def spread(e0, weights, b):
    z = np.exp(-b * (np.asarray(e0) - np.average(e0, weights=weights)))
    mean = np.average(z, weights=weights)
    return float(np.average((z / mean - 1.0) ** 2, weights=weights)), float(z.max() / z.min())


def part1():
    with open(os.path.join(ROOT, "dog_lifespan", "results.json")) as fh:
        fits = json.load(fh)["part1_mortality_fits"]
    b_uk = fits["RECOMMENDED_dog_parameters"]["b_per_year"]
    b_us = fits["RECOMMENDED_dog_parameters"]["alternative_US_model"]["b_per_year"]
    _, rows = read_csv("teng2022_breed_e0.csv")
    breeds = [(r[0], float(r[1]), float(r[4])) for r in rows]
    v_breed, ratio_breed = spread([b_[1] for b_ in breeds], [b_[2] for b_ in breeds], b_uk)
    pure = [b_ for b_ in breeds if b_[0].lower() != "crossbred"]
    v_pure, _ = spread([b_[1] for b_ in pure], [b_[2] for b_ in pure], b_uk)
    cols, rows = read_csv("montoya2023_size_e0.csv")
    e0 = [float(r[cols.index("e0")]) for r in rows]
    n = [float(r[cols.index("mean_annual_n")]) for r in rows]
    v_size, ratio_size = spread(e0, n, b_us)
    return dict(note="Variance of risk multipliers divided by the squared mean. Differences inside a breed or a "
                     "size class are not counted, so these are lower bounds on the true spread between dogs.",
                breeds_in_table=len(breeds), shortest_lived=min(breeds, key=lambda x: x[1])[:2],
                longest_lived=max(breeds, key=lambda x: x[1])[:2],
                variance_between_breeds_uk=v_breed, variance_between_pure_breeds_uk=v_pure,
                highest_to_lowest_breed_risk=ratio_breed, variance_between_size_classes_us=v_size,
                highest_to_lowest_size_risk=ratio_size, fitted_frailty_variance_experiment_004=0.102,
                variance_that_broke_the_test=0.3)


def size_classes():
    cols, rows = read_csv("size_class_gompertz_calibration.csv")
    b = np.array([float(r[cols.index("A_b")]) for r in rows])
    a = float(rows[0][cols.index("A_a")])
    cols2, rows2 = read_csv("montoya2023_size_e0.csv")
    n = np.array([float(r[cols2.index("mean_annual_n")]) for r in rows2])
    return a, b, n / n.sum(), [r[0] for r in rows]


def simulate(n, lo, hi, follow, r, s, a, b_class, shares, rng):
    cls = rng.choice(len(shares), size=n, p=shares)
    x0 = rng.uniform(lo, hi, n)
    arm = (np.arange(n) % 2).astype(float)
    b = b_class[cls]
    h0 = a * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died, cls


def negll_and_grad(theta, x0, arm, t, died, cls, K, has_r, has_s):
    """theta = [ln a for each stratum, ln b for each stratum, (ln r), (ln s)]."""
    lna = theta[:K][cls]
    b = np.exp(theta[K:2 * K])[cls]
    k = 2 * K
    lnr = theta[k] if has_r else 0.0
    k += int(has_r)
    lns = theta[k] if has_s else 0.0
    s = np.exp(lns * arm)
    A = lna + b * x0 + lnr * arm
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    H = eA * (E - 1.0) / bs
    ll = np.sum(died * (A + bs * t) - H)
    ga = died - H
    gb = b * (died * (x0 + s * t) - H * x0 - eA * t * E / b + H / b)
    grad = [np.bincount(cls, weights=ga, minlength=K), np.bincount(cls, weights=gb, minlength=K)]
    if has_r:
        grad.append([np.sum(arm * ga)])
    if has_s:
        grad.append([np.sum(arm * (died * bs * t - eA * t * E + H))])
    return -ll, -np.concatenate(grad)


def fit(model, data, cls, K, a0, b0, extra=None):
    has_r, has_s = MODELS[model]
    n_extra = int(has_r) + int(has_s)
    start = [np.log(a0)] * K + [np.log(b0)] * K + (list(extra) if extra is not None else [0.0] * n_extra)
    bounds = ([(np.log(a0) - 4, np.log(a0) + 4)] * K + [(np.log(b0) - 1, np.log(b0) + 1)] * K
              + ([(np.log(0.05), np.log(5.0))] if has_r else []) + ([(np.log(0.05), np.log(3.0))] if has_s else []))
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, start, args=(*data, cls, K, has_r, has_s), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-7})
    return -res.fun, res.x


def fit_all(data, cls, K, a0, b0):
    ll = {m: fit(m, data, cls, K, a0, b0) for m in ("null", "level", "rate")}
    f1 = fit("full", data, cls, K, a0, b0, extra=[ll["level"][1][2 * K], 0.0])
    f2 = fit("full", data, cls, K, a0, b0, extra=[0.0, ll["rate"][1][2 * K]])
    ll["full"] = f1 if f1[0] >= f2[0] else f2
    return ll


def check_gradient(a, b_class, shares):
    rng = np.random.default_rng(3)
    *data, cls5 = simulate(600, 7, 12, 4, 0.8, 0.9, a, b_class, shares, rng)
    cls, K = STRATUM[cls5], 4
    theta = np.concatenate([np.full(K, np.log(a) + 0.1), np.full(K, np.log(0.35)), [-0.2, -0.1]])
    _, g = negll_and_grad(theta, *data, cls, K, True, True)
    worst = 0.0
    for i in range(theta.size):
        d = np.zeros(theta.size)
        d[i] = 1e-6
        num = (negll_and_grad(theta + d, *data, cls, K, True, True)[0]
               - negll_and_grad(theta - d, *data, cls, K, True, True)[0]) / 2e-6
        worst = max(worst, abs(g[i] - num) / (abs(num) + 1e-9))
    return float(worst)


def one_case(args):
    design, effect, nsim, a, b_class, shares, seed = args
    n, lo, hi, follow = DESIGNS[design]
    r, s = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    b0 = float(np.average(b_class, weights=shares))
    out = {"ignoring size": dict(any=0, flag=0, s=[]), "size-aware": dict(any=0, flag=0, s=[])}
    for _ in range(nsim):
        *data, cls5 = simulate(n, lo, hi, follow, r, s, a, b_class, shares, rng)
        for label, cls, K in (("ignoring size", np.zeros(n, int), 1), ("size-aware", STRATUM[cls5], 4)):
            ll = fit_all(data, cls, K, a, b0)
            out[label]["any"] += 2 * (ll["level"][0] - ll["null"][0]) > c1
            out[label]["flag"] += 2 * (ll["full"][0] - ll["level"][0]) > c1
            out[label]["s"].append(float(np.exp(ll["full"][1][-1])))
    return [dict(design=design, effect=effect, analysis=label, n_trials=nsim,
                 detected_by_usual_test=v["any"] / nsim, slope_change_flagged=v["flag"] / nsim,
                 median_slope_estimate=float(np.median(v["s"]))) for label, v in out.items()]


def main():
    quick = "--quick" in sys.argv
    nsim = 20 if quick else 300
    a, b_class, shares, names = size_classes()
    grad_err = check_gradient(a, b_class, shares)
    assert grad_err < 1e-5, grad_err
    cases = [(d, e, nsim, a, b_class, shares, SEED + 1000 * i + j)
             for i, d in enumerate(DESIGNS) for j, e in enumerate(EFFECTS)]
    if quick:
        cases = cases[:2]
    with Pool(2) as pool:
        rows = [r for part in pool.map(one_case, cases) for r in part]
    out = dict(seed=SEED, part1_spread_of_risk=part1(),
               part2_trials_with_mixed_sizes=dict(
                   note="Size-class rates of aging are a calibration to published life expectancies (scenario A), "
                        "not a fit to full life tables. Class shares among trial dogs are assumed equal to US "
                        "population shares.",
                   size_classes=names, rate_of_aging_per_year=[float(x) for x in b_class],
                   population_share=[float(x) for x in shares], common_level_a=a, trials_per_case=nsim,
                   gradient_check_max_relative_error=grad_err, cases=rows))
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    p1 = out["part1_spread_of_risk"]
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in p1.items() if k != "note"})
    for r in rows:
        print(f"{r['design']:21s} | {r['effect']:17s} | {r['analysis']:13s} | usual {r['detected_by_usual_test']:.2f} | "
              f"slope flagged {r['slope_change_flagged']:.2f} | s-hat {r['median_slope_estimate']:.2f}")


if __name__ == "__main__":
    main()
