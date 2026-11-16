"""Experiment 012: a slope test that allows for hidden differences between dogs.

Each dog has a hidden personal risk factor z. For a dog enrolled at age x0 and followed for t years
    hazard(t) = z * a * exp(b * x0) * r * exp(b * s * t)
with r = s = 1 in the placebo arm (r < 1: risk cut; s < 1: aging slowed), as in experiment 001.

Worlds (how z is spread among dogs alive at age 2; mean 1, variance v)
  A  no hidden differences (v = 0)
  B  gamma, v = 0.1   about the size fitted to the UK life table in experiment 004 (0.102)
  C  gamma, v = 0.3   ASSUMPTION: the control that broke the plain test in experiment 008
  D  log-normal, v = 0.3   ASSUMPTION: a shape the analysis does not assume
A dog enters the trial only if alive at its entry age, so frail dogs are rarer among older entrants.

Two analyses of every simulated trial
  plain   the four models of experiment 001 (no hidden differences)
  new     the same models with one extra setting v: z is taken to be gamma with variance v at age 2.
          With H0 = (a/b)(exp(b*x0) - exp(2b)) the risk built up between age 2 and entry, and
          L(t) = a exp(b*x0) r (exp(b*s*t) - 1)/(b*s) the risk built up in the trial,
              chance of surviving to t = ((1 + v*H0 + v*L(t)) / (1 + v*H0)) ** (-1/v)
              hazard seen in the group = a exp(b*x0) r exp(b*s*t) / (1 + v*H0 + v*L(t))
          The term v*H0 allows for frail dogs already being rarer at entry.

Inputs: a, b from the UK fit in dog_lifespan/results.json. Design (1,300 dogs, entry ages 7 to 12, 10 years, no
dropout): entry ages are an ASSUMPTION, as in experiments 001, 010 and 011.

Run:  python3 hidden_differences.py            (about 5 minutes on 2 cores)
      python3 hidden_differences.py --quick
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
od = _load("old_dogs", "004_old_dogs", "old_dogs.py")

SEED = 20261120
AGE0 = od.AGE0                                    # hidden differences are defined among dogs alive at age 2
N, LO, HI, FOLLOW = 1300, 7.0, 12.0, 10.0         # entry ages: ASSUMPTION
WORLDS = {
    "A: no hidden differences": ("none", 0.0),
    "B: gamma, spread 0.1": ("gamma", 0.1),
    "C: gamma, spread 0.3": ("gamma", 0.3),        # ASSUMPTION
    "D: log-normal, spread 0.3": ("lognormal", 0.3),   # ASSUMPTION
}
EFFECTS = {"no effect": (1.0, 1.0), "risk cut 50%": (0.5, 1.0), "aging slowed 20%": (1.0, 0.8)}
MODELS = fp.MODELS
V_LO, V_HI = 1e-4, 5.0                            # bounds on the spread setting in the new analysis


def simulate_lognormal(n, lo, hi, follow, r, s, a, b, v, rng):
    """z is log-normal with mean 1 and variance v among dogs alive at age 2. A dog with factor z is still alive
    at entry age x0 with chance exp(-z * H0(x0)), so z is drawn by accept/reject with that chance."""
    x0 = rng.uniform(lo, hi, n)
    arm = (np.arange(n) % 2).astype(float)
    sig2 = np.log1p(v)
    H0 = od.cum_hazard(x0, a, b)
    z = np.empty(n)
    todo = np.arange(n)
    while todo.size:
        cand = rng.lognormal(-sig2 / 2, np.sqrt(sig2), todo.size)
        ok = rng.random(todo.size) < np.exp(-cand * H0[todo])
        z[todo[ok]] = cand[ok]
        todo = todo[~ok]
    h0 = z * a * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def simulate(world, r, s, a, b, rng):
    kind, v = WORLDS[world]
    if kind == "none":
        return fp.simulate(N, LO, HI, FOLLOW, r, s, a, b, rng)
    if kind == "gamma":
        return od.simulate_frail(N, LO, HI, FOLLOW, r, s, a, b, v, rng)
    return simulate_lognormal(N, LO, HI, FOLLOW, r, s, a, b, v, rng)


def negll_and_grad(theta, x0, arm, t, died, has_r, has_s):
    """New analysis. theta = [ln a, ln b, ln v, (ln r), (ln s)]. Returns minus log-likelihood and gradient."""
    lna, lnb, lnv = theta[0], theta[1], theta[2]
    k = 3
    lnr = theta[k] if has_r else 0.0
    k += int(has_r)
    lns = theta[k] if has_s else 0.0
    a, b, v = np.exp(lna), np.exp(lnb), np.exp(lnv)
    s = np.exp(lns * arm)
    A = lna + b * x0 + lnr * arm
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    L = eA * (E - 1.0) / bs                        # risk built up in the trial
    ebx, e2b = np.exp(b * x0), np.exp(b * AGE0)
    H0 = (a / b) * (ebx - e2b)                     # risk built up from age 2 to entry
    D0 = 1.0 + v * H0
    D = D0 + v * L
    lnratio = np.log1p(v * L / D0)                 # ln D - ln D0
    ll = np.sum(died * (A + bs * t - np.log(D)) - lnratio / v)
    c = (1.0 + v * died) / D
    dL_b = b * x0 * L + eA * t * E - L             # d L / d ln b
    dH0_b = -H0 + a * (x0 * ebx - AGE0 * e2b)      # d H0 / d ln b
    grad = [np.sum(died - c * (H0 + L) + H0 / D0),
            np.sum(died * b * (x0 + s * t) - c * (dH0_b + dL_b) + dH0_b / D0),
            np.sum(-died * v * (H0 + L) / D + lnratio / v - (H0 + L) / D + H0 / D0)]
    if has_r:
        grad.append(np.sum(arm * (died - c * L)))
    if has_s:
        grad.append(np.sum(arm * (died * bs * t - c * (eA * t * E - L))))
    return -ll, -np.array(grad)


def fit(model, data, a, b, start):
    has_r, has_s = MODELS[model]
    bounds = [(np.log(a) - 3, np.log(a) + 3), (np.log(b) - 1, np.log(b) + 1), (np.log(V_LO), np.log(V_HI))]
    if has_r:
        bounds.append((np.log(0.05), np.log(5.0)))
    if has_s:
        bounds.append((np.log(0.05), np.log(3.0)))
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, start, args=(*data, has_r, has_s), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-7})
    return -res.fun, res.x, bool(res.success)


def best(fits):
    return max(fits, key=lambda f: f[0])


def fit_all(data, a, b, plain):
    """The four nested models with the extra spread setting. Each model is started twice: from the plain fit of
    the same model with a small spread, and with a large spread; the better result is kept. The full model is
    also started from the level and rate solutions, so its likelihood is never below a model nested inside it."""
    ll = {}
    for m in ("null", "level", "rate"):
        p = plain[m][1]
        ll[m] = best([fit(m, data, a, b, [p[0], p[1], np.log(v0), *p[2:]]) for v0 in (0.02, 0.5)])
    base_l, base_r = ll["level"][1], ll["rate"][1]
    p = plain["full"][1]
    ll["full"] = best([fit("full", data, a, b, [*base_l[:4], 0.0]),
                       fit("full", data, a, b, [*base_r[:3], 0.0, base_r[3]]),
                       fit("full", data, a, b, [p[0], p[1], np.log(0.5), p[2], p[3]])])
    return ll


def check_gradient(a, b):
    rng = np.random.default_rng(1)
    data = od.simulate_frail(400, 7, 12, 6, 0.8, 0.9, a, b, 0.3, rng)
    theta = np.array([np.log(a) + 0.1, np.log(b) - 0.05, np.log(0.2), -0.2, -0.1])
    _, g = negll_and_grad(theta, *data, True, True)
    num = []
    for i in range(5):
        d = np.zeros(5)
        d[i] = 1e-6
        num.append((negll_and_grad(theta + d, *data, True, True)[0]
                    - negll_and_grad(theta - d, *data, True, True)[0]) / 2e-6)
    return float(np.max(np.abs(g - np.array(num)) / (np.abs(num) + 1e-9)))


def one_chunk(args):
    world, effect, nsim, a, b, seed = args
    r, s = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    out = {k: dict(any=0, flag=0, slower=0, s=[], v=[]) for k in ("plain", "new")}
    deaths = 0.0
    for _ in range(nsim):
        data = simulate(world, r, s, a, b, rng)
        deaths += data[3].sum()
        plain = fp.fit_all(data, a, b)
        new = fit_all(data, a, b, plain)
        for label, ll, i_s in (("plain", plain, 3), ("new", new, 4)):
            flag = 2 * (ll["full"][0] - ll["level"][0]) > c1
            out[label]["any"] += int(2 * (ll["level"][0] - ll["null"][0]) > c1)
            out[label]["flag"] += int(flag)
            out[label]["slower"] += int(flag and ll["full"][1][i_s] < 0)
            out[label]["s"].append(float(np.exp(ll["full"][1][i_s])))
        out["new"]["v"].append(float(np.exp(new["level"][1][2])))
    return world, effect, nsim, deaths, out


def main():
    quick = "--quick" in sys.argv
    chunks, per = (1, 12) if quick else (2, 400)
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    jobs = [(w, e, per, a, b, SEED + 10000 * i + 100 * j + k)
            for i, w in enumerate(WORLDS) for j, e in enumerate(EFFECTS) for k in range(chunks)]
    if quick:
        jobs = [j for j in jobs if j[0].startswith(("A", "D")) and j[1] != "no effect"]
    with Pool(2) as pool:
        parts = pool.map(one_chunk, jobs, chunksize=1)
    rows = []
    for w in WORLDS:
        for e in EFFECTS:
            mine = [p for p in parts if p[0] == w and p[1] == e]
            if not mine:
                continue
            n = sum(p[2] for p in mine)
            for label in ("plain", "new"):
                s_all = [x for p in mine for x in p[4][label]["s"]]
                row = dict(world=w, effect=e, analysis=label, n_trials=n,
                           share_of_dogs_dead=sum(p[3] for p in mine) / n / N,
                           detected_by_usual_test=sum(p[4][label]["any"] for p in mine) / n,
                           slope_change_flagged=sum(p[4][label]["flag"] for p in mine) / n,
                           flagged_as_slower_aging=sum(p[4][label]["slower"] for p in mine) / n,
                           median_slope_estimate=float(np.median(s_all)))
                if label == "new":
                    v_all = np.array([x for p in mine for x in p[4]["new"]["v"]])
                    row["median_spread_estimate"] = float(np.median(v_all))
                    row["spread_estimate_at_lower_bound"] = float(np.mean(v_all < 2 * V_LO))
                rows.append(row)
    out = dict(note="Simulated trials: 1,300 dogs, entry ages 7 to 12 (ASSUMPTION), half treated from enrolment, "
                    "10 years, no dropout. UK dog mortality fit. The spreads in worlds C and D are ASSUMPTIONS. "
                    "'slope_change_flagged' is the share of trials in which a 5% likelihood-ratio test rejects "
                    "'the slope did not change'. 'plain' ignores hidden differences; 'new' adds one setting for "
                    "the spread of a gamma-shaped hidden risk factor.",
               seed=SEED, dog_a_per_year=a, dog_b_per_year=b, trials_per_case=chunks * per,
               gradient_check_max_relative_error=grad_err, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['world'][:26]:26s} | {r['effect']:16s} | {r['analysis']:5s} | usual {r['detected_by_usual_test']:.2f}"
              f" | slope flagged {r['slope_change_flagged']:.3f} (slower {r['flagged_as_slower_aging']:.3f})"
              f" | s-hat {r['median_slope_estimate']:.2f}"
              + (f" | v-hat {r['median_spread_estimate']:.2f}" if "median_spread_estimate" in r else ""))


if __name__ == "__main__":
    main()
