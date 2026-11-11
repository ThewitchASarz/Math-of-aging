"""Experiment 008: do differences in risk between breeds break the slope test, and does a baseline per breed fix it?

Model for one dog of breed k, enrolled at age x0 and followed for t years:
    hazard(t) = a_k * exp(b * x0) * r * exp(b * s * t)
with r = s = 1 in the placebo arm (r < 1: risk cut; s < 1: aging slowed), as in experiment 001.

Inputs
  b, a      Gompertz fit to the UK dog life table (dog_lifespan/results.json).
  breeds    Life expectancy at birth and number of dogs for 19 groups, Teng et al. 2022
            (dog_lifespan/data/teng2022_breed_e0.csv).
  a_k       ASSUMPTION: breeds differ only in level. a_k = a * z_k / mean(z), z_k = exp(-b * (e0_k - mean e0)),
            the conversion used in experiment 007. All breeds share the rate of aging b.
  mix       ASSUMPTION: the breed mix at birth follows the table's counts. A dog enters a trial only if alive
            at its entry age, so breed k is drawn with weight p_k * exp(-H_k(x0)).
  designs   Entry-age ranges are ASSUMPTIONS, as in experiments 001 and 007.

Each simulated trial is analysed twice: ignoring breed (the four models of experiment 001) and breed-aware
(one level of risk per breed, one shared b).

Run:  python3 breed_levels.py            (about 3 minutes on 2 cores)
      python3 breed_levels.py --quick
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
ROOT = os.path.join(HERE, "..", "..")
_spec = importlib.util.spec_from_file_location(
    "fingerprint", os.path.join(HERE, "..", "001_trial_fingerprint", "fingerprint.py"))
fp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fp)

SEED = 20261116
MODELS = fp.MODELS
DESIGNS = {"1,300 dogs, 4 years": (1300, 10.0, 14.0, 4.0), "1,300 dogs, 10 years": (1300, 7.0, 12.0, 10.0)}
EFFECTS = {"no effect": (1.0, 1.0), "risk cut 50%": (0.5, 1.0), "aging slowed 20%": (1.0, 0.8)}
# The three versions of the breed list measured in experiment 007.
DROPPED = {
    "all 19 groups": (),
    "without the French Bulldog": ("French Bulldog",),
    "without the four shortest-lived breeds": ("French Bulldog", "English Bulldog", "Pug", "American Bulldog"),
}


def breed_table():
    path = os.path.join(ROOT, "dog_lifespan", "data", "teng2022_breed_e0.csv")
    rows = [l.strip().split(",") for l in open(path) if not l.startswith("#")]
    cols, rows = rows[0], rows[1:]
    return [(r[cols.index("breed")], float(r[cols.index("e0")]), float(r[cols.index("n")])) for r in rows]


def scenario(name, a, b):
    """Breed levels a_k and birth shares p_k for one version of the breed list."""
    keep = [x for x in breed_table() if x[0] not in DROPPED[name]]
    e0 = np.array([x[1] for x in keep])
    p = np.array([x[2] for x in keep])
    p = p / p.sum()
    z = np.exp(-b * (e0 - np.sum(p * e0)))
    z = z / np.sum(p * z)                          # mean multiplier 1 at birth
    spread_birth = float(np.sum(p * (z - 1.0) ** 2))
    return dict(names=[x[0] for x in keep], a_k=a * z, p=p, spread_at_birth=spread_birth,
                highest_to_lowest_risk=float(z.max() / z.min()))


def entry_weights(x0, a_k, p, b):
    """Breed mix among dogs alive at each entry age: p_k * exp(-H_k(x0)), normalised. Shape (n, K)."""
    H = np.outer(np.expm1(b * np.asarray(x0)) / b, a_k)
    w = p * np.exp(-H)
    return w / w.sum(axis=1, keepdims=True)


def spread_at_entry(lo, hi, a_k, p, b):
    """Spread of the risk multipliers among dogs alive at entry, entry age uniform on [lo, hi]."""
    w = entry_weights(np.linspace(lo, hi, 2001), a_k, p, b).mean(axis=0)
    mean = np.sum(w * a_k)
    return float(np.sum(w * (a_k / mean - 1.0) ** 2))


def simulate(n, lo, hi, follow, r, s, a_k, p, b, rng):
    x0 = rng.uniform(lo, hi, n)
    cdf = np.cumsum(entry_weights(x0, a_k, p, b), axis=1)
    breed = np.minimum((rng.random(n)[:, None] > cdf).sum(axis=1), len(p) - 1)
    arm = (np.arange(n) % 2).astype(float)
    h0 = a_k[breed] * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died, breed


def negll_and_grad(theta, x0, arm, t, died, breed, K, has_r, has_s):
    """Breed-aware model. theta = [ln a for each breed (K values), ln b, (ln r), (ln s)]."""
    lna = theta[:K][breed]
    b = np.exp(theta[K])
    k = K + 1
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
    grad = [np.bincount(breed, weights=ga, minlength=K),
            [b * np.sum(died * (x0 + s * t) - H * x0 - eA * t * E / b + H / b)]]
    if has_r:
        grad.append([np.sum(arm * ga)])
    if has_s:
        grad.append([np.sum(arm * (died * bs * t - eA * t * E + H))])
    return -ll, -np.concatenate(grad)


def profiled(theta, x0, arm, t, died, breed, K, D, has_r, has_s):
    """Breed-aware likelihood with the K breed levels solved exactly for the given b, r, s.
    For fixed b, r, s the best level for breed k is (deaths in k) / (exposure sum in k), so only b, r and s are
    searched. By the envelope theorem the gradient is the last part of the full gradient at those levels.
    A breed with no deaths in the trial gets level zero (ln a = -60 here) and adds nothing to the likelihood."""
    b = np.exp(theta[0])
    k = 1
    lnr = theta[k] if has_r else 0.0
    k += int(has_r)
    lns = theta[k] if has_s else 0.0
    bs = b * np.exp(lns * arm)
    G = np.bincount(breed, weights=np.exp(b * x0 + lnr * arm) * np.expm1(bs * t) / bs, minlength=K)
    lna = np.full(K, -60.0)
    ok = D > 0
    lna[ok] = np.log(D[ok] / G[ok])
    f, g = negll_and_grad(np.concatenate([lna, theta]), x0, arm, t, died, breed, K, has_r, has_s)
    return f, g[K:]


def fit(model, data, breed, K, b, extra=None):
    """Returns (log-likelihood, [ln b, (ln r), (ln s)]). Bounds on b, r, s as in experiment 001."""
    has_r, has_s = MODELS[model]
    n_extra = int(has_r) + int(has_s)
    D = np.bincount(breed, weights=data[3], minlength=K)
    start = [np.log(b)] + (list(extra) if extra is not None else [0.0] * n_extra)
    bounds = ([(np.log(b) - 1, np.log(b) + 1)]
              + ([(np.log(0.05), np.log(5.0))] if has_r else []) + ([(np.log(0.05), np.log(3.0))] if has_s else []))
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(profiled, start, args=(*data, breed, K, D, has_r, has_s), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
    return -res.fun, res.x


def fit_all(data, breed, K, b):
    ll = {m: fit(m, data, breed, K, b) for m in ("null", "level", "rate")}
    f1 = fit("full", data, breed, K, b, extra=[ll["level"][1][1], 0.0])
    f2 = fit("full", data, breed, K, b, extra=[0.0, ll["rate"][1][1]])
    ll["full"] = f1 if f1[0] >= f2[0] else f2
    return ll


def check_gradient(a, b):
    sc = scenario("all 19 groups", a, b)
    K = len(sc["p"])
    rng = np.random.default_rng(5)
    *data, breed = simulate(800, 7, 12, 4, 0.8, 0.9, sc["a_k"], sc["p"], b, rng)
    theta = np.concatenate([np.log(sc["a_k"]) + 0.1, [np.log(b) - 0.05, -0.2, -0.1]])
    _, g = negll_and_grad(theta, *data, breed, K, True, True)
    worst = 0.0
    for i in range(theta.size):
        d = np.zeros(theta.size)
        d[i] = 1e-6
        num = (negll_and_grad(theta + d, *data, breed, K, True, True)[0]
               - negll_and_grad(theta - d, *data, breed, K, True, True)[0]) / 2e-6
        if abs(num) > 1e-6:                        # a breed with no dogs in this sample has gradient 0
            worst = max(worst, abs(g[i] - num) / abs(num))
    D = np.bincount(breed, weights=data[3], minlength=K)
    small = theta[K:]
    _, g = profiled(small, *data, breed, K, D, True, True)
    for i in range(3):
        d = np.zeros(3)
        d[i] = 1e-6
        num = (profiled(small + d, *data, breed, K, D, True, True)[0]
               - profiled(small - d, *data, breed, K, D, True, True)[0]) / 2e-6
        worst = max(worst, abs(g[i] - num) / abs(num))
    return float(worst)


def one_case(args):
    sc_name, design, effect, nsim, a, b, seed = args
    n, lo, hi, follow = DESIGNS[design]
    r, s = EFFECTS[effect]
    sc = scenario(sc_name, a, b)
    K = len(sc["p"])
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    out = {"ignoring breed": dict(any=0, flag=0, s=[]), "breed-aware": dict(any=0, flag=0, s=[])}
    deaths = []
    for _ in range(nsim):
        *data, breed = simulate(n, lo, hi, follow, r, s, sc["a_k"], sc["p"], b, rng)
        deaths.append(float(data[3].sum()))
        plain = fp.fit_all(tuple(data), a, b)
        aware = fit_all(data, breed, K, b)
        for label, ll, i_s in (("ignoring breed", plain, 3), ("breed-aware", aware, 2)):
            out[label]["any"] += 2 * (ll["level"][0] - ll["null"][0]) > c1
            out[label]["flag"] += 2 * (ll["full"][0] - ll["level"][0]) > c1
            out[label]["s"].append(float(np.exp(ll["full"][1][i_s])))
    return [dict(breed_list=sc_name, design=design, effect=effect, analysis=label, n_trials=nsim,
                 share_of_dogs_dead=float(np.mean(deaths)) / n,
                 detected_by_usual_test=v["any"] / nsim, slope_change_flagged=v["flag"] / nsim,
                 median_slope_estimate=float(np.median(v["s"]))) for label, v in out.items()]


def gamma_control(args):
    """Added after the main result was seen; not one of the predictions. Hidden differences drawn from a gamma
    distribution with spread 0.3 (the case that broke the test in experiment 004's check), simulated with
    experiment 004's own function, 1,300 dogs for 10 years, true effect a 50% risk cut, analysed ignoring them."""
    nsim, a, b, seed = args
    spec = importlib.util.spec_from_file_location("old_dogs", os.path.join(HERE, "..", "004_old_dogs", "old_dogs.py"))
    od = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(od)
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    flag = 0
    for _ in range(nsim):
        data = od.simulate_frail(1300, 7.0, 12.0, 10.0, 0.5, 1.0, a, b, 0.3, rng)
        ll = fp.fit_all(tuple(data[:4]), a, b)
        flag += 2 * (ll["full"][0] - ll["level"][0]) > c1
    return int(flag)


def entry_mix(a, b):
    """Risk multiplier and share of each group at birth and among dogs entering the 10-year design."""
    sc = scenario("all 19 groups", a, b)
    _, lo, hi, _ = DESIGNS["1,300 dogs, 10 years"]
    w = entry_weights(np.linspace(lo, hi, 2001), sc["a_k"], sc["p"], b).mean(axis=0)
    return [dict(group=n, risk_multiplier=float(z), share_at_birth=float(p), share_at_entry=float(x))
            for n, z, p, x in zip(sc["names"], sc["a_k"] / a, sc["p"], w)]


def main():
    quick = "--quick" in sys.argv
    nsim = 20 if quick else 1000
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    scen = {}
    for name in DROPPED:
        sc = scenario(name, a, b)
        scen[name] = dict(groups=len(sc["p"]), spread_at_birth=sc["spread_at_birth"],
                          highest_to_lowest_risk=sc["highest_to_lowest_risk"],
                          spread_among_dogs_at_entry={d: spread_at_entry(lo, hi, sc["a_k"], sc["p"], b)
                                                      for d, (_, lo, hi, _) in DESIGNS.items()})
    cases = [(sn, d, e, nsim, a, b, SEED + 10000 * i + 100 * j + k)
             for i, sn in enumerate(DROPPED) for j, d in enumerate(DESIGNS) for k, e in enumerate(EFFECTS)]
    if quick:
        cases = cases[:2]
    with Pool(2) as pool:
        rows = [r for part in pool.map(one_case, cases, chunksize=1) for r in part]
        half = nsim // 2
        control = sum(pool.map(gamma_control, [(half, a, b, SEED + 900000), (half, a, b, SEED + 900001)]))
    out = dict(note="Simulated trials. ASSUMPTIONS: breeds differ only in level of risk and share the UK rate of "
                    "aging; breed levels come from published life expectancies (Teng et al. 2022) by the "
                    "conversion of experiment 007; breed mix at birth follows the table's counts; entry ages. "
                    "'slope_change_flagged' is the share of trials in which a 5% likelihood-ratio test rejects "
                    "'the slope did not change'.",
               seed=SEED, dog_a_per_year=a, dog_b_per_year=b, trials_per_case=nsim,
               gradient_check_max_relative_error=grad_err, breed_lists=scen, cases=rows,
               added_after_the_result=dict(
                   note="Not predictions. A control with gamma-distributed hidden differences of spread 0.3 "
                        "(experiment 004's simulation), and the breed mix among dogs entering the 10-year design.",
                   gamma_spread_0_3_risk_cut_50_slope_change_flagged=control / (2 * half), n_trials=2 * half,
                   breed_mix_10_year_design=entry_mix(a, b)))
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k, v in scen.items():
        print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()})
    print("gamma control:", control / (2 * half))
    for r in rows:
        print(f"{r['breed_list'][:26]:26s} | {r['design']:20s} | {r['effect']:16s} | {r['analysis']:14s} | "
              f"usual {r['detected_by_usual_test']:.2f} | slope flagged {r['slope_change_flagged']:.3f} | "
              f"s-hat {r['median_slope_estimate']:.2f}")


if __name__ == "__main__":
    main()
