"""Experiment 013: if breeds differ in their rate of aging, does the slope test still work, and what repairs it?

Model for one dog of breed k, enrolled at age x0 and followed for t years:
    hazard(t) = a_k * exp(b_k * x0) * r * exp(b_k * s * t)
with r = s = 1 in the placebo arm (r < 1: risk cut; s < 1: the rise of risk slowed), as in experiment 001.

Inputs
  a, b      Gompertz fit to the UK dog life table (dog_lifespan/results.json). b is the "UK rate".
  breeds    Life expectancy at birth and number of dogs for 19 groups, Teng et al. 2022
            (dog_lifespan/data/teng2022_breed_e0.csv).
  a_k, b_k  ASSUMPTION. Every breed follows the doubling law and matches its published life expectancy e0_k.
            A share f of the difference between breeds is put on the rate:
              a_ref    solves  e0(a_ref, b) = count-weighted mean of e0_k
              b_k_full solves  e0(a_ref, b_k_full) = e0_k
              b_k = b * (b_k_full / b) ** f,   a_k solves e0(a_k, b_k) = e0_k
            where e0(a, b) = integral from 0 of exp(-(a / b) * (exp(b * x) - 1)) dx.
            f = 1 puts the whole difference on the rate (the extreme case); f = 0.5 half of it.
            No breed-by-breed rates were found in a source that could be opened, so f is an ASSUMPTION.
  mix       ASSUMPTION: the breed mix at birth follows the table's counts. A dog enters a trial only if alive
            at its entry age, so breed k is drawn with weight p_k * exp(-H_k(x0)).
  designs   Entry-age ranges are ASSUMPTIONS, as in experiments 001, 007 and 008.

Each simulated trial is analysed three ways:
  "ignoring breed"            the four models of experiment 001
  "level per breed"           one level per breed, one shared rate (experiment 008's breed-aware analysis)
  "level and rate per breed"  every breed expected to supply at least 3 percent of the dogs at entry gets its own
                              level and rate; rarer breeds are pooled into one group with its own level and rate.
                              The treatment multiplies the level by r and the rate by s, shared across groups.

Run:  python3 breed_slopes.py            (about 4 minutes on 2 cores)
      python3 breed_slopes.py --quick
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
from scipy.optimize import brentq, minimize
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "001_trial_fingerprint", "fingerprint.py")
bl = _load("breed_levels", "008_breed_levels", "breed_levels.py")

SEED = 20261118
MODELS = fp.MODELS
MIN_SHARE = 0.03                                   # a breed below this share of dogs at entry is pooled
DESIGNS = {"1,300 dogs, 10 years": (1300, 7.0, 12.0, 10.0), "1,300 dogs, 4 years": (1300, 10.0, 14.0, 4.0)}
EFFECTS = {"no effect": (1.0, 1.0), "risk cut 50%": (0.5, 1.0), "aging slowed 20%": (1.0, 0.8)}
DROPPED = {"without the French Bulldog": ("French Bulldog",), "all 19 groups": ()}
# (breed list, share f of the breed difference put on the rate, design)
WORLDS = [("without the French Bulldog", 1.0, "1,300 dogs, 10 years"),
          ("without the French Bulldog", 0.5, "1,300 dogs, 10 years"),
          ("all 19 groups", 1.0, "1,300 dogs, 10 years"),
          ("without the French Bulldog", 1.0, "1,300 dogs, 4 years")]


def e0(a, b):
    """Life expectancy at birth under the doubling law."""
    return quad(lambda x: np.exp(-(a / b) * np.expm1(b * x)), 0, 60, limit=200)[0]


def scenario(list_name, f, a, b):
    """Levels a_k, rates b_k and birth shares p_k for one breed list and one share f on the rate."""
    keep = [x for x in bl.breed_table() if x[0] not in DROPPED[list_name]]
    e = np.array([x[1] for x in keep])
    p = np.array([x[2] for x in keep])
    p = p / p.sum()
    a_ref = brentq(lambda la: e0(np.exp(la), b) - np.sum(p * e), np.log(a) - 4, np.log(a) + 4)
    a_ref = float(np.exp(a_ref))
    b_full = np.array([brentq(lambda bb: e0(a_ref, bb) - ek, 0.05, 5.0) for ek in e])
    b_k = b * (b_full / b) ** f
    a_k = np.array([np.exp(brentq(lambda la: e0(np.exp(la), bk) - ek, np.log(a) - 12, np.log(a) + 6))
                    for bk, ek in zip(b_k, e)])
    return dict(names=[x[0] for x in keep], e0=e, p=p, a_ref=a_ref, a_k=a_k, b_k=b_k)


def entry_weights(x0, a_k, b_k, p):
    """Breed mix among dogs alive at each entry age: p_k * exp(-H_k(x0)), normalised. Shape (n, K)."""
    H = a_k * np.expm1(np.outer(np.asarray(x0), b_k)) / b_k
    w = p * np.exp(-H)
    return w / w.sum(axis=1, keepdims=True)


def groups_for(sc, lo, hi):
    """Analysis groups: a breed with at least MIN_SHARE of the dogs at entry stands alone; the rest are pooled.
    Returns (group index for each breed, number of groups, share of each breed at entry)."""
    w = entry_weights(np.linspace(lo, hi, 2001), sc["a_k"], sc["b_k"], sc["p"]).mean(axis=0)
    own = w >= MIN_SHARE
    g = np.full(len(w), int(own.sum()))
    g[own] = np.arange(own.sum())
    return g, int(own.sum()) + int((~own).any()), w


def simulate(n, lo, hi, follow, r, s, sc, rng):
    x0 = rng.uniform(lo, hi, n)
    cdf = np.cumsum(entry_weights(x0, sc["a_k"], sc["b_k"], sc["p"]), axis=1)
    breed = np.minimum((rng.random(n)[:, None] > cdf).sum(axis=1), len(sc["p"]) - 1)
    arm = (np.arange(n) % 2).astype(float)
    bk = sc["b_k"][breed]
    h0 = sc["a_k"][breed] * np.exp(bk * x0) * np.where(arm == 1, r, 1.0)
    bs = bk * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died, breed


def strat_negll(theta, x0, arm, t, died, grp, G, D, has_r, has_s):
    """Level and rate per group. theta = [ln b for each group (G values), (ln r), (ln s)].
    For fixed rates, r and s, the best level for group g is (deaths in g) / (exposure sum in g), so the levels
    are solved exactly and only the rates, r and s are searched. By the envelope theorem the gradient is the
    partial gradient at those levels."""
    b = np.exp(theta[:G])[grp]
    k = G
    lnr = theta[k] if has_r else 0.0
    k += int(has_r)
    lns = theta[k] if has_s else 0.0
    s = np.exp(lns * arm)
    bs = b * s
    E = np.exp(bs * t)
    core = np.exp(b * x0 + lnr * arm)
    expo = core * (E - 1.0) / bs
    S = np.bincount(grp, weights=expo, minlength=G)
    lna = np.log(D / S)
    eA = np.exp(lna[grp]) * core
    H = eA * (E - 1.0) / bs
    ll = np.sum(died * (lna[grp] + b * x0 + lnr * arm + bs * t) - H)
    gb = b * (died * (x0 + s * t) - H * x0 - eA * t * E / b + H / b)
    grad = [np.bincount(grp, weights=gb, minlength=G)]
    if has_r:
        grad.append([np.sum(arm * (died - H))])
    if has_s:
        grad.append([np.sum(arm * (died * bs * t - eA * t * E + H))])
    return -ll, -np.concatenate(grad)


def strat_fit(model, data, grp, G, D, b, start):
    has_r, has_s = MODELS[model]
    bounds = ([(np.log(b) - 1.5, np.log(b) + 1.5)] * G
              + ([(np.log(0.05), np.log(5.0))] if has_r else []) + ([(np.log(0.05), np.log(3.0))] if has_s else []))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        res = minimize(strat_negll, start, args=(*data, grp, G, D, has_r, has_s), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-7})
    return -res.fun, res.x


def strat_fit_all(data, grp, G, b):
    """The four nested models. Each larger model is started from the smaller ones, so its likelihood is never
    below a model nested inside it. Returns None if a group has no deaths (its level cannot be estimated)."""
    D = np.bincount(grp, weights=data[3], minlength=G)
    if np.any(D < 2):
        return None
    lb = [np.log(b)] * G
    ll = {"null": strat_fit("null", data, grp, G, D, b, lb)}
    base = list(ll["null"][1])
    ll["level"] = strat_fit("level", data, grp, G, D, b, base + [0.0])
    ll["rate"] = strat_fit("rate", data, grp, G, D, b, base + [0.0])
    f1 = strat_fit("full", data, grp, G, D, b, list(ll["level"][1]) + [0.0])
    f2 = strat_fit("full", data, grp, G, D, b, list(ll["rate"][1][:G]) + [0.0, ll["rate"][1][G]])
    ll["full"] = f1 if f1[0] >= f2[0] else f2
    return ll


def check_gradient(a, b):
    sc = scenario("without the French Bulldog", 1.0, a, b)
    g, G, _ = groups_for(sc, 7.0, 12.0)
    rng = np.random.default_rng(5)
    *data, breed = simulate(1300, 7, 12, 6, 0.8, 0.9, sc, rng)
    grp = g[breed]
    D = np.bincount(grp, weights=data[3], minlength=G)
    theta = np.concatenate([np.log(b) + rng.normal(0, 0.05, G), [-0.2, -0.1]])
    _, gr = strat_negll(theta, *data, grp, G, D, True, True)
    worst = 0.0
    for i in range(theta.size):
        d = np.zeros(theta.size)
        d[i] = 1e-6
        num = (strat_negll(theta + d, *data, grp, G, D, True, True)[0]
               - strat_negll(theta - d, *data, grp, G, D, True, True)[0]) / 2e-6
        worst = max(worst, abs(gr[i] - num) / (abs(num) + 1e-6))
    return float(worst)


ANALYSES = ("ignoring breed", "level per breed", "level and rate per breed")


def one_case(args):
    list_name, f, design, effect, nsim, a, b, seed = args
    n, lo, hi, follow = DESIGNS[design]
    r, s = EFFECTS[effect]
    sc = scenario(list_name, f, a, b)
    K = len(sc["p"])
    g, G, _ = groups_for(sc, lo, hi)
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    out = {k: dict(any=0, flag=0, n=0, s=[]) for k in ANALYSES}
    deaths = []
    for _ in range(nsim):
        *data, breed = simulate(n, lo, hi, follow, r, s, sc, rng)
        deaths.append(float(data[3].sum()))
        grp = g[breed]
        fits = (("ignoring breed", fp.fit_all(tuple(data), a, b), 3),
                ("level per breed", bl.fit_all(data, breed, K, b), 2),
                ("level and rate per breed", strat_fit_all(data, grp, G, b), G + 1))
        for label, ll, i_s in fits:
            if ll is None:
                continue
            o = out[label]
            o["n"] += 1
            o["any"] += 2 * (ll["level"][0] - ll["null"][0]) > c1
            o["flag"] += 2 * (ll["full"][0] - ll["level"][0]) > c1
            o["s"].append(float(np.exp(ll["full"][1][i_s])))
    return [dict(breed_list=list_name, share_on_rate=f, design=design, effect=effect, analysis=label,
                 n_trials=nsim, trials_analysed=v["n"], share_of_dogs_dead=float(np.mean(deaths)) / n,
                 detected_by_usual_test=v["any"] / v["n"], slope_change_flagged=v["flag"] / v["n"],
                 median_slope_estimate=float(np.median(v["s"]))) for label, v in out.items()]


def describe(list_name, f, design, a, b):
    sc = scenario(list_name, f, a, b)
    _, lo, hi, _ = DESIGNS[design]
    g, G, w = groups_for(sc, lo, hi)
    own = [n for n, gi in zip(sc["names"], g) if gi < G - 1 or G == len(g)]
    return dict(breed_list=list_name, share_on_rate=f, design=design, reference_level=sc["a_ref"],
                analysis_groups=G, breeds_with_own_group=own,
                lowest_rate=float(sc["b_k"].min()), highest_rate=float(sc["b_k"].max()),
                rate_mean_at_entry=float(np.sum(w * sc["b_k"])),
                rate_sd_at_entry=float(np.sqrt(np.sum(w * (sc["b_k"] - np.sum(w * sc["b_k"])) ** 2))),
                breeds=[dict(group=n, life_expectancy=float(e), level=float(ak), rate=float(bk),
                             share_at_birth=float(p), share_at_entry=float(x))
                        for n, e, ak, bk, p, x in zip(sc["names"], sc["e0"], sc["a_k"], sc["b_k"], sc["p"], w)])


def main():
    quick = "--quick" in sys.argv
    nsim = 12 if quick else 1000
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    cases = [(ln, f, d, e, nsim, a, b, SEED + 10000 * i + k)
             for i, (ln, f, d) in enumerate(WORLDS) for k, e in enumerate(EFFECTS)]
    if quick:
        cases = cases[:2]
    with Pool(2) as pool:
        rows = [r for part in pool.map(one_case, cases, chunksize=1) for r in part]
    out = dict(note="Simulated trials. ASSUMPTIONS: every breed follows the doubling law and matches its published "
                    "life expectancy (Teng et al. 2022); a share f of the breed difference is put on the rate of "
                    "aging (f is an assumption; 1 is the extreme case); breed mix at birth follows the table's "
                    "counts; entry ages. 'slope_change_flagged' is the share of trials in which a 5% "
                    "likelihood-ratio test rejects 'the slope did not change'.",
               seed=SEED, dog_a_per_year=a, dog_b_per_year=b, trials_per_case=nsim, min_share_for_own_group=MIN_SHARE,
               gradient_check_max_relative_error=grad_err,
               worlds=[describe(ln, f, d, a, b) for ln, f, d in WORLDS], cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for w in out["worlds"]:
        print(w["breed_list"], w["share_on_rate"], w["design"], "groups", w["analysis_groups"],
              "rates %.3f-%.3f sd at entry %.3f" % (w["lowest_rate"], w["highest_rate"], w["rate_sd_at_entry"]))
    for r in rows:
        print(f"{r['breed_list'][:20]:20s} f={r['share_on_rate']} | {r['design'][12:]:8s} | {r['effect']:16s} | "
              f"{r['analysis']:24s} | dead {r['share_of_dogs_dead']:.2f} | usual {r['detected_by_usual_test']:.2f} | "
              f"slope flagged {r['slope_change_flagged']:.3f} | s-hat {r['median_slope_estimate']:.3f} | "
              f"n {r['trials_analysed']}")


if __name__ == "__main__":
    main()
