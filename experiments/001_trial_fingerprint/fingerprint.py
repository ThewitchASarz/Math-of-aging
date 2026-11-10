"""Experiment 001: can a dog trial tell slowed aging (a slope change) from delayed death (a level shift)?

Model for one dog enrolled at age x0 and followed for t years:
    hazard(t) = a * exp(b * x0) * r * exp(b * s * t)
with r = s = 1 in the placebo arm. In the treated arm
    r < 1 : the risk of death is cut by a fixed fraction from enrolment (level effect)
    s < 1 : the yearly growth of risk is slowed from enrolment (rate effect)

a, b: Gompertz fit to the UK dog life table (dog_lifespan/results.json). Entry-age ranges are ASSUMPTIONS.
The effect is assumed to start at enrolment and persist for the whole follow-up.

Run:  python3 fingerprint.py            (about 5 minutes on 2 cores)
      python3 fingerprint.py --quick    (small run for tests)
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
SEED = 20261110

# name: (dogs, youngest entry age, oldest entry age, years of follow-up)
DESIGNS = {
    "580 dogs, 3 years (like the rapamycin trial)": (580, 7.0, 12.0, 3.0),
    "1,300 dogs, 4 years (like Loyal's trial)": (1300, 10.0, 14.0, 4.0),
    "580 dogs, 3 years, younger at entry": (580, 7.0, 9.0, 3.0),
    "5,000 dogs, 3 years (hypothetical)": (5000, 7.0, 12.0, 3.0),
    "1,300 dogs, 10 years (hypothetical)": (1300, 7.0, 12.0, 10.0),
}
# name: (kind, r, s)
EFFECTS = {
    "no effect": ("none", 1.0, 1.0),
    "risk cut 20%": ("level", 0.8, 1.0),
    "risk cut 30%": ("level", 0.7, 1.0),
    "risk cut 50%": ("level", 0.5, 1.0),
    "aging slowed 10%": ("rate", 1.0, 0.9),
    "aging slowed 20%": ("rate", 1.0, 0.8),
    "aging slowed 30%": ("rate", 1.0, 0.7),
}
MODELS = {"null": (False, False), "level": (True, False), "rate": (False, True), "full": (True, True)}


def dog_parameters():
    with open(os.path.join(ROOT, "dog_lifespan", "results.json")) as fh:
        d = json.load(fh)["part1_mortality_fits"]["RECOMMENDED_dog_parameters"]
    return d["a_per_year"], d["b_per_year"]


def simulate(n, lo, hi, follow, r, s, a, b, rng):
    """Entry ages uniform on [lo, hi]; alternate arms; exact death times by inverse transform."""
    x0 = rng.uniform(lo, hi, n)
    arm = (np.arange(n) % 2).astype(float)
    h0 = a * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def negll_and_grad(theta, x0, arm, t, died, has_r, has_s):
    """Negative log-likelihood and its gradient. theta = [ln a, ln b, (ln r), (ln s)]."""
    lna, lnb = theta[0], theta[1]
    k = 2
    lnr = theta[k] if has_r else 0.0
    k += int(has_r)
    lns = theta[k] if has_s else 0.0
    b = np.exp(lnb)
    s = np.exp(lns * arm)
    A = lna + b * x0 + lnr * arm
    bs = b * s
    E = np.exp(bs * t)
    eA = np.exp(A)
    H = eA * (E - 1.0) / bs                      # cumulative hazard over follow-up
    ll = np.sum(died * (A + bs * t) - H)
    grad = [np.sum(died - H),
            b * np.sum(died * (x0 + s * t) - H * x0 - eA * t * E / b + H / b)]
    if has_r:
        grad.append(np.sum(arm * (died - H)))
    if has_s:
        grad.append(np.sum(arm * (died * bs * t - eA * t * E + H)))
    return -ll, -np.array(grad)


def fit(model, data, a, b, start_extra=None):
    """Maximum likelihood with box bounds, so the optimiser cannot wander into overflow.
    Bounds: a within a factor 20, b within a factor 2.7, r in [0.05, 5], s in [0.05, 3]."""
    has_r, has_s = MODELS[model]
    start = [np.log(a), np.log(b)] + (list(start_extra) if start_extra is not None
                                      else [0.0] * (int(has_r) + int(has_s)))
    bounds = [(np.log(a) - 3, np.log(a) + 3), (np.log(b) - 1, np.log(b) + 1)]
    if has_r:
        bounds.append((np.log(0.05), np.log(5.0)))
    if has_s:
        bounds.append((np.log(0.05), np.log(3.0)))
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, start, args=(*data, has_r, has_s), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
    return -res.fun, res.x, bool(res.success)


def fit_all(data, a, b):
    """Fit the four nested models. The full model is started from each one-effect solution and the better
    result kept, so its likelihood is never below a model nested inside it."""
    ll = {m: fit(m, data, a, b) for m in ("null", "level", "rate")}
    f1 = fit("full", data, a, b, start_extra=[ll["level"][1][2], 0.0])
    f2 = fit("full", data, a, b, start_extra=[0.0, ll["rate"][1][2]])
    ll["full"] = f1 if f1[0] >= f2[0] else f2
    return ll


def one_case(args):
    design, effect, nsim, a, b, seed = args
    n, lo, hi, follow = DESIGNS[design]
    kind, r, s = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1, c2 = chi2.ppf(0.95, 1), chi2.ppf(0.95, 2)
    any_level = any_2df = slope_flag = picks_rate = failed = 0
    s_hat, deaths = [], []
    for _ in range(nsim):
        data = simulate(n, lo, hi, follow, r, s, a, b, rng)
        ll = fit_all(data, a, b)
        failed += not all(v[2] for v in ll.values())
        any_level += 2 * (ll["level"][0] - ll["null"][0]) > c1       # the usual "did it cut deaths" test
        any_2df += 2 * (ll["full"][0] - ll["null"][0]) > c2
        slope_flag += 2 * (ll["full"][0] - ll["level"][0]) > c1      # fingerprint: did the slope change?
        picks_rate += ll["rate"][0] > ll["level"][0]                 # which one-parameter story fits better
        s_hat.append(float(np.exp(ll["full"][1][3])))
        deaths.append(float(data[3].sum()))
    q = np.percentile(s_hat, [5, 50, 95])
    return dict(design=design, effect=effect, kind=kind, true_r=r, true_s=s, n_trials=nsim,
                mean_deaths=float(np.mean(deaths)), fits_not_converged=failed / nsim,
                power_any_effect=any_level / nsim, power_any_effect_2df=any_2df / nsim,
                slope_change_flagged=slope_flag / nsim, rate_story_preferred=picks_rate / nsim,
                slope_estimate_p5=float(q[0]), slope_estimate_median=float(q[1]), slope_estimate_p95=float(q[2]))


def years_gained(x0, r, s, a, b):
    """Extra remaining life expectancy at entry age x0 under an effect (r, s)."""
    def e(rr, ss):
        h0, bs = a * np.exp(b * x0) * rr, b * ss
        return quad(lambda t: np.exp(-h0 * (np.exp(bs * t) - 1.0) / bs), 0, 80, limit=200)[0]
    return e(r, s) - e(1.0, 1.0)


def check_gradient(a, b):
    """Analytic gradient against finite differences on one simulated trial."""
    rng = np.random.default_rng(1)
    data = simulate(400, 7, 12, 3, 0.8, 0.9, a, b, rng)
    theta = np.array([np.log(a) + 0.1, np.log(b) - 0.05, -0.2, -0.1])
    _, g = negll_and_grad(theta, *data, True, True)
    num = []
    for i in range(4):
        d = np.zeros(4)
        d[i] = 1e-6
        num.append((negll_and_grad(theta + d, *data, True, True)[0]
                    - negll_and_grad(theta - d, *data, True, True)[0]) / 2e-6)
    return float(np.max(np.abs(g - np.array(num)) / (np.abs(num) + 1e-9)))


def main():
    quick = "--quick" in sys.argv
    nsim = 40 if quick else 400
    a, b = dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    cases = [(d, e, nsim, a, b, SEED + 1000 * i + j)
             for i, d in enumerate(DESIGNS) for j, e in enumerate(EFFECTS)]
    if quick:
        cases = [c for c in cases if c[0].startswith(("580", "1,300 dogs, 4"))]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases)
    gains = {e: {f"entry age {x0:g}": years_gained(x0, r, s, a, b) for x0 in (7.0, 10.0, 12.0)}
             for e, (_, r, s) in EFFECTS.items() if e != "no effect"}
    out = dict(
        note="Simulated trials. Dog mortality: Gompertz fit to the UK life table. Entry ages are assumptions. "
             "Effects start at enrolment and persist. 'slope_change_flagged' is the share of simulated trials in "
             "which a 5% likelihood-ratio test rejects 'the slope did not change'.",
        dog_a_per_year=a, dog_b_per_year=b, seed=SEED, trials_per_case=nsim,
        gradient_check_max_relative_error=grad_err,
        years_gained_by_effect=gains, cases=rows)
    name = "results_quick.json" if quick else "results.json"
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(out, fh, indent=1)
    for r_ in rows:
        print(f"{r_['design'][:28]:28s} | {r_['effect']:17s} | deaths {r_['mean_deaths']:6.0f} | any effect "
              f"{r_['power_any_effect']:.2f} | slope flagged {r_['slope_change_flagged']:.2f} | "
              f"rate story {r_['rate_story_preferred']:.2f} | s-hat {r_['slope_estimate_median']:.2f} "
              f"({r_['slope_estimate_p5']:.2f}-{r_['slope_estimate_p95']:.2f})")


if __name__ == "__main__":
    main()
