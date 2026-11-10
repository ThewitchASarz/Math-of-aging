"""Experiment 006: a two-period analysis for trials in which treatment stops.

Data are simulated exactly as in experiment 003 (one year of treatment, then observation).
Analysis model for a treated dog, t = years since enrolment, tau = 1:
    hazard = a * exp(b * (x0 + t)) * r_on      while t <  tau
    hazard = a * exp(b * (x0 + t)) * r_off     once  t >= tau
Placebo dogs have r_on = r_off = 1. a and b are estimated from the data.

Lasting-benefit test: likelihood-ratio test of r_off = 1 (1 degree of freedom, 5% level).

Run:  python3 two_period.py            (about 1 minute on 2 cores)
      python3 two_period.py --quick
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
_spec = importlib.util.spec_from_file_location(
    "fading", os.path.join(HERE, "..", "003_fading_effect", "fading.py"))
fading = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fading)
fp, TAU, DESIGNS, EFFECTS = fading.fp, fading.TAU, fading.DESIGNS, fading.EFFECTS

SEED = 20261113
MODELS = {"none": (False, False), "on_only": (True, False), "both": (True, True)}


def negll_and_grad(theta, x0, arm, t, died, has_on, has_off):
    """theta = [ln a, ln b, (ln r_on), (ln r_off)]."""
    lna, lnb = theta[0], theta[1]
    k = 2
    ln_on = theta[k] if has_on else 0.0
    k += int(has_on)
    ln_off = theta[k] if has_off else 0.0
    b = np.exp(lnb)
    eA = np.exp(lna + b * x0)
    after = t > TAU
    t1 = np.minimum(t, TAU)
    E1, Et, Etau = np.exp(b * t1), np.exp(b * t), np.exp(b * TAU)
    C1 = eA * (E1 - 1.0) / b                                   # cumulative baseline hazard while on treatment
    C2 = np.where(after, eA * (Et - Etau) / b, 0.0)            # and after treatment
    r_on, r_off = np.exp(ln_on * arm), np.exp(ln_off * arm)
    H = r_on * C1 + r_off * C2
    logh = lna + b * (x0 + t) + np.where(after, ln_off, ln_on) * arm
    ll = np.sum(died * logh - H)
    dC1 = x0 * C1 + eA * (t1 * E1 / b - (E1 - 1.0) / b ** 2)
    dC2 = np.where(after, x0 * C2 + eA * ((t * Et - TAU * Etau) / b - (Et - Etau) / b ** 2), 0.0)
    grad = [np.sum(died - H), b * np.sum(died * (x0 + t) - r_on * dC1 - r_off * dC2)]
    if has_on:
        grad.append(np.sum(arm * (died * (~after) - r_on * C1)))
    if has_off:
        grad.append(np.sum(arm * (died * after - r_off * C2)))
    return -ll, -np.array(grad)


def fit(model, data, a, b, extra=None):
    has_on, has_off = MODELS[model]
    n_extra = int(has_on) + int(has_off)
    start = [np.log(a), np.log(b)] + (list(extra) if extra is not None else [0.0] * n_extra)
    bounds = [(np.log(a) - 3, np.log(a) + 3), (np.log(b) - 1, np.log(b) + 1)] + [(np.log(0.05), np.log(5.0))] * n_extra
    with np.errstate(over="ignore", invalid="ignore"):
        res = minimize(negll_and_grad, start, args=(*data, has_on, has_off), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
    return -res.fun, res.x


def check_gradient(a, b):
    rng = np.random.default_rng(2)
    data = fading.simulate(500, 7, 12, 3, 0.7, 1.0, 0.8, a, b, rng)
    theta = np.array([np.log(a) + 0.1, np.log(b) - 0.05, -0.2, -0.1])
    _, g = negll_and_grad(theta, *data, True, True)
    num = []
    for i in range(4):
        d = np.zeros(4)
        d[i] = 1e-6
        num.append((negll_and_grad(theta + d, *data, True, True)[0]
                    - negll_and_grad(theta - d, *data, True, True)[0]) / 2e-6)
    return float(np.max(np.abs(g - np.array(num)) / (np.abs(num) + 1e-9)))


def one_case(args):
    design, effect, nsim, a, b, seed = args
    n, lo, hi, follow = DESIGNS[design]
    r1, s1, r2 = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1, c2 = chi2.ppf(0.95, 1), chi2.ppf(0.95, 2)
    lasting = overall = 0
    off_hat, on_hat = [], []
    for _ in range(nsim):
        data = fading.simulate(n, lo, hi, follow, r1, s1, r2, a, b, rng)
        ll_none = fit("none", data, a, b)
        ll_on = fit("on_only", data, a, b)
        ll_both = fit("both", data, a, b, extra=[ll_on[1][2], 0.0])
        lasting += 2 * (ll_both[0] - ll_on[0]) > c1
        overall += 2 * (ll_both[0] - ll_none[0]) > c2
        on_hat.append(float(np.exp(ll_both[1][2])))
        off_hat.append(float(np.exp(ll_both[1][3])))
    return dict(design=design, effect=effect, n_trials=nsim, lasting_benefit_reported=lasting / nsim,
                any_effect_reported=overall / nsim, median_ratio_on_treatment=float(np.median(on_hat)),
                median_ratio_after_treatment=float(np.median(off_hat)))


def main():
    quick = "--quick" in sys.argv
    nsim = 40 if quick else 400
    a, b = fp.dog_parameters()
    grad_err = check_gradient(a, b)
    assert grad_err < 1e-5, grad_err
    cases = [(d, e, nsim, a, b, SEED + 1000 * i + j)
             for i, d in enumerate(DESIGNS) for j, e in enumerate(EFFECTS)]
    if quick:
        cases = [c for c in cases if c[0].startswith("580 dogs, 3")]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases)
    out = dict(note="Trials simulated as in experiment 003, analysed with separate risk ratios on and after "
                    "treatment. 'lasting_benefit_reported' is the share of trials in which a 5% likelihood-ratio "
                    "test rejects 'the after-treatment ratio is 1'.",
               seed=SEED, trials_per_case=nsim, gradient_check_max_relative_error=grad_err, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['design'][:20]:20s} | {r['effect']:42s} | lasting {r['lasting_benefit_reported']:.2f} | "
              f"any {r['any_effect_reported']:.2f} | on {r['median_ratio_on_treatment']:.2f} "
              f"after {r['median_ratio_after_treatment']:.2f}")


if __name__ == "__main__":
    main()
