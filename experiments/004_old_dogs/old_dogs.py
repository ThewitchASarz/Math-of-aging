"""Experiment 004: does the doubling law hold for old dogs, and does hidden frailty break the fingerprint test?

Part 1  UK life table (Teng et al. 2022; deaths d_x among l_x dogs alive at age x, ages 2-19), fitted by
        binomial maximum likelihood with
          plain doubling law      hazard = a * exp(b * x)
          doubling law + frailty  each dog's hazard is z * a * exp(b * x), z ~ gamma(mean 1, variance v) among
                                  dogs alive at age 2; population hazard then levels off at old ages.
Part 2  Plain law fitted on ages 2-13, forecast compared with observed rates at ages 14-19.
Part 3  US table (rates only): curvature of log death rate against age.
Part 4  Trials simulated with hidden frailty of the fitted size, analysed with the plain models of
        experiment 001. How often is a pure risk cut misread as a slope change?

Run:  python3 old_dogs.py            (about 4 minutes on 2 cores)
      python3 old_dogs.py --quick
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

SEED = 20261112
AGE0 = 2.0                                        # frailty is defined among dogs alive at this age


def load(name):
    path = os.path.join(ROOT, "dog_lifespan", "data", name)
    rows = [l.strip().split(",") for l in open(path) if not l.startswith("#")]
    return {c: np.array([float(r[i]) for r in rows[1:]]) for i, c in enumerate(rows[0])}


def cum_hazard(x, a, b):
    """Baseline cumulative hazard from AGE0 to age x."""
    return (a / b) * (np.exp(b * x) - np.exp(b * AGE0))


def q_plain(ages, a, b):
    return 1.0 - np.exp(-(cum_hazard(ages + 1, a, b) - cum_hazard(ages, a, b)))


def q_frailty(ages, a, b, v):
    surv = lambda x: (1.0 + v * cum_hazard(x, a, b)) ** (-1.0 / v)
    return 1.0 - surv(ages + 1) / surv(ages)


def binom_negll(q, d, l):
    q = np.clip(q, 1e-12, 1 - 1e-12)
    return -np.sum(d * np.log(q) + (l - d) * np.log(1 - q))


def fit_plain(ages, d, l):
    f = lambda th: binom_negll(q_plain(ages, np.exp(th[0]), np.exp(th[1])), d, l)
    res = minimize(f, [np.log(0.005), np.log(0.3)], method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-8, "maxiter": 4000})
    return np.exp(res.x), -res.fun


def fit_frailty(ages, d, l):
    f = lambda th: binom_negll(q_frailty(ages, np.exp(th[0]), np.exp(th[1]), np.exp(th[2])), d, l)
    best = None
    for v0 in (0.05, 0.3, 1.0):
        res = minimize(f, [np.log(0.004), np.log(0.35), np.log(v0)], method="Nelder-Mead",
                       options={"xatol": 1e-10, "fatol": 1e-8, "maxiter": 8000})
        if best is None or res.fun < best.fun:
            best = res
    return np.exp(best.x), -best.fun


def simulate_frail(n, lo, hi, follow, r, s, a, b, v, rng):
    """Like experiment 001, with a hidden personal risk factor z. Among dogs alive at entry age x0 the factor
    is gamma with shape 1/v and mean 1 / (1 + v * H(x0)): frailer dogs have already been thinned out."""
    x0 = rng.uniform(lo, hi, n)
    arm = (np.arange(n) % 2).astype(float)
    z = rng.gamma(1.0 / v, v / (1.0 + v * cum_hazard(x0, a, b)))
    h0 = z * a * np.exp(b * x0) * np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    t = np.log1p(-np.log(rng.random(n)) * bs / h0) / bs
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


DESIGNS = {
    "580 dogs, 3 years": (580, 7.0, 12.0, 3.0),
    "1,300 dogs, 4 years": (1300, 10.0, 14.0, 4.0),
    "1,300 dogs, 10 years": (1300, 7.0, 12.0, 10.0),
    "5,000 dogs, 3 years": (5000, 7.0, 12.0, 3.0),
}
EFFECTS = {"no effect": (1.0, 1.0), "risk cut 30%": (0.7, 1.0), "risk cut 50%": (0.5, 1.0),
           "aging slowed 20%": (1.0, 0.8)}


def one_case(args):
    design, effect, nsim, a, b, v, start_a, start_b, seed = args
    n, lo, hi, follow = DESIGNS[design]
    r, s = EFFECTS[effect]
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    any_level = flagged = 0
    s_hat = []
    for _ in range(nsim):
        data = simulate_frail(n, lo, hi, follow, r, s, a, b, v, rng)
        ll = fp.fit_all(data, start_a, start_b)
        any_level += 2 * (ll["level"][0] - ll["null"][0]) > c1
        flagged += 2 * (ll["full"][0] - ll["level"][0]) > c1
        s_hat.append(float(np.exp(ll["full"][1][3])))
    return dict(design=design, effect=effect, n_trials=nsim, detected_by_usual_test=any_level / nsim,
                slope_change_flagged=flagged / nsim, median_slope_estimate=float(np.median(s_hat)))


def main():
    quick = "--quick" in sys.argv
    uk = load("teng2022_uk_dog_life_table.csv")
    keep = (uk["age_start"] >= 2) & (uk["age_start"] <= 19)
    ages, d, l = uk["age_start"][keep], uk["dx"][keep], uk["lx"][keep]

    # Part 1
    (a_p, b_p), ll_p = fit_plain(ages, d, l)
    (a_f, b_f, v_f), ll_f = fit_frailty(ages, d, l)
    lr = 2 * (ll_f - ll_p)
    part1 = dict(plain=dict(a=float(a_p), b=float(b_p), doubling_time_years=float(np.log(2) / b_p), loglik=float(ll_p)),
                 with_frailty=dict(a=float(a_f), b=float(b_f), frailty_variance=float(v_f),
                                   individual_doubling_time_years=float(np.log(2) / b_f),
                                   population_plateau_per_year=float(b_f / v_f), loglik=float(ll_f)),
                 likelihood_ratio=float(lr), p_value=float(0.5 * chi2.sf(lr, 1)),
                 p_value_note="Halved because the frailty variance cannot be negative.",
                 deaths_used=float(d.sum()))
    obs = -np.log(1 - d / l)
    part1["by_age"] = [dict(age=int(x), observed=float(o), plain=float(-np.log(1 - qp)), with_frailty=float(-np.log(1 - qf)))
                       for x, o, qp, qf in zip(ages, obs, q_plain(ages, a_p, b_p), q_frailty(ages, a_f, b_f, v_f))]

    # Part 2
    young = ages <= 13
    (a_y, b_y), _ = fit_plain(ages[young], d[young], l[young])
    forecast = -np.log(1 - q_plain(ages, a_y, b_y))
    part2 = dict(fit_ages="2-13", a=float(a_y), b=float(b_y), doubling_time_years=float(np.log(2) / b_y),
                 forecast_over_observed={str(int(x)): float(f / o) for x, f, o in zip(ages, forecast, obs) if x >= 14})

    # Part 3
    us = load("montoya2023_us_dog_life_table.csv")
    k = (us["age_start"] >= 2) & (us["age_start"] <= 16)
    x_us, h_us = us["age_start"][k] + 0.5, -np.log(1 - us["qx"][k])
    quad = np.polyfit(x_us, np.log(h_us), 2)
    yo = us["age_start"][k] <= 13
    lin = np.polyfit(x_us[yo], np.log(h_us[yo]), 1)
    part3 = dict(quadratic_coefficient_per_year2=float(quad[0]),
                 local_doubling_time_at_age_4=float(np.log(2) / (2 * quad[0] * 4 + quad[1])),
                 local_doubling_time_at_age_15=float(np.log(2) / (2 * quad[0] * 15 + quad[1])),
                 forecast_over_observed_from_fit_2_13={str(int(a_)): float(np.exp(np.polyval(lin, a_ + 0.5)) / h)
                                                       for a_, h in zip(us["age_start"][k], h_us) if a_ >= 14})

    # Part 4
    nsim = 30 if quick else 400
    cases = [(dsg, eff, nsim, a_f, b_f, v_f, a_p, b_p, SEED + 1000 * i + j)
             for i, dsg in enumerate(DESIGNS) for j, eff in enumerate(EFFECTS)]
    if quick:
        cases = [c for c in cases if c[0].startswith(("580", "1,300 dogs, 10"))]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases)
    out = dict(note="Parts 1-3 use the published dog life tables. Part 4 simulates trials with hidden frailty of the "
                    "size fitted in part 1 and analyses them with the plain models of experiment 001.",
               seed=SEED, part1_uk_fit=part1, part2_forecast_old_ages=part2, part3_us_curvature=part3,
               part4_trials_with_frailty=dict(trials_per_case=nsim, cases=rows))
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("plain:   b %.4f doubling %.2f" % (b_p, np.log(2) / b_p))
    print("frailty: b %.4f doubling %.2f variance %.3f plateau %.2f/yr  LR %.1f p %.2g" %
          (b_f, np.log(2) / b_f, v_f, b_f / v_f, lr, part1["p_value"]))
    print("forecast/observed (UK, fit 2-13):", {k_: round(v, 2) for k_, v in part2["forecast_over_observed"].items()})
    print("US quad coef %.4f; local doubling at 4: %.2f, at 15: %.2f; forecast/obs %s" %
          (quad[0], part3["local_doubling_time_at_age_4"], part3["local_doubling_time_at_age_15"],
           {k_: round(v, 2) for k_, v in part3["forecast_over_observed_from_fit_2_13"].items()}))
    for r_ in rows:
        print(f"{r_['design']:22s} | {r_['effect']:17s} | usual test {r_['detected_by_usual_test']:.2f} | "
              f"slope flagged {r_['slope_change_flagged']:.2f} | s-hat {r_['median_slope_estimate']:.2f}")


if __name__ == "__main__":
    main()
