"""Experiment 2: fit the senescent-cell (saturating-removal) model to real dog life tables.

Question. In the model of Karin et al. 2019, can a dog be a human whose damage production simply
climbs faster ("speed alone"), or must something else differ?

Model. Senescent-cell load X:  dX = [eta*t - beta*X/(kappa+X)] dt + sqrt(2*eps) dW, death at X = Xc.
Hazard at age t is the quasi-static first-passage rate for production eta*t, plus a constant
background hazard lam. In hypotheses A to D, kappa and Xc are held at the published human values
(0.5 and 17), as Karin et al. do when moving between mouse and human. Published human values: eta = 0.00135 /day/yr,
beta = 0.15 /day, eps = 0.142 /day.

Hypotheses, fitted to the UK dog life table (Teng et al. 2022, deaths and numbers at risk, ages 2-18)
by binomial maximum likelihood:
    A  speed alone:            eta free; beta, eps at human values
    B  speed + noise:          eta, eps free
    C  speed + clearance:      eta, beta free
    D  speed + both:           eta, beta, eps free
    E  speed + threshold:      eta, Xc free      (added after the independent check)
    F  speed + saturation:     eta, kappa free   (added after the independent check)
Every hypothesis also has a free background hazard lam.
The same four are refitted to the US table (Montoya et al. 2023) by least squares on log hazard
as a robustness check.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "damage_repair"))
from saturating_removal import HUMAN_PUB, MOUSE_PUB, qs_hazard  # noqa: E402

KAPPA, XC = HUMAN_PUB["kappa"], HUMAN_PUB["Xc"]
H = HUMAN_PUB
STEP = 0.1
GRID = np.arange(2.0, 19.0001, STEP)          # ages at which the model hazard is evaluated
AGES = np.arange(2, 19)                        # one-year intervals [x, x+1), x = 2..18

HYPOTHESES = {
    "A speed alone": ("eta",),
    "B speed + noise": ("eta", "eps"),
    "C speed + clearance": ("eta", "beta"),
    "D speed + clearance + noise": ("eta", "beta", "eps"),
    "E speed + lethal threshold": ("eta", "Xc"),
    "F speed + saturation point": ("eta", "kappa"),
}


def unpack(theta, free):
    """theta = log multipliers of the free parameters relative to human, then log(lam)."""
    mult = dict(eta=1.0, beta=1.0, eps=1.0, Xc=1.0, kappa=1.0)
    for name, v in zip(free, theta[:-1]):
        mult[name] = float(np.exp(v))
    return mult, float(np.exp(theta[-1]))


def interval_hazard(mult, lam, ages=AGES):
    """Integrated hazard over each one-year interval x..x+1 for x in ages."""
    h = qs_hazard(GRID, H["eta"] * mult["eta"], H["beta"] * mult["beta"], KAPPA * mult["kappa"],
                  H["eps"] * mult["eps"], XC * mult["Xc"], lam)
    cum = cumulative_trapezoid(h, GRID, initial=0.0)
    idx = np.round((ages - GRID[0]) / STEP).astype(int)
    idx1 = np.round((ages + 1 - GRID[0]) / STEP).astype(int)
    return cum[idx1] - cum[idx]


def nll_binomial(theta, free, dx, lx):
    mult, lam = unpack(theta, free)
    with np.errstate(all="ignore"):
        Hx = interval_hazard(mult, lam)
        q = -np.expm1(-Hx)
        if not np.all(np.isfinite(q)) or np.any(q <= 0) or np.any(q >= 1):
            return 1e12
        return float(-np.sum(dx * np.log(q) + (lx - dx) * np.log1p(-q)))


def sse_loghazard(theta, free, m_obs, ages):
    mult, lam = unpack(theta, free)
    with np.errstate(all="ignore"):
        Hx = interval_hazard(mult, lam, ages)
        if not np.all(np.isfinite(Hx)) or np.any(Hx <= 0):
            return 1e12
        return float(np.sum((np.log(Hx) - np.log(m_obs)) ** 2))


def fit(objective, free, args):
    best = None
    for eta0 in (3.0, 6.0, 10.0):
        for other0 in (0.0, 0.6, -0.6):
            x0 = [np.log(eta0)] + [other0] * (len(free) - 1) + [np.log(0.005)]
            res = minimize(objective, x0, args=(free,) + args, method="Nelder-Mead",
                           options=dict(xatol=1e-5, fatol=1e-6, maxiter=4000, maxfev=4000))
            if best is None or res.fun < best.fun:
                best = res
            if len(free) == 1:
                break
    return best


def numerical_se(objective, theta, free, args, h=1e-3):
    """Standard errors of the log multipliers from the observed information (sampling error only)."""
    k = len(theta)
    hess = np.zeros((k, k))
    f0 = objective(theta, free, *args)
    for i in range(k):
        for j in range(i, k):
            e_i, e_j = np.eye(k)[i] * h, np.eye(k)[j] * h
            if i == j:
                hess[i, i] = (objective(theta + e_i, free, *args) - 2 * f0
                              + objective(theta - e_i, free, *args)) / h ** 2
            else:
                hess[i, j] = hess[j, i] = (objective(theta + e_i + e_j, free, *args)
                                           - objective(theta + e_i - e_j, free, *args)
                                           - objective(theta - e_i + e_j, free, *args)
                                           + objective(theta - e_i - e_j, free, *args)) / (4 * h ** 2)
    try:
        cov = np.linalg.inv(hess)
        return np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        return np.full(k, np.nan)


def describe(mult, lam):
    """Derived quantities a lab could measure, and the mortality the fit implies."""
    ages = np.arange(2.0, 19.01, 0.05)
    h = qs_hazard(ages, H["eta"] * mult["eta"], H["beta"] * mult["beta"], KAPPA * mult["kappa"],
                  H["eps"] * mult["eps"], XC * mult["Xc"], lam)
    sel = (ages >= 4) & (ages <= 14)
    slope = np.polyfit(ages[sel], np.log(h[sel] - lam), 1)[0]
    half_life_days = np.log(2) * KAPPA * mult["kappa"] / (H["beta"] * mult["beta"])       # clearance half-life at low load
    return dict(doubling_time_years_4_to_14=float(np.log(2) / slope),
                senescent_cell_half_life_days_at_low_load=float(half_life_days),
                human_half_life_days_at_low_load=float(np.log(2) * KAPPA / H["beta"]),
                age_when_production_reaches_clearance_ceiling=float(H["beta"] * mult["beta"] / (H["eta"] * mult["eta"])))


def main():
    out = dict(human_published=H, mouse_published=MOUSE_PUB, kappa=KAPPA, Xc=XC,
               note=("Multipliers are relative to the published human values. Standard errors reflect sampling "
                     "error in the life table only; the choice of model matters far more."))

    # ---------------- UK table: deaths and numbers at risk -> binomial likelihood
    uk = pd.read_csv(os.path.join(HERE, "..", "dog_lifespan", "data", "teng2022_uk_dog_life_table.csv"), comment="#")
    uk = uk[(uk.age_start >= 2) & (uk.age_start <= 18)]
    dx, lx = uk.dx.to_numpy(float), uk.lx.to_numpy(float)
    q_obs = dx / lx
    fits = {}
    for name, free in HYPOTHESES.items():
        res = fit(nll_binomial, free, (dx, lx))
        mult, lam = unpack(res.x, free)
        se = numerical_se(nll_binomial, res.x, free, (dx, lx))
        Hx = interval_hazard(mult, lam)
        q_fit = -np.expm1(-Hx)
        k = len(free) + 1
        fits[name] = dict(
            multipliers=mult, background_hazard_per_year=lam,
            log_multiplier_se={n: float(s) for n, s in zip(free, se[:-1])},
            neg_loglik=float(res.fun), parameters=k, aic=float(2 * k + 2 * res.fun),
            max_relative_error_in_death_probability=float(np.max(np.abs(q_fit / q_obs - 1))),
            median_relative_error_in_death_probability=float(np.median(np.abs(q_fit / q_obs - 1))),
            fitted_q=[float(v) for v in q_fit], **describe(mult, lam))
        print(f"UK  {name:30s} mult {dict((k_, round(v, 3)) for k_, v in mult.items())} lam {lam:.4f} "
              f"nll {res.fun:.1f} maxerr {fits[name]['max_relative_error_in_death_probability']:.2f}")
    a_nll = fits["A speed alone"]["neg_loglik"]
    for name in list(HYPOTHESES)[1:]:
        df = fits[name]["parameters"] - fits["A speed alone"]["parameters"]
        stat = 2 * (a_nll - fits[name]["neg_loglik"])
        fits[name]["lr_test_vs_speed_alone"] = dict(statistic=float(stat), df=int(df), p_value=float(chi2.sf(stat, df)))
    out["uk_dogs_binomial_fit"] = dict(ages=[int(a) for a in AGES], observed_q=[float(v) for v in q_obs],
                                       deaths=int(dx.sum()), fits=fits)

    # ---------------- US table: death probabilities only -> least squares on log hazard
    us = pd.read_csv(os.path.join(HERE, "..", "dog_lifespan", "data", "montoya2023_us_dog_life_table.csv"), comment="#")
    us = us[(us.age_start >= 2) & (us.age_start <= 15)]
    ages_us = us.age_start.to_numpy(int)
    m_obs = -np.log1p(-us.qx.to_numpy(float))
    us_fits = {}
    for name, free in HYPOTHESES.items():
        res = fit(sse_loghazard, free, (m_obs, ages_us))
        mult, lam = unpack(res.x, free)
        us_fits[name] = dict(multipliers=mult, background_hazard_per_year=lam,
                             rms_log_error=float(np.sqrt(res.fun / len(m_obs))), **describe(mult, lam))
        print(f"US  {name:30s} mult {dict((k_, round(v, 3)) for k_, v in mult.items())} lam {lam:.4f} "
              f"rms {us_fits[name]['rms_log_error']:.3f}")
    out["us_dogs_least_squares_fit"] = dict(ages=[int(a) for a in us.age_start], fits=us_fits)

    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    # ---------------- figure and its data
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(AGES + 0.5, q_obs, "ko", label="UK dogs, observed (Teng 2022)")
    table = {"age": AGES, "observed_q": q_obs}
    for name in HYPOTHESES:
        ax.semilogy(AGES + 0.5, fits[name]["fitted_q"], label=name)
        table[name] = fits[name]["fitted_q"]
    ax.set_xlabel("age (years)")
    ax.set_ylabel("probability of dying within the year")
    ax.set_title("Senescent-cell model fitted to dogs: faster production alone does not fit")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "dog_sr_fits.png"), dpi=150)
    with open(os.path.join(HERE, "data", "dog_sr_fits.csv"), "w") as fh:
        fh.write("# Observed UK dog death probabilities (Teng et al. 2022, https://pmc.ncbi.nlm.nih.gov/articles/PMC9050668) "
                 "and four fits of the saturating-removal model\n")
        pd.DataFrame(table).to_csv(fh, index=False)


if __name__ == "__main__":
    main()
