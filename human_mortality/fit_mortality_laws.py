"""Part 1: fit Gompertz, Gompertz-Makeham, gamma-Gompertz-Makeham (logistic/Perks family) and Kannisto
models to the SSA 2023 US period life table, by sex. Verify closed forms for survival, life expectancy,
modal age and doubling time. Writes results.json['fits'], figures/fit_*.png, data/fit_*.csv.

Fitting method: nonlinear least squares on log of the one-year hazard m_x = -ln(1-q_x), where the model
prediction is the exact integrated hazard H(x+1)-H(x). Confidence intervals: moving-block bootstrap of
log-residuals (block length 5 years, 2000 reps, seed 1) because residuals of a smoothed table are
autocorrelated. These CIs measure how well-determined the parameters are given the shape of this table;
they are NOT sampling uncertainty from death counts (the table does not report counts).
"""
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import least_squares, brentq

from common import (load_ssa, update_results, H_gm, mu_gm, H_ggm, mu_ggm, H_kannisto, mu_kannisto,
                    life_expectancy_from_mu, lifetable_e0, SSA_URL, SSA_DESC, FIG, DATA)

rng = np.random.default_rng(1)
d = load_ssa()
ages_all = d.index.values.astype(float)

MODELS = {
    # name: (param names, initial guess, log-param fit?, age range, interval-hazard function)
    "gompertz": (["a", "b"], [3e-5, 0.09], (40, 95),
                 lambda x, p: H_gm(x + 1, p[0], p[1]) - H_gm(x, p[0], p[1])),
    "gompertz_makeham": (["a", "b", "c"], [2e-5, 0.095, 5e-4], (30, 95),
                         lambda x, p: H_gm(x + 1, *p) - H_gm(x, *p)),
    # gamma-Gompertz (Vaupel frailty / Perks-type logistic) fitted to old ages only: below ~90 the US
    # 2023 log-hazard is concave-up (slope rises from ~0.05 at 40 to ~0.12 at 90), which frailty cannot
    # produce, so a 30-105 fit drives s2 to 0. Makeham term omitted (small relative to hazard at 70+).
    "gamma_gompertz": (["a", "b", "s2"], [2e-5, 0.11, 0.05], (70, 110),
                       lambda x, p: H_ggm(x + 1, p[0], p[1], 0.0, p[2]) - H_ggm(x, p[0], p[1], 0.0, p[2])),
    "kannisto": (["a", "b"], [0.05, 0.11], (80, 110),
                 lambda x, p: H_kannisto(x + 1, *p) - H_kannisto(x, *p)),
}


def fit(model, sex, x=None, y=None):
    names, p0, (lo, hi), f = MODELS[model]
    if x is None:
        sel = (ages_all >= lo) & (ages_all <= hi)
        x = ages_all[sel]
        y = np.log(d[f"mu_{sex}"].values[sel])
    def resid(lp):
        return np.log(f(x, np.exp(lp))) - y
    r = least_squares(resid, np.log(p0), method="lm", x_scale="jac", max_nfev=20000)
    return np.exp(r.x), r.fun, x, y


def block_bootstrap(model, sex, p_hat, res, x, y, reps=2000, L=5):
    f = MODELS[model][3]
    fitted = np.log(f(x, p_hat))
    n = len(res)
    out = []
    for _ in range(reps):
        idx = np.concatenate([np.arange(s, s + L) for s in rng.integers(0, n - L + 1, size=n // L + 1)])[:n]
        yb = fitted + res[idx]
        pb, _, _, _ = fit(model, sex, x, yb)
        out.append(pb)
    return np.array(out)


# ---------- closed-form checks ----------
def e0_gompertz_closed(a, b):
    """e0 = (1/b) e^{a/b} E1(a/b)."""
    z = a / b
    return float(mp.e ** z * mp.e1(z) / b)


def e0_gm_closed(a, b, c):
    """e0 = (1/b) e^{a/b} (a/b)^{c/b} Gamma(-c/b, a/b)  (upper incomplete gamma)."""
    z = a / b
    return float(mp.e ** z * z ** (c / b) * mp.gammainc(-c / b, z) / b)


def mode_gm(a, b, c=0.0):
    """Modal age of death density f = mu S: solve mu' = mu^2, i.e. a b e^{bx} = (c + a e^{bx})^2."""
    g = lambda x: a * b * np.exp(b * x) - (c + a * np.exp(b * x)) ** 2
    return brentq(g, 20, 150)


results = {"source": SSA_DESC, "source_url": SSA_URL, "data_year": 2023,
           "units": {"a": "per year (hazard at age 0 implied by the adult fit)", "b": "per year of age",
                     "c": "per year", "s2": "dimensionless (variance of gamma frailty)",
                     "ages": "years"},
           "method": "NLS on log one-year hazard; moving-block bootstrap CIs (L=5, 2000 reps, seed 1). "
                     "CIs reflect fit determinacy to a smoothed table, not death-count sampling error.",
           "models": {}}

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
curves = {"age": ages_all}
for k, sex in enumerate(["male", "female"]):
    ax = axes[k]
    mu = d[f"mu_{sex}"].values
    curves[f"observed_{sex}"] = mu
    ax.semilogy(ages_all, mu, "k.", ms=4, label="SSA 2023 period table (observed/smoothed)")
    results["models"][sex] = {}
    for model, style in zip(MODELS, ["-", "--", "-.", ":"]):
        names, _, (lo, hi), f = MODELS[model]
        p, res, x, y = fit(model, sex)
        boot = block_bootstrap(model, sex, p, res, x, y)
        ci = np.percentile(boot, [2.5, 97.5], axis=0)
        # goodness of fit inside range and across all ages
        pred_all = np.log(f(ages_all, p))
        rmse_in = float(np.sqrt(np.mean(res ** 2)))
        n = len(res); kpar = len(p)
        aic = float(n * np.log(np.mean(res ** 2)) + 2 * kpar)
        entry = {"fit_ages": [lo, hi], "n_ages": n,
                 "params": {nm: float(v) for nm, v in zip(names, p)},
                 "ci95": {nm: [float(ci[0, i]), float(ci[1, i])] for i, nm in enumerate(names)},
                 "rmse_log_hazard_in_range": rmse_in,
                 "max_abs_pct_error_in_range": float(100 * np.max(np.abs(np.expm1(res)))),
                 "aic_gaussian_log": aic}
        # error by age band outside range (where the model fails)
        bands = {"0": (0, 0), "1-14": (1, 14), "15-29": (15, 29), "30-39": (30, 39), "40-79": (40, 79),
                 "80-94": (80, 94), "95-104": (95, 104), "105-119": (105, 119)}
        entry["mean_pct_error_by_age_band"] = {
            bnm: float(100 * np.mean(np.expm1(pred_all[(ages_all >= b0) & (ages_all <= b1)]
                                              - np.log(mu[(ages_all >= b0) & (ages_all <= b1)]))))
            for bnm, (b0, b1) in bands.items()}
        if model in ("gompertz", "gompertz_makeham", "gamma_gompertz"):
            b = p[1]
            entry["doubling_time_years"] = float(np.log(2) / b)
            entry["doubling_time_ci95"] = [float(np.log(2) / ci[1, 1]), float(np.log(2) / ci[0, 1])]
        if model == "gompertz":
            a, b = p
            e_closed = e0_gompertz_closed(a, b)
            e_num, _, _ = life_expectancy_from_mu(lambda t: mu_gm(t, a, b))
            entry["e0_if_law_held_from_birth"] = e_closed
            entry["e0_numeric_check"] = e_num
            entry["e0_approx_(ln(b/a)-gamma)/b"] = float((np.log(b / a) - np.euler_gamma) / b)
            entry["modal_age_closed_ln(b/a)/b"] = float(np.log(b / a) / b)
            entry["modal_age_numeric"] = mode_gm(a, b)
        if model == "gompertz_makeham":
            a, b, c = p
            e_closed = e0_gm_closed(a, b, c)
            e_num, _, _ = life_expectancy_from_mu(lambda t: mu_gm(t, a, b, c))
            entry["e0_if_law_held_from_birth"] = e_closed
            entry["e0_numeric_check"] = e_num
            entry["modal_age_numeric"] = mode_gm(a, b, c)
        if model == "gamma_gompertz":
            a, b, s2 = p
            entry["late_life_plateau_b/s2_per_year"] = float(b / s2)
            pl = boot[:, 1] / boot[:, 2]
            entry["plateau_ci95"] = [float(np.percentile(pl, 2.5)), float(np.percentile(pl, 97.5))]
            entry["plateau_annual_death_prob"] = float(-np.expm1(-b / s2))
            entry["note"] = ("Fitted 70-110; SSA rates above ~100 are actuarial extrapolations, so the implied "
                             "plateau describes SSA's extrapolation, not observed data.")
        if model == "kannisto":
            entry["asymptote_per_year"] = 1.0
            entry["note"] = "Kannisto hazard tends to 1/yr; fitted only to 80-110 where SSA values are largely actuarial extrapolation above ~100."
        results["models"][sex][model] = entry
        xs = np.linspace(0, 119, 600)
        if model == "kannisto":
            xs = np.linspace(60, 119, 300)
            ys = mu_kannisto(xs, *p)
        elif model == "gamma_gompertz":
            xs = np.linspace(40, 119, 400)
            ys = mu_ggm(xs, p[0], p[1], 0.0, p[2])
        else:
            ys = mu_gm(xs, *p)
        ax.semilogy(xs, ys, style, lw=1.6, label=f"{model} (fit {lo}-{hi})")
        curves[f"{model}_{sex}"] = np.exp(pred_all)
    ax.set_title(f"US {sex}s, 2023 period hazard and fitted mortality laws")
    ax.set_xlabel("age (years)")
    ax.axvspan(95, 119, color="0.9", zorder=0)
    ax.text(97, 2e-5, "SSA values here are\nlargely extrapolated", fontsize=8)
    ax.grid(alpha=0.3, which="both")
axes[0].set_ylabel("one-year hazard  -ln(1-q_x)  (per year, log scale)")
axes[0].legend(fontsize=8, loc="upper left")
plt.tight_layout()
plt.savefig(f"{FIG}/fit_mortality_laws.png", dpi=150)
pd.DataFrame(curves).to_csv(f"{DATA}/fit_mortality_laws.csv", index=False)

# residual plot for Gompertz-Makeham over all ages
fig, ax = plt.subplots(figsize=(8, 4.5))
rows = {"age": ages_all}
for sex, col in [("male", "C0"), ("female", "C3")]:
    p = [results["models"][sex]["gompertz_makeham"]["params"][k] for k in "abc"]
    pred = mu_gm(ages_all, *p)
    pred_int = H_gm(ages_all + 1, *p) - H_gm(ages_all, *p)
    r = 100 * (d[f"mu_{sex}"].values / pred_int - 1)
    rows[f"pct_obs_minus_GM_{sex}"] = r
    ax.plot(ages_all, r, color=col, label=sex)
ax.axhline(0, color="k", lw=0.8)
ax.set_ylim(-100, 400)
ax.set_xlabel("age (years)"); ax.set_ylabel("observed / Gompertz-Makeham - 1  (%)")
ax.set_title("Where Gompertz-Makeham (fit 30-95) fails: infancy, young-adult hump, oldest ages")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(f"{FIG}/fit_residuals_gm.png", dpi=150)
pd.DataFrame(rows).to_csv(f"{DATA}/fit_residuals_gm.csv", index=False)

# observed life-table quantities
for sex in ["male", "female"]:
    q = d[f"{sex}_qx"].values
    e0_calc = lifetable_e0(q)
    lx = np.concatenate([[1], np.cumprod(1 - q)])
    dx = -np.diff(lx)
    mode_obs = int(np.argmax(dx[10:]) + 10)
    median = float(np.interp(0.5, lx[::-1], np.arange(len(lx))[::-1]))
    results.setdefault("observed", {})[sex] = {
        "e0_published": float(d.loc[0, f"{sex}_ex"]), "e0_recomputed_from_qx": e0_calc,
        "e30_published": float(d.loc[30, f"{sex}_ex"]), "e65_published": float(d.loc[65, f"{sex}_ex"]),
        "modal_age_at_death_adult": mode_obs, "median_age_at_death": median,
        "hazard_age20_per_year": float(d.loc[20, f"mu_{sex}"]),
        "hazard_age30_per_year": float(d.loc[30, f"mu_{sex}"]),
        "minimum_hazard_per_year": float(d[f"mu_{sex}"].min()),
        "age_of_minimum_hazard": int(d[f"mu_{sex}"].idxmin()),
    }
    # empirical doubling time from local log-slope 40-90
    sel = (ages_all >= 40) & (ages_all <= 90)
    slope = np.polyfit(ages_all[sel], np.log(d[f"mu_{sex}"].values[sel]), 1)[0]
    results["observed"][sex]["loglinear_slope_40_90"] = float(slope)
    results["observed"][sex]["doubling_time_40_90_years"] = float(np.log(2) / slope)
    m_ = d[f"mu_{sex}"].values
    results["observed"][sex]["local_log_slope_10y_centered"] = {
        str(a_): float(np.log(m_[a_ + 5] / m_[a_ - 5]) / 10) for a_ in range(30, 115, 10)}

# Strehler-Mildvan style convergence between sexes (Gompertz fits): x* where lines cross
gm_ = results["models"]
am, bm = gm_["male"]["gompertz"]["params"]["a"], gm_["male"]["gompertz"]["params"]["b"]
af, bf = gm_["female"]["gompertz"]["params"]["a"], gm_["female"]["gompertz"]["params"]["b"]
results["sex_convergence_age_gompertz"] = float((np.log(af) - np.log(am)) / (bm - bf))
pm, pf = gm_["male"]["gompertz_makeham"]["params"], gm_["female"]["gompertz_makeham"]["params"]
results["sex_convergence_age_gm_senescent_part"] = float((np.log(pf["a"]) - np.log(pm["a"])) / (pm["b"] - pf["b"]))

update_results("fits", results)
if __name__ == "__main__":
    import json
    for sex in ["male", "female"]:
        for m, e in results["models"][sex].items():
            print(sex, m, {k: round(v, 7) for k, v in e["params"].items()},
                  "CI", e["ci95"], "rmse", round(e["rmse_log_hazard_in_range"], 4))
            print("   ", {k: v for k, v in e.items() if k not in ("params", "ci95", "fit_ages")})
    print(json.dumps(results["observed"], indent=1))
    print("sex convergence age", results["sex_convergence_age_gompertz"])
