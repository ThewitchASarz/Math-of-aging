"""Part 4: consequences if aging were solved.

(a) Lotka renewal equation / Leslie matrix, female-only, 1-year steps.
(b) Life expectancy with only age-independent (background) risk.
(c) Dog equivalents.
Fertility schedules and TFR scenarios are ILLUSTRATIVE (world TFR 2.25 per UN WPP 2024 shown for reference).
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from common import (C, DATA, INK2, MUTED, US2024_INJURY_AADR, US2024_RATE_5_14, US2024_RATE_15_24,
                    C_DOG, cohort_e, dog_params, finish, gm_hazard, human_params, period_e, plt, save_json)

np.random.seed(0)
out = {}
P = human_params(0.085, US2024_INJURY_AADR)
MAXAGE = 1500          # last class is open-ended
ages = np.arange(MAXAGE + 1)
FEM = 100 / 205        # share of births that are girls (sex ratio 105:100, conventional value, unverified here)


def survival_probs(regime):
    """One-year survival p_x for x = 0..MAXAGE."""
    if regime == "current (Gompertz-Makeham, US 2024 level)":
        mu = gm_hazard(ages + 0.5, **P)
    elif regime == "negligible senescence after 25":
        mu = gm_hazard(np.minimum(ages + 0.5, 25.0), **P)
    else:
        raise ValueError(regime)
    return np.exp(-np.minimum(mu, 50))


def fertility(tfr, lo=15, hi=50):
    """ILLUSTRATIVE fertility schedule: triangular-ish (beta-shaped) on [lo, hi)."""
    m = np.zeros(MAXAGE + 1)
    x = np.arange(lo, hi)
    u = (x + 0.5 - lo) / (hi - lo)
    w = u ** 1.8 * (1 - u) ** 3.0
    m[lo:hi] = tfr * FEM * w / w.sum()   # daughters per woman per year
    return m


def lotka(p, m):
    l = np.concatenate([[1.0], np.cumprod(p[:-1])])  # l_x survival to exact age x
    lm = l * m
    nrr = lm.sum()
    f = lambda r: np.sum(np.exp(-r * (ages + 0.5)) * lm) - 1.0
    r = brentq(f, -0.2, 0.2)
    T = np.sum((ages + 0.5) * np.exp(-r * (ages + 0.5)) * lm)  # mean age of childbearing in stable pop
    e0 = l[:-1].sum() + l[-1] / (1 - p[-1]) - 0.5   # open-ended last class, mid-year correction
    return r, nrr, T, l, e0


def project(p, m, n0, years):
    n = n0.copy()
    tot = [n.sum()]
    for _ in range(years):
        births = np.sum(m * n)
        new = np.zeros_like(n)
        new[1:] = n[:-1] * p[:-1]
        new[-1] += n[-1] * p[-1]
        new[0] = births
        n = new
        tot.append(n.sum())
    return np.array(tot)


regimes = ["current (Gompertz-Makeham, US 2024 level)", "negligible senescence after 25"]
tfrs = [1.0, 1.5, 2.1, 2.25, 3.0]
rows, proj = [], {}
for reg in regimes:
    p = survival_probs(reg)
    for tfr in tfrs:
        r, nrr, T, l, e0 = lotka(p, fertility(tfr))
        rows.append(dict(regime=reg, TFR=tfr, NRR=nrr, r_per_yr=r, gen_time=T,
                         doubling_or_halving_yr=np.log(2) / abs(r) if abs(r) > 1e-9 else np.inf,
                         e0=float(e0), l45=l[45]))
# extended fertility window (speculative): same TFR spread over 15-75
p_ns = survival_probs(regimes[1])
for tfr in (1.5, 2.1):
    r, nrr, T, l, e0 = lotka(p_ns, fertility(tfr, 15, 75))
    rows.append(dict(regime="negligible senescence + fertility window 15-75 (SPECULATIVE)", TFR=tfr, NRR=nrr,
                     r_per_yr=r, gen_time=T, doubling_or_halving_yr=np.log(2) / abs(r), e0=float(e0), l45=l[45]))
lot = pd.DataFrame(rows)
lot.to_csv(f"{DATA}/lotka_growth_rates.csv", index=False)
out["lotka"] = lot.to_dict(orient="records")

# projection: start from stationary structure of current mortality (1 birth/yr), switch regime at t=0
p_cur = survival_probs(regimes[0])
l_cur = np.concatenate([[1.0], np.cumprod(p_cur[:-1])])
n0 = l_cur / l_cur.sum()
YEARS = 300
pr = {"year": np.arange(YEARS + 1)}
for reg in regimes:
    p = survival_probs(reg)
    for tfr in (1.5, 2.1, 3.0):
        pr[f"{reg[:22]} | TFR {tfr}"] = project(p, fertility(tfr), n0, YEARS)
pr = pd.DataFrame(pr)
pr.to_csv(f"{DATA}/fig_population_projection.csv", index=False)
out["projection_pop_multiple"] = {k: {"y100": float(pr[k].iloc[100]), "y300": float(pr[k].iloc[300])}
                                  for k in pr.columns if k != "year"}

# stationary population size per annual birth = e0 (renewal identity), check
out["stationary_size_equals_e0_check"] = {
    reg: dict(e0_lifetable=float(lotka(survival_probs(reg), fertility(2.1))[4]),
              stationary_pop_per_birth=float(project(survival_probs(reg), np.zeros(MAXAGE + 1),
                                                     np.r_[1.0, np.zeros(MAXAGE)], 6000).sum()))
    for reg in regimes}

# ---------------------------------------------------------- (b) background risk only
bg = []
for lab, c in [("US 2024 all-cause death rate, ages 5-14", US2024_RATE_5_14),
               ("US 2024 unintentional-injury death rate (age-adjusted)", US2024_INJURY_AADR),
               ("US 2024 all-cause death rate, ages 15-24", US2024_RATE_15_24),
               ("Model hazard frozen at age 25 (GM fit)", gm_hazard(25, **P))]:
    bg.append(dict(source=lab, hazard_per_yr=c, mean_life_yr=1 / c, median_yr=np.log(2) / c,
                   p_survive_500y=np.exp(-500 * c), p_survive_1000y=np.exp(-1000 * c)))
bg = pd.DataFrame(bg)
bg.to_csv(f"{DATA}/background_risk_lifespans.csv", index=False)
out["background_only"] = bg.to_dict(orient="records")

# ---------------------------------------------------------- (c) dogs
dog = []
for b in (0.30, 0.32, 0.34):
    dp = dog_params(b)
    hp = human_params(0.085, US2024_INJURY_AADR)
    for x0 in (2, 6, 10):
        dog.append(dict(b_dog=b, age=x0,
                        baseline_remaining=period_e(x0, **dp),
                        # same calendar rate as the human individual-LEV threshold (rho = 0.085/yr)
                        with_rho_0p085=cohort_e(x0, dp, 0.085),
                        # same proportional slowing as a 50% cut in human slope
                        slope_halved=period_e(x0, dp["a"] * np.exp(dp["b"] * x0) / np.exp(dp["b"] / 2 * x0), dp["b"] / 2, dp["c"]),
                        dog_LEV_rho=dp["b"],
                        negligible_senescence_from_now=1 / gm_hazard(x0, **dp),
                        background_only=1 / dp["c"]))
dog = pd.DataFrame(dog)
dog.to_csv(f"{DATA}/dog_consequences.csv", index=False)
out["dog"] = dog.to_dict(orient="records")
# Note on slope_halved: hazard at current age kept, future rise at half the slope.

out["economics_scott_2021"] = {
    "value_1yr_LE_gain_via_slowed_aging_USD": 38e12, "value_10yr_USD": 367e12,
    "source": "https://www.london.edu/faculty-and-research/academic-research/a/alls-well-that-ages-well-the-economic-value-of-targeting-aging"}

# ---------------------------------------------------------- figures
fig, ax = plt.subplots(figsize=(8, 4.2))
cols = {1.5: C[0], 2.1: C[1], 3.0: C[2]}
for k in pr.columns[1:]:
    tfr = float(k.split("TFR ")[1])
    ls = "-" if k.startswith("negligible") else "--"
    ax.plot(pr.year, pr[k], color=cols[tfr], ls=ls, lw=1.8)
    ax.text(YEARS + 3, pr[k].iloc[-1], f"TFR {tfr}{' NS' if ls == '-' else ''}", fontsize=7, color=INK2, va="center")
ax.set_yscale("log"); ax.set_xlim(0, YEARS + 40)
ax.set_xlabel("years after aging stops (solid) or not (dashed)")
ax.set_ylabel("population relative to start (log)")
ax.set_title("Fertility, not death rate, sets the long-run direction")
finish(fig, "fig_population_projection", "Leslie projection, female-only; fertility schedule ILLUSTRATIVE; NS = negligible senescence after age 25.")

save_json("results_population.json", out)
print(lot.round(4).to_string())
print(bg.round(4).to_string())
print(dog.round(2).to_string())
print({k: v for k, v in out["projection_pop_multiple"].items()})
print(out["stationary_size_equals_e0_check"])
