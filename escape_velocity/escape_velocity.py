"""Parts 1-2: definitions of 'solved' and the escape-velocity threshold.

Model: Gompertz-Makeham hazard in calendar year t at age x
    mu(x, t) = a(t) * exp(b(t) * x) + c
Escape velocity conditions are derived in findings.md and checked numerically here.
"""
import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

from common import (B_DOG, cohort_e, B_HUMAN, C, DATA, INK2, MUTED, SRC, US2024_INJURY_AADR,
                    dog_params, finish, fit_us2024, gm_hazard, human_params, period_e,
                    plt, save_json)

np.random.seed(0)
out = {}

# ------------------------------------------------------------ 0. calibration
fit = fit_us2024()
out["us2024_gm_fit"] = fit
C_H = US2024_INJURY_AADR  # Makeham term for the escape-velocity calculations (see findings)
HUMAN = {b: human_params(b, C_H) for b in B_HUMAN}
DOG = {b: dog_params(b) for b in B_DOG}

# ------------------------------------------------------------ 1. symbolic check
# Claim: if a(t) = a0 exp(-rho t) with b, c fixed, then
#   de/dt = -(rho/b) de/dx   and   de/dx = mu e - 1.
x, t, s, a0, b, c, rho = sp.symbols("x t s a0 b c rho", positive=True)
z = a0 * sp.exp(-rho * t) * sp.exp(b * x)
H = c * s + z / b * (sp.exp(b * s) - 1)   # cumulative hazard from x to x+s
integrand = sp.exp(-H)
lhs = sp.diff(integrand, t)
rhs = -(rho / b) * sp.diff(integrand, x)
sym_ok = sp.simplify(lhs - rhs) == 0
out["symbolic_check_de_dt_equals_minus_rho_over_b_de_dx"] = bool(sym_ok)


def de_dt_numeric(x, p, rho, h=1e-3):
    """Finite-difference d e(x,t)/dt when a falls at exponential rate rho."""
    e1 = period_e(x, p["a"] * np.exp(-rho * h), p["b"], p["c"])
    e0 = period_e(x, p["a"], p["b"], p["c"])
    return (e1 - e0) / h


def de_dt_slope(x, p, sigma, h=1e-4):
    """d e(x,t)/dt when b falls by sigma per year (absolute), a fixed (age-0 anchor)."""
    e1 = period_e(x, p["a"], p["b"] - sigma * h, p["c"])
    e0 = period_e(x, p["a"], p["b"], p["c"])
    return (e1 - e0) / h


def requirements(x, p):
    mu = gm_hazard(x, **p)
    e = period_e(x, **p)
    de_dx = (period_e(x + 1e-3, **p) - e) / 1e-3
    rho_ind = p["b"]                       # individual: own remaining LE does not fall
    rho_strong = p["b"] / (1 - mu * e)     # fixed-age remaining LE rises >= 1 yr/yr
    # numeric verification of the strong threshold
    chk = de_dt_numeric(x, p, rho_strong)
    # slope route: de/dt = -sigma * de/db is linear in sigma (instantaneous rate), so
    # strong: sigma = 1/(-de/db); individual: sigma = -de_dx/(-de/db)
    neg_de_db = de_dt_slope(x, p, 1.0)   # = -de/db (with sigma=1)
    sig_strong = 1.0 / neg_de_db
    sig_ind = -de_dx / neg_de_db
    return dict(age=x, mu=mu, e_period=e, de_dx=de_dx, mu_e_minus_1=mu * e - 1,
                rho_individual=rho_ind, pct_mortality_fall_per_yr_individual=100 * (1 - np.exp(-rho_ind)),
                rho_strong=rho_strong, pct_mortality_fall_per_yr_strong=100 * (1 - np.exp(-rho_strong)),
                check_de_dt_at_rho_strong=chk,
                slope_fall_per_yr_strong=sig_strong, slope_fall_per_yr_individual=sig_ind,
                slope_fall_pct_of_b_strong=100 * sig_strong / p["b"],
                slope_fall_pct_of_b_individual=100 * sig_ind / p["b"])


rows = []
for b_, p in HUMAN.items():
    for age in (30, 50, 70):
        r = requirements(age, p)
        r.update(species="human", b=b_)
        rows.append(r)
for b_, p in DOG.items():
    for age in (2, 6, 10):
        r = requirements(age, p)
        r.update(species="dog", b=b_)
        rows.append(r)
req = pd.DataFrame(rows)
req.to_csv(f"{DATA}/ev_requirements.csv", index=False)
out["ev_requirements"] = req.to_dict(orient="records")
out["identity_check_max_abs_err_de_dx"] = float(np.max(np.abs(req.de_dx - req.mu_e_minus_1)))
out["strong_threshold_check_max_abs_err"] = float(np.max(np.abs(req.check_de_dt_at_rho_strong - 1)))

# ------------------------------------------------------------ 2. cohort view
ratios = np.round(np.linspace(0, 1.2, 49), 4)
coh = {"rho_over_b": ratios}
pick = {"human age 50 (b=0.085)": (50, HUMAN[0.085]), "human age 70 (b=0.085)": (70, HUMAN[0.085]),
        "dog age 6 (b=0.32)": (6, DOG[0.32])}
for lab, (x0, p) in pick.items():
    coh[lab] = [cohort_e(x0, p, rr * p["b"]) for rr in ratios]
coh = pd.DataFrame(coh)
coh.to_csv(f"{DATA}/cohort_remaining_le_vs_rho.csv", index=False)
out["cohort_remaining_le"] = {lab: {"rho=0": float(coh[lab].iloc[0]),
                                     "rho=0.5b": float(coh[lab][np.isclose(coh.rho_over_b, 0.5)].iloc[0]),
                                     "rho=0.9b": float(coh[lab][np.isclose(coh.rho_over_b, 0.9)].iloc[0]),
                                     "rho=b": float(coh[lab][np.isclose(coh.rho_over_b, 1.0)].iloc[0])}
                              for lab in pick}

# ------------------------------------------------------------ 3. history
hist = [
    # label, age, value_start, value_end, year_start, year_end, source
    ("Record female e0, 1840-2000 (Oeppen-Vaupel slope)", 0, None, None, 1840, 2000, SRC["oeppen_vaupel"], 0.243),
    ("US e0, 2010-2024", 0, 78.7, 79.0, 2010, 2024, SRC["hus2016"] + " ; " + SRC["cdc2024"], None),
    ("US e65, 1950-2010", 65, 13.9, 19.1, 1950, 2010, SRC["hus2016"], None),
    ("US e65, 2010-2024", 65, 19.1, 19.7, 2010, 2024, SRC["hus2016"] + " ; " + SRC["cdc2024"], None),
    ("Japan female e0, 2010-2024", 0, 86.30, 87.13, 2010, 2024, SRC["japan2024_hist"], None),
    ("Japan male e0, 2010-2024", 0, 79.55, 81.09, 2010, 2024, SRC["japan2024_hist"], None),
    ("Japan female e65, 2017-2024", 65, 24.43, 24.38, 2017, 2024, SRC["japan2017"] + " ; " + SRC["japan2024"], None),
    ("Japan female e90, 2017-2024", 90, 5.61, 5.55, 2017, 2024, SRC["japan2017"] + " ; " + SRC["japan2024"], None),
]
H = []
p = HUMAN[0.085]
for lab, age, v0, v1, y0, y1, src, slope in hist:
    rate = slope if slope is not None else (v1 - v0) / (y1 - y0)
    mu, e = gm_hazard(age, **p), period_e(age, **p)
    rho_implied = p["b"] * rate / (1 - mu * e)
    H.append(dict(series=lab, age=age, gain_yr_per_yr=rate, years=f"{y0}-{y1}",
                  gap_factor_to_1yr_per_yr=(1 / rate) if rate > 0 else float("inf"),
                  rho_implied_model=rho_implied,
                  pct_mortality_fall_per_yr_implied=100 * (1 - np.exp(-rho_implied)),
                  source=src))
H = pd.DataFrame(H)
H.to_csv(f"{DATA}/historical_gains_vs_required.csv", index=False)
out["historical"] = H.drop(columns="source").to_dict(orient="records")

# ------------------------------------------------------------ figures
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for i, (sp_, ages) in enumerate((("human", np.arange(20, 96, 2)), ("dog", np.arange(1, 16, 0.5)))):
    ax = axs[i]
    params = HUMAN[0.085] if sp_ == "human" else DOG[0.32]
    strong = [100 * (1 - np.exp(-params["b"] / (1 - gm_hazard(a_, **params) * period_e(a_, **params)))) for a_ in ages]
    ind = [100 * (1 - np.exp(-params["b"]))] * len(ages)
    ax.plot(ages, strong, color=C[0], label="remaining LE at fixed age rises 1 yr/yr")
    ax.plot(ages, ind, color=C[1], label="person's own remaining LE stops falling")
    pd.DataFrame({"age": ages, "pct_strong": strong, "pct_individual": ind}).to_csv(
        f"{DATA}/fig_required_mortality_fall_{sp_}.csv", index=False)
    ax.set_xlabel("age (years)")
    ax.set_ylabel("required fall in senescent death rate, % per calendar year")
    ax.set_title(f"{sp_}: b = {params['b']} /yr (fitted)")
    ax.set_ylim(0, 100 if sp_ == "dog" else 40)
axs[0].legend(loc="upper left", fontsize=8)
finish(fig, "fig_required_mortality_fall",
       "Gompertz-Makeham; human level calibrated to US 2024 rate at 70 (CDC db548); dog a fitted to median 12.5 y, Makeham term 0.01/yr ILLUSTRATIVE.")

fig, ax = plt.subplots(figsize=(7, 4))
for i, lab in enumerate(pick):
    ax.plot(coh.rho_over_b, coh[lab], color=C[i], label=lab)
ax.set_yscale("log")
ax.axvline(1.0, color=MUTED, lw=1, ls="--")
ax.text(1.01, 30, "rate = b\n(hazard stops rising)", color=INK2, fontsize=8)
ax.set_xlabel("rate at which death rates fall each year, as a fraction of Gompertz slope b")
ax.set_ylabel("expected remaining years (log scale)")
ax.set_title("Remaining life if improvement continues at a constant rate forever")
ax.legend(fontsize=8)
finish(fig, "fig_cohort_remaining_le", "Model calculation. Human background hazard 53.3/100k/yr (US 2024 injury AADR); dog background 0.01/yr ILLUSTRATIVE.")

fig, ax = plt.subplots(figsize=(8, 4))
hh = H.sort_values("gain_yr_per_yr")
cols = [C[0] if g > 0 else C[7] for g in hh.gain_yr_per_yr]
ax.barh(hh.series, hh.gain_yr_per_yr, color=cols, height=0.6)
ax.axvline(1.0, color=INK2, lw=1.5, ls="--")
ax.text(0.98, -0.6, "escape velocity\n(1 yr per yr)", ha="right", fontsize=8, color=INK2)
for yv, g in enumerate(hh.gain_yr_per_yr):
    ax.text(max(g, 0) + 0.01, yv, f"{g:+.3f}", va="center", fontsize=8, color=INK2)
ax.set_xlabel("observed gain in life expectancy, years per calendar year")
ax.set_xlim(-0.05, 1.1)
ax.set_title("Observed gains vs the 1-year-per-year escape-velocity line")
finish(fig, "fig_history_vs_required", "Sources: Oeppen & Vaupel 2002; CDC Health US 2016 Table 15; CDC db548; MHLW Japan life tables 2017, 2024.")

save_json("results_escape_velocity.json", out)
print(req[["species", "b", "age", "mu", "e_period", "pct_mortality_fall_per_yr_individual",
           "pct_mortality_fall_per_yr_strong", "check_de_dt_at_rho_strong",
           "slope_fall_pct_of_b_strong", "slope_fall_pct_of_b_individual"]].round(4).to_string())
print(H.drop(columns="source").round(3).to_string())
print(fit, "symbolic:", sym_ok)
print(coh.iloc[[0, 20, 36, 40, 44]].round(1).to_string())
