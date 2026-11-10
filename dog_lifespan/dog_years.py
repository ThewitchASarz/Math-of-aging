"""Part 3: dog years to human years.

Mappings h(d) from dog age d (years) to human age h (years):
  R7      h = 7 d                                   (folk rule)
  AVMA    h = 15 d (d<=1); 15 + 9(d-1) (1<d<=2); 24 + 5(d-2) (d>2)   (piecewise vet rule, via AKC page)
  WANG    h = 16 ln d + 31                          (Wang et al. 2020, Labrador methylation clock)
  SQ      S_dog(d) = S_human(h)                     (equal fraction of the birth cohort still alive)
  EH      mu_dog(d) = mu_human(h)                   (equal annual risk of death)
"""
import numpy as np
import pandas as pd
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from common import read, load_results, update_results, save_csv, surv, hazard, FIG

R = load_results()["part1_mortality_fits"]
P = R["parameters_for_reuse"]
teng = read("teng2022_uk_dog_life_table.csv")
mont = read("montoya2023_us_dog_life_table.csv")
ssa = read("ssa2023_us_human_life_table.csv")
lx_h = (0.5 * ssa.male_lx + 0.5 * ssa.female_lx).values / 1e5   # ASSUMPTION: equal sexes
age_h = ssa.age.values.astype(float)

# ---------------------------------------------------------------- symbolic checks
d, h, ad, bd, ah, bh = sp.symbols("d h a_d b_d a_h b_h", positive=True)
# survival-quantile mapping under pure Gompertz
H_d = ad / bd * (sp.exp(bd * d) - 1)
h_sq = sp.log(1 + (ad * bh) / (ah * bd) * (sp.exp(bd * d) - 1)) / bh
assert sp.simplify(ah / bh * (sp.exp(bh * h_sq) - 1) - H_d) == 0
slope_limit = sp.limit(sp.diff(h_sq, d), d, sp.oo)
assert sp.simplify(slope_limit - bd / bh) == 0
# equal-hazard mapping under pure Gompertz
h_eh = sp.solve(sp.Eq(ad * sp.exp(bd * d), ah * sp.exp(bh * h)), h)[0]
assert sp.simplify(sp.diff(h_eh, d) - bd / bh) == 0
# Wang slope 16/d equals mortality slope b_d/b_h at d* = 16 b_h / b_d
print("SQ mapping:", h_sq, "\nasymptotic slope:", slope_limit, "\nEH mapping:", sp.simplify(h_eh))


# ---------------------------------------------------------------- mapping functions
def r7(x): return 7 * x
def avma(x):
    x = np.asarray(x, float)
    return np.where(x <= 1, 15 * x, np.where(x <= 2, 15 + 9 * (x - 1), 24 + 5 * (x - 2)))
def wang(x): return 16 * np.log(x) + 31


def interp_logS(age, S):
    lS = np.log(np.clip(S, 1e-12, 1))
    return lambda x: np.interp(x, age, lS)


S_teng = interp_logS(teng.age_start.values.astype(float), teng.lx.values / teng.lx.values[0])
S_mont_vals = np.concatenate([[1.0], np.cumprod(1 - mont.qx.values)])
S_mont = interp_logS(np.arange(len(S_mont_vals), dtype=float), S_mont_vals)
lS_h = np.log(np.clip(lx_h, 1e-12, 1))


def human_age_with_logS(target):
    # human log-survival is decreasing in age; invert by interpolation
    if target > lS_h[0] or target < lS_h[-2]:
        return np.nan
    return float(np.interp(-target, -lS_h, age_h))


def sq_empirical(x, Sfun):
    return np.array([human_age_with_logS(float(Sfun(xi))) for xi in np.atleast_1d(x)])


def sq_model(x, pd_, ph=None, x_join=2.0):
    """Survival-quantile mapping with a hybrid dog survival curve: empirical Montoya survival up to
    age x_join (juvenile deaths are not in the adult model), fitted/calibrated adult model after it.
    Human side uses the empirical SSA 2023 curve."""
    out = []
    for xi in np.atleast_1d(x):
        if xi <= x_join:
            logS = float(S_mont(xi))
        else:
            logS = float(S_mont(x_join)) + np.log(surv(xi, pd_["a"], pd_["b"], pd_["c"], x0=x_join))
        out.append(human_age_with_logS(logS))
    return np.array(out)


# human annual hazard from table (ages >= 10, where it is increasing)
mu_h_tab = -np.log(1 - (0.5 * ssa.male_qx + 0.5 * ssa.female_qx).values)   # approx. 50:50 hazard
sel = age_h >= 10


def eh_empirical(x, dog_age, dog_q):
    mu_d = np.interp(x, dog_age + 0.5, -np.log(1 - dog_q))
    return np.interp(np.log(mu_d), np.log(mu_h_tab[sel]), age_h[sel] + 0.5, left=np.nan, right=np.nan)


def eh_model(x, pd_, ph):
    # requires mu_d >= mu_h(0)
    mu_d = hazard(x, pd_["a"], pd_["b"], pd_["c"])
    res = []
    for m in np.atleast_1d(mu_d):
        f = lambda hh: hazard(hh, ph["a"], ph["b"], ph["c"]) - m
        res.append(brentq(f, 0, 140) if f(0) < 0 < f(140) else np.nan)
    return np.array(res)


hum = P["human_both_makeham_30to94"]
dogM = P["montoya_makeham_1to16"]
dogT = P["teng_gompertz_2to19"]
cal = {r["size"]: r for r in R["size_class_calibration"]["rows"]}
small = dict(a=cal["small"]["A_a"], b=cal["small"]["A_b"], c=cal["small"]["c"])
giant = dict(a=cal["giant"]["A_a"], b=cal["giant"]["A_b"], c=cal["giant"]["c"])

ages = np.array([0.5, 1, 2, 3, 4, 5, 7, 10, 12, 14, 16, 18])
tab = pd.DataFrame(dict(dog_age=ages, times7=r7(ages), AVMA=avma(ages), Wang2020=wang(ages),
                        SQ_US_empirical=sq_empirical(ages, S_mont),
                        SQ_UKdog_USHuman_empirical=sq_empirical(ages, S_teng),
                        SQ_US_model=sq_model(ages, dogM, hum),
                        SQ_small_model=sq_model(ages, small, hum),
                        SQ_giant_model=sq_model(ages, giant, hum),
                        EH_US_empirical=eh_empirical(ages, mont.age_start.values[:-1].astype(float), mont.qx.values[:-1]),
                        EH_US_model=eh_model(ages, dogM, hum)))
pd.set_option("display.width", 250)
print(tab.round(1).to_string(index=False))

# slopes: years of human age per dog year in adulthood (ages 4..12)
def slope(col):
    s = tab[(tab.dog_age >= 4) & (tab.dog_age <= 12)][["dog_age", col]].dropna()
    return float(np.polyfit(s.dog_age, s[col], 1)[0])
slopes = {c: slope(c) for c in tab.columns[1:]}
bd_v, bh_v = dogM["b"], hum["b"]
crossover = 16 * bh_v / bd_v
print("adult slopes:", {k: round(v, 2) for k, v in slopes.items()}, "crossover age", crossover)

disagreements = dict(
    at_dog_age_1={c: float(tab.loc[tab.dog_age == 1, c].iloc[0]) for c in tab.columns[1:]},
    at_dog_age_12={c: float(tab.loc[tab.dog_age == 12, c].iloc[0]) for c in tab.columns[1:]},
    at_dog_age_16={c: float(tab.loc[tab.dog_age == 16, c].iloc[0]) for c in tab.columns[1:]})
update_results("part3_dog_years", dict(
    table=tab.to_dict(orient="records"),
    adult_slope_human_years_per_dog_year_ages4to12=slopes,
    gompertz_asymptotic_slope_b_dog_over_b_human=bd_v / bh_v,
    wang_slope_equals_mortality_slope_at_dog_age=crossover,
    wang_formula=dict(formula="human_age = 16*ln(dog_age) + 31", log="natural log", breed="Labrador retrievers (104 dogs)",
                      humans=320, source="https://www.eurekalert.org/news-releases/464089 ; https://www.biorxiv.org/content/10.1101/829192.full.pdf"),
    human_reference="US SSA 2023 period table, 50:50 sex mix (assumption)",
    disagreements=disagreements))
save_csv(tab, "dog_years_mappings.csv", ["Series for figures/dog_years_mappings.png; produced by dog_years.py",
                                          "SQ = survival-quantile matching; EH = equal annual hazard; US = Montoya 2023 dogs vs SSA 2023 humans"])

# figure
xx = np.linspace(0.1, 18, 200)
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(xx, r7(xx), label="×7 folk rule")
ax.plot(xx, avma(xx), label="AVMA piecewise (15, +9, +5/y)")
ax.plot(xx, wang(xx), label="Wang 2020 epigenetic: 16 ln(d)+31 (Labradors)")
ax.plot(xx, sq_empirical(xx, S_mont), lw=2.5, label="same % still alive (US dogs vs US humans, data)")
ax.plot(xx, sq_model(xx, small, hum), "--", label="same % alive, small dogs (calibrated model)")
ax.plot(xx, sq_model(xx, giant, hum), "--", label="same % alive, giant dogs (calibrated model)")
ax.plot(xx, eh_empirical(xx, mont.age_start.values[:-1].astype(float), mont.qx.values[:-1]), ":", lw=2.5, label="same annual death risk (data)")
ax.set_xlabel("dog age (years)"); ax.set_ylabel("'equivalent' human age (years)"); ax.set_ylim(0, 125)
ax.set_title("Five ways to convert dog years, and where they disagree"); ax.legend(fontsize=7.5); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(f"{FIG}/dog_years_mappings.png", dpi=150); plt.close()
print("done")
