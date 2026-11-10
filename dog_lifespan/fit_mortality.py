"""Part 1: fit Gompertz and Gompertz-Makeham models to dog and human life tables.

Data (transcribed, see CSV headers for URLs):
  dogs  : Teng et al. 2022 (UK VetCompass, deaths-only life table with counts)
          Montoya et al. 2023 (US Banfield, 13.3 M dogs, qx only)
  humans: US SSA 2023 period life table

Method: each 1-year interval [x, x+1) has model death probability
    q_x(theta) = 1 - exp(-H(x, x+1; theta)).
Teng: binomial maximum likelihood with d_x deaths out of l_x at risk.
Montoya / SSA: no usable counts, so we minimise squared error of log interval
cumulative hazard  log(-log(1-q_x))  (each age weighted equally).
Open-ended last intervals are excluded.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from common import (read, cumhaz, hazard, surv, life_expectancy, median_age,
                    update_results, save_csv, FIG, LN2)

rng = np.random.default_rng(20261110)

# ---------------------------------------------------------------- data checks
teng = read("teng2022_uk_dog_life_table.csv")
mont = read("montoya2023_us_dog_life_table.csv")
ssa = read("ssa2023_us_human_life_table.csv")

# Teng internal consistency: l_{x+1} = l_x - d_x and q_x = d_x / l_x (to 3 dp)
assert np.all(teng.lx.values[1:] == (teng.lx - teng.dx).values[:-1]), "Teng lx/dx mismatch"
assert np.allclose((teng.dx / teng.lx).round(3), teng.qx, atol=0.0015), "Teng qx mismatch"
# SSA internal consistency: l_{x+1} ~= l_x (1 - q_x)
for sex in ("male", "female"):
    pred = ssa[f"{sex}_lx"].values[:-1] * (1 - ssa[f"{sex}_qx"].values[:-1])
    obs = ssa[f"{sex}_lx"].values[1:]
    ok = np.abs(pred - obs) <= np.maximum(2.0, 1e-4 * obs)
    assert ok.all(), f"SSA {sex} lx inconsistency at ages {ssa.age.values[1:][~ok]}"
# Montoya: check e_x recursion roughly holds (e_x ~ (1-q_x)(1+e_{x+1}) + q_x*0.5)
m_e = (1 - mont.qx.values[:-1]) * (1 + mont.ex.values[1:]) + mont.qx.values[:-1] * 0.5
mont_e_err = np.max(np.abs(m_e - mont.ex.values[:-1]))
print(f"Montoya e_x recursion max abs error = {mont_e_err:.3f} y (transcription check)")
assert mont_e_err < 0.05


# ---------------------------------------------------------------- fitting
def unpack(p, makeham):
    a, b = np.exp(p[0]), np.exp(p[1])
    c = np.exp(p[2]) if makeham else 0.0
    return a, b, c


def q_model(x, a, b, c):
    return 1 - np.exp(-cumhaz(x, x + 1, a, b, c))


def fit_binomial(x, d, n, makeham):
    def nll(p):
        a, b, c = unpack(p, makeham)
        q = np.clip(q_model(x, a, b, c), 1e-12, 1 - 1e-12)
        return -np.sum(d * np.log(q) + (n - d) * np.log(1 - q))
    best = None
    for b0 in (0.1, 0.2, 0.3, 0.4):
        p0 = [np.log(0.005), np.log(b0)] + ([np.log(0.005)] if makeham else [])
        r = minimize(nll, p0, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-10, fatol=1e-10))
        if best is None or r.fun < best.fun:
            best = r
    return unpack(best.x, makeham), best.fun


def fit_logH(x, q, makeham):
    y = np.log(-np.log(1 - q))
    def sse(p):
        a, b, c = unpack(p, makeham)
        return np.sum((np.log(cumhaz(x, x + 1, a, b, c)) - y) ** 2)
    best = None
    for b0 in (0.05, 0.1, 0.2, 0.3, 0.4):
        for c0 in ((1e-4, 1e-3, 1e-2) if makeham else (None,)):
            p0 = [np.log(1e-3), np.log(b0)] + ([np.log(c0)] if makeham else [])
            r = minimize(sse, p0, method="Nelder-Mead", options=dict(maxiter=40000, xatol=1e-11, fatol=1e-12))
            if best is None or r.fun < best.fun:
                best = r
    a, b, c = unpack(best.x, makeham)
    rmse = np.sqrt(best.fun / len(x))
    return (a, b, c), rmse


def describe(a, b, c, x_e=2.0, table_e=None):
    return dict(a_per_year=a, b_per_year=b, c_per_year=c,
                MRDT_years=LN2 / b,
                hazard_at_age_2=hazard(2.0, a, b, c) if b > 0.2 else None,
                model_e_at=x_e, model_remaining_LE=life_expectancy(a, b, c, x0=x_e),
                table_remaining_LE=table_e)


results = {}
fits = {}

# --- Teng (UK) dogs, adult ages 2..19 for Gompertz, 1..19 for G-M
t = teng[teng.age_start <= 19]
for name, lo, mk in (("teng_gompertz_2to19", 2, False), ("teng_makeham_1to19", 1, True)):
    s = t[t.age_start >= lo]
    (a, b, c), nll = fit_binomial(s.age_start.values.astype(float), s.dx.values, s.lx.values, mk)
    fits[name] = (a, b, c)
    results[name] = describe(a, b, c, 2.0, float(teng.loc[teng.age_start == 2, "ex"].iloc[0]))
    results[name]["neg_log_lik"] = nll

# parametric bootstrap (resample deaths multinomially) for Teng Gompertz 2..19
s = teng[(teng.age_start >= 2)]
p_death = s.dx.values / s.dx.sum()
boot_b, boot_a = [], []
for _ in range(300):
    d = rng.multinomial(s.dx.sum(), p_death)
    lx = d[::-1].cumsum()[::-1]
    keep = s.age_start.values <= 19
    (a, b, c), _ = fit_binomial(s.age_start.values[keep].astype(float), d[keep], lx[keep], False)
    boot_a.append(a); boot_b.append(b)
boot_b = np.array(boot_b)
results["teng_gompertz_2to19"]["b_95CI_bootstrap"] = list(np.percentile(boot_b, [2.5, 97.5]))
results["teng_gompertz_2to19"]["MRDT_95CI_bootstrap"] = list(np.percentile(LN2 / boot_b, [2.5, 97.5]))

# --- Montoya (US) dogs, ages 2..16 Gompertz, 1..16 G-M
m = mont[mont.age_start <= 16]
for name, lo, mk in (("montoya_gompertz_2to16", 2, False), ("montoya_makeham_1to16", 1, True)):
    s = m[m.age_start >= lo]
    (a, b, c), rmse = fit_logH(s.age_start.values.astype(float), s.qx.values, mk)
    fits[name] = (a, b, c)
    results[name] = describe(a, b, c, 2.0, float(mont.loc[mont.age_start == 2, "ex"].iloc[0]))
    results[name]["rmse_log_interval_hazard"] = rmse

# --- Humans, SSA 2023, ages 40..94 Gompertz, 30..94 G-M; male, female, and 50:50 combined
lx_both = 0.5 * ssa.male_lx + 0.5 * ssa.female_lx   # ASSUMPTION: equal numbers of male and female births
q_both = 1 - lx_both.shift(-1) / lx_both
ssa = ssa.assign(both_qx=q_both)
ex_both_40 = None
for sex in ("male", "female", "both"):
    for name, lo, mk in ((f"human_{sex}_gompertz_40to94", 40, False), (f"human_{sex}_makeham_30to94", 30, True)):
        s = ssa[(ssa.age >= lo) & (ssa.age <= 94)]
        (a, b, c), rmse = fit_logH(s.age.values.astype(float), s[f"{sex}_qx"].values, mk)
        fits[name] = (a, b, c)
        d = describe(a, b, c, 40.0, float(ssa.loc[ssa.age == 40, f"{sex}_ex"].iloc[0]) if sex != "both" else None)
        d["hazard_at_age_2"] = None
        d["hazard_at_age_40"] = hazard(40.0, a, b, c)
        d["rmse_log_interval_hazard"] = rmse
        results[name] = d

for k, v in results.items():
    print(f"{k:32s} a={v['a_per_year']:.3e} b={v['b_per_year']:.4f} c={v['c_per_year']:.2e} "
          f"MRDT={v['MRDT_years']:.2f}y  model e({v['model_e_at']:.0f})={v['model_remaining_LE']:.2f} table={v['table_remaining_LE']}")
print("Teng b bootstrap 95% CI", results["teng_gompertz_2to19"]["b_95CI_bootstrap"])

# ---------------------------------------------------------------- headline comparison
dog_b = results["teng_gompertz_2to19"]["b_per_year"]
dog_b_us = results["montoya_gompertz_2to16"]["b_per_year"]
hum_b = results["human_both_gompertz_40to94"]["b_per_year"]
comparison = dict(
    dog_MRDT_UK_years=LN2 / dog_b, dog_MRDT_US_years=LN2 / dog_b_us,
    human_MRDT_US2023_years=LN2 / hum_b,
    ratio_b_dog_over_human_UK=dog_b / hum_b, ratio_b_dog_over_human_US=dog_b_us / hum_b,
    ratio_e0_human_over_dog=None,
    note=("Dog Gompertz slope b is ~3.5-4x the human slope, but dogs live ~6x shorter: the rest of the "
          "difference sits in the level parameter a (dogs start adult life at a far higher hazard)."))
comparison["ratio_e0_human_over_dog"] = float((0.5 * 75.79 + 0.5 * 81.06) / 12.69)
results["comparison"] = comparison
print(comparison)

# ---------------------------------------------------------------- size classes (calibration, not a fit)
size = read("montoya2023_size_e0.csv")
mc = read("mcmillan2024_size_face.csv")
a0, b0, c0 = fits["montoya_makeham_1to16"]
e0_model_overall = life_expectancy(a0, b0, c0, 0.0)
offset = e0_model_overall - 12.69   # model error at birth (juvenile mortality not modelled)
from scipy.optimize import brentq
size_rows = []
for _, r in size.iterrows():
    target = r.e0 + offset
    bA = brentq(lambda bb: life_expectancy(a0, bb, c0, 0.0) - target, 0.05, 1.5)   # scenario A: aging rate varies
    aB = brentq(lambda aa: life_expectancy(aa, b0, c0, 0.0) - target, 1e-7, 1.0)   # scenario B: level varies
    size_rows.append(dict(size=r["size"], kg_lo=r.kg_lo, kg_hi=r.kg_hi, e0_observed=r.e0,
                          A_b=bA, A_MRDT=LN2 / bA, A_a=a0, A_median=median_age(a0, bA, c0),
                          B_a=aB, B_b=b0, B_MRDT=LN2 / b0, B_median=median_age(aB, b0, c0), c=c0))
size_df = pd.DataFrame(size_rows)
print(size_df.round(4).to_string())
# Check against McMillan: proportional-hazards HR small vs large implied by scenario B
hr_B_small_large = size_df.set_index("size").loc["small", "B_a"] / size_df.set_index("size").loc["large", "B_a"]
results["size_class_calibration"] = dict(
    method=("Montoya US Gompertz-Makeham (ages 1-16) parameters held fixed except one; that one parameter is "
            "solved so model e0 matches each size group's published e0 (after removing the overall model-vs-table "
            "e0 offset of %.3f y). Scenario A: Gompertz slope b varies (size changes the rate of aging). "
            "Scenario B: level a varies (size multiplies hazard at all ages)." % offset),
    model_e0_offset_years=offset,
    rows=size_df.to_dict(orient="records"),
    check_McMillan=dict(observed_HR_small_vs_large=0.83, scenarioB_HR_small_vs_large=hr_B_small_large,
                        observed_median_gap_large_minus_small_years=12.7 - 11.9),
    note="Calibration of one parameter to one number per group; not an independent fit. Kraus et al. 2013 found size acts mainly through the rate of aging (scenario A).")
save_csv(size_df, "size_class_gompertz_calibration.csv",
         ["Derived by fit_mortality.py: one-parameter calibration of Gompertz-Makeham to Montoya 2023 e0 by size group",
          "Source e0: https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2023.1082102/pdf"])

# ---------------------------------------------------------------- figures
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
x = teng.age_start.values[:-1] + 0.5
ax[0].semilogy(x, -np.log(1 - teng.qx.values[:-1]), "o", label="UK VetCompass (Teng 2022)")
xm = mont.age_start.values[:-1] + 0.5
ax[0].semilogy(xm, -np.log(1 - mont.qx.values[:-1]), "s", mfc="none", label="US Banfield (Montoya 2023)")
xx = np.linspace(0, 20, 200)
a, b, c = fits["teng_gompertz_2to19"]; ax[0].semilogy(xx, hazard(xx, a, b, c), "-", label=f"Gompertz UK, b={b:.3f}/y")
a, b, c = fits["montoya_makeham_1to16"]; ax[0].semilogy(xx, hazard(xx, a, b, c), "--", label=f"G-Makeham US, b={b:.3f}/y")
ax[0].set_xlabel("dog age (years)"); ax[0].set_ylabel("annual hazard  -ln(1-q)"); ax[0].set_title("Dogs: mortality rises exponentially after ~age 3")
ax[0].legend(fontsize=8)
for sex, mk in (("male", "^"), ("female", "v")):
    ax[1].semilogy(ssa.age[:-1] + 0.5, -np.log(1 - ssa[f"{sex}_qx"][:-1]), mk, ms=3, label=f"US {sex} 2023 (SSA)")
xh = np.linspace(0, 110, 300)
a, b, c = fits["human_both_gompertz_40to94"]; ax[1].semilogy(xh, hazard(xh, a, b, c), "-", label=f"Gompertz 40-94, b={b:.3f}/y")
ax[1].set_xlabel("human age (years)"); ax[1].set_title("Humans (for comparison)"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{FIG}/dog_human_hazard_fits.png", dpi=150); plt.close()

rows = []
for src, df, xs in (("teng_uk", teng[:-1], x), ("montoya_us", mont[:-1], xm)):
    for xi, q in zip(xs, df.qx.values):
        rows.append(dict(series=src, age_mid=xi, observed_hazard=-np.log(1 - q)))
for nm in ("teng_gompertz_2to19", "montoya_makeham_1to16"):
    a, b, c = fits[nm]
    for xi in xx[::5]:
        rows.append(dict(series=nm, age_mid=xi, observed_hazard=np.nan, model_hazard=hazard(xi, a, b, c)))
save_csv(pd.DataFrame(rows), "dog_human_hazard_fits.csv", ["Series behind figures/dog_human_hazard_fits.png (dog panel); human panel uses ssa2023_us_human_life_table.csv"])

# Survival curves with time scaled by each species' doubling time
fig, ax = plt.subplots(figsize=(6.5, 4.5))
S_teng = teng.lx / teng.lx.iloc[0]
ax.plot(teng.age_start, S_teng, "o-", label="UK dogs (Teng 2022)")
S_m = np.concatenate([[1], np.cumprod(1 - mont.qx.values)])[:-1]
ax.plot(mont.age_start, S_m, "s-", mfc="none", label="US dogs (Montoya 2023)")
ax.plot(ssa.age / 6.3, lx_both / 1e5, "-", label="US humans 2023, age / 6.3")
ax.set_xlabel("dog age (years); human age divided by 6.3"); ax.set_ylabel("fraction surviving")
ax.set_xlim(0, 21); ax.legend(); ax.set_title("Survival curves on a common time scale")
plt.tight_layout(); plt.savefig(f"{FIG}/survival_curves.png", dpi=150); plt.close()
save_csv(pd.DataFrame(dict(age=teng.age_start, S_teng=S_teng.values)), "survival_curves_teng.csv", ["figure survival_curves.png"])
save_csv(pd.DataFrame(dict(age=mont.age_start, S_montoya=S_m)), "survival_curves_montoya.csv", ["figure survival_curves.png"])

# ---------------------------------------------------------------- store params for other streams
reuse = {}
for k, (a, b, c) in fits.items():
    reuse[k] = dict(a=a, b=b, c=c, units="a, c in deaths per year (hazard); b per year; age x in years; mu(x)=c+a*exp(b*x)")
results["parameters_for_reuse"] = reuse
a, b, c = fits["teng_gompertz_2to19"]
a2, b2, c2 = fits["montoya_makeham_1to16"]
results["RECOMMENDED_dog_parameters"] = {
    "headline_model": "Gompertz mu(x)=a*exp(b*x), adult dogs ages 2-19, UK VetCompass deaths 2016-2020 (Teng et al. 2022), binomial MLE",
    "a_per_year": a, "b_per_year": b, "MRDT_years": LN2 / b,
    "MRDT_95CI_sampling_only": results["teng_gompertz_2to19"]["MRDT_95CI_bootstrap"],
    "MRDT_range_across_datasets_and_models_years": [min(LN2 / v[1] for k, v in fits.items() if not k.startswith("human")),
                                                     max(LN2 / v[1] for k, v in fits.items() if not k.startswith("human"))],
    "life_expectancy_at_birth_years_UK_table": 11.23, "life_expectancy_at_birth_years_US_table": 12.69,
    "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9050668",
    "alternative_US_model": {"form": "Gompertz-Makeham mu(x)=c+a*exp(b*x), ages 1-16, US Banfield 2013-2019 (Montoya et al. 2023)",
                             "a_per_year": a2, "b_per_year": b2, "c_per_year": c2, "MRDT_years": LN2 / b2,
                             "source_url": "https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2023.1082102/pdf"},
    "human_comparator": {"form": "Gompertz ages 40-94, US SSA 2023 period table, 50:50 sexes",
                         "a_per_year": fits["human_both_gompertz_40to94"][0], "b_per_year": fits["human_both_gompertz_40to94"][1],
                         "MRDT_years": LN2 / fits["human_both_gompertz_40to94"][1], "source_url": "https://www.ssa.gov/oact/STATS/table4c6.html"},
    "units": "age in years; a and c are annual hazards (per year); b per year; MRDT = ln2/b in years",
    "caveat": "UK table is built from recorded deaths only (no censoring) and misses early puppy deaths; US table uses clinic survivorship (Sullivan method). Both are clinic populations."}
results["sources"] = {
    "teng2022": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9050668",
    "montoya2023": "https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2023.1082102/pdf",
    "ssa2023": "https://www.ssa.gov/oact/STATS/table4c6.html",
    "mcmillan2024": "https://www.nature.com/articles/s41598-023-50458-w.pdf"}
update_results("part1_mortality_fits", results)
print("done")
