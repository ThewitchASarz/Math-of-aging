"""Part 1: converting a hazard ratio (or a slowed rate of ageing) into years of life.

All claims are checked numerically here. Parameters are fitted (see common.py).
Run:  python3 translation_math.py
"""
import numpy as np
import pandas as pd
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import special

from common import (MORTALITY, life_expectancy, le_gompertz_closed, update_results,
                    FIG, DATA)

out = {}
EG = float(np.euler_gamma)

# ---------------------------------------------------------------- 1. symbolic checks
x, a, b, r, t, M, s = sp.symbols("x a b r t M s", positive=True)
mu = a * sp.exp(b * x)
shift = sp.log(r) / b
# r*mu(x) == mu(x + ln(r)/b): a proportional hazard ratio is an exact age shift
ok_shift = sp.simplify(r * mu - mu.subs(x, x + shift)) == 0
# cumulative hazard from x0 over t:  H = mu0 (e^{bt}-1)/b ; check derivative
mu0 = sp.symbols("mu0", positive=True)
H = mu0 * (sp.exp(b * t) - 1) / b
ok_H = sp.simplify(sp.diff(H, t) - mu0 * sp.exp(b * t)) == 0
# series of e^M E1(M) for small M: -gamma - ln M + M(1 - gamma - ln M) + O(M^2 ln M)
ser = sp.series(sp.exp(M) * sp.expint(1, M), M, 0, 2).removeO()
target = -sp.EulerGamma - sp.log(M) + M * (1 - sp.EulerGamma - sp.log(M))
ok_series = sp.simplify(sp.expand(ser - target)) == 0
out["symbolic_checks"] = {
    "hazard_ratio_equals_age_shift": bool(ok_shift),
    "cumulative_hazard_derivative": bool(ok_H),
    "small_M_series_of_expM_E1M": bool(ok_series),
}
print("symbolic:", out["symbolic_checks"])


def gain_exact_gompertz(x0, a_, b_, hr):
    return le_gompertz_closed(x0, a_, b_, hr) - le_gompertz_closed(x0, a_, b_)


def gain_first_order(hr, b_):
    return -np.log(hr) / b_


def gain_second_order(x0, a_, b_, hr):
    Mv = a_ * np.exp(b_ * x0) / b_
    f = lambda m: -EG - np.log(m) + m * (1 - EG - np.log(m))
    return (f(hr * Mv) - f(Mv)) / b_


# ---------------------------------------------------------------- 2. closed form vs integral
checks = []
for spn, ages in (("human", [30, 50, 65, 80]), ("dog", [2, 5, 7, 10])):
    p = MORTALITY[spn]
    for x0 in ages:
        for hr in (0.95, 0.8, 0.5):
            e_num = life_expectancy(x0, p["a"], p["b"], 0.0, hr) - life_expectancy(x0, p["a"], p["b"], 0.0)
            e_cf = gain_exact_gompertz(x0, p["a"], p["b"], hr)
            e_mk = life_expectancy(x0, p["a"], p["b"], p["c"], hr) - life_expectancy(x0, p["a"], p["b"], p["c"])
            Mv = p["a"] * np.exp(p["b"] * x0) / p["b"]
            checks.append(dict(species=spn, start_age=x0, HR=hr, M=Mv,
                               gain_numeric=e_num, gain_closed_form=e_cf,
                               gain_first_order=gain_first_order(hr, p["b"]),
                               gain_second_order=gain_second_order(x0, p["a"], p["b"], hr),
                               gain_with_makeham=e_mk))
df = pd.DataFrame(checks)
df["rel_err_first_order"] = df.gain_first_order / df.gain_numeric - 1
df.to_csv(f"{DATA}/translation_checks.csv", index=False)
max_cf_err = float((df.gain_closed_form - df.gain_numeric).abs().max())
out["closed_form_max_abs_error_years"] = max_cf_err
print(df.round(3).to_string())
assert max_cf_err < 1e-4

# when does -ln(HR)/b hold?  relative error vs M (HR = 0.8)
Ms = np.geomspace(1e-4, 3, 200)
b0 = 1.0  # error depends only on M and HR when expressed relative to -ln(HR)/b
rel = []
for Mv in Ms:
    exact = (np.exp(0.8 * Mv) * special.exp1(0.8 * Mv) - np.exp(Mv) * special.exp1(Mv)) / b0
    rel.append(exact / (-np.log(0.8) / b0) - 1)
rel = np.array(rel)
M_10pct = float(Ms[np.argmax(np.abs(rel) > 0.10)])
M_25pct = float(Ms[np.argmax(np.abs(rel) > 0.25)])
out["first_order_error_HR0.8"] = {"M_at_10pct_error": M_10pct, "M_at_25pct_error": M_25pct}
pd.DataFrame({"M": Ms, "relative_error_of_minus_lnHR_over_b": rel}).to_csv(
    f"{DATA}/fig_approx_error.csv", index=False)

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.semilogx(Ms, 100 * rel, color="#2a6f97")
ax.axhline(-10, ls=":", color="grey")
ax.set_xlabel("M = Gompertz hazard at start age / b  (dimensionless)")
ax.set_ylabel("error of -ln(HR)/b  (%)")
ax.set_title("When the shortcut -ln(HR)/b works (HR = 0.8, pure Gompertz)")
for spn, ages in (("human", [50, 80]), ("dog", [7, 10])):
    p = MORTALITY[spn]
    for x0 in ages:
        Mv = p["a"] * np.exp(p["b"] * x0) / p["b"]
        ax.axvline(Mv, color="#c44" if spn == "dog" else "#4a4", alpha=.5)
        ax.text(Mv, -55, f"{spn} {x0}", rotation=90, fontsize=7, va="bottom")
ax.text(1.2e-4, 3, "fitted Gompertz parameters", fontsize=7)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_approx_error.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 3. gain vs HR
hrs = np.linspace(0.3, 1.0, 71)
rows = []
for spn, x0 in (("human", 50), ("human", 80), ("dog", 7), ("dog", 10)):
    p = MORTALITY[spn]
    base = life_expectancy(x0, p["a"], p["b"], p["c"])
    for hr in hrs:
        rows.append(dict(species=spn, start_age=x0, HR=hr,
                         gain_exact_GM=life_expectancy(x0, p["a"], p["b"], p["c"], hr) - base,
                         gain_approx=-np.log(hr) / p["b"]))
g = pd.DataFrame(rows)
g.to_csv(f"{DATA}/fig_gain_vs_hr.csv", index=False)
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, spn in zip(axs, ("human", "dog")):
    for x0, col in zip(sorted(g[g.species == spn].start_age.unique()), ("#2a6f97", "#c0392b")):
        d = g[(g.species == spn) & (g.start_age == x0)]
        ax.plot(d.HR, d.gain_exact_GM, color=col, label=f"exact, start age {x0}")
    d = g[(g.species == spn)]
    ax.plot(hrs, -np.log(hrs) / MORTALITY[spn]["b"], "k--", label="-ln(HR)/b")
    ax.set_xlabel("hazard ratio applied from start age"); ax.set_ylabel("years gained")
    ax.set_title(f"{spn} (b = {MORTALITY[spn]['b']}/yr, fitted)")
    ax.legend(fontsize=8); ax.invert_xaxis()
fig.tight_layout(); fig.savefig(f"{FIG}/fig_gain_vs_hr.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 4. slowing the rate of ageing
# From x0 on, the Gompertz slope becomes s*b.
# e(x0) = e^{M_s} E1(M_s)/(b s), M_s = mu0/(b s)  ~  [ln(b s/mu0) - gamma]/(b s)
slow = []
for spn, x0s in (("human", [0, 30, 50, 65, 80]), ("dog", [0, 3, 7, 10])):
    p = MORTALITY[spn]
    for x0 in x0s:
        mu0v = p["a"] * np.exp(p["b"] * x0)
        for sv in (0.95, 0.9, 0.8, 0.5):
            e_num = life_expectancy(x0, p["a"], p["b"], 0.0, 1.0, sv)
            e_cf = le_gompertz_closed(x0, p["a"], p["b"], 1.0, sv)
            e_ap = (np.log(p["b"] * sv / mu0v) - EG) / (p["b"] * sv)
            base = life_expectancy(x0, p["a"], p["b"], 0.0)
            gm = life_expectancy(x0, p["a"], p["b"], p["c"], 1.0, sv) - life_expectancy(x0, p["a"], p["b"], p["c"])
            slow.append(dict(species=spn, start_age=x0, s=sv, LE_base=base,
                             gain_numeric=e_num - base, gain_closed=e_cf - le_gompertz_closed(x0, p["a"], p["b"]),
                             gain_log_approx=e_ap - (np.log(p["b"] / mu0v) - EG) / p["b"],
                             gain_with_makeham=gm))
sl = pd.DataFrame(slow)
sl.to_csv(f"{DATA}/slowing_rate_gains.csv", index=False)
assert (sl.gain_closed - sl.gain_numeric).abs().max() < 1e-3
print(sl.round(2).to_string())

# ---------------------------------------------------------------- 5. gain vs start age
starts_h = np.arange(20, 91, 2); starts_d = np.arange(1, 14.01, 0.5)
rows = []
for spn, starts in (("human", starts_h), ("dog", starts_d)):
    p = MORTALITY[spn]
    for x0 in starts:
        base = life_expectancy(x0, p["a"], p["b"], p["c"])
        rows.append({"species": spn, "start_age": float(x0),
                     "gain_HR_0.8": life_expectancy(x0, p["a"], p["b"], p["c"], 0.8) - base,
                     "gain_slow_s_0.9": life_expectancy(x0, p["a"], p["b"], p["c"], 1.0, 0.9) - base})
st = pd.DataFrame(rows)
st.to_csv(f"{DATA}/fig_gain_vs_start_age.csv", index=False)
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, spn in zip(axs, ("human", "dog")):
    d = st[st.species == spn]
    ax.plot(d.start_age, d["gain_HR_0.8"], label="hazard x 0.8 from start age", color="#2a6f97")
    ax.plot(d.start_age, d["gain_slow_s_0.9"], label="ageing rate x 0.9 from start age", color="#c0392b")
    ax.axhline(-np.log(0.8) / MORTALITY[spn]["b"], ls="--", color="grey", label="-ln(0.8)/b")
    ax.set_xlabel("age treatment starts (years)"); ax.set_ylabel("remaining-life years gained")
    ax.set_title(f"{spn}: fitted Gompertz parameters"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_gain_vs_start_age.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 6. what 'solving ageing' needs
need = {}
for spn, x0 in (("human", 50), ("dog", 7)):
    p = MORTALITY[spn]
    base = life_expectancy(x0, p["a"], p["b"], p["c"])
    mu0v = p["a"] * np.exp(p["b"] * x0)
    # constant HR needed to add 10, 20, 50 years (first order) and exact (Makeham)
    from scipy.optimize import brentq
    tgt = {}
    for add in ((10, 20, 50) if spn == "human" else (1, 2, 5)):
        hr_ex = brentq(lambda h: life_expectancy(x0, p["a"], p["b"], p["c"], h) - base - add, 1e-6, 1.0)
        tgt[f"HR_to_add_{add}y"] = {"first_order": float(np.exp(-add * p["b"])), "exact_GM": float(hr_ex)}
    # stop ageing completely at x0 (hazard frozen at mu(x0)+c)
    frozen = 1.0 / (mu0v + p["c"])
    # slowing needed to add the same amounts
    for add in ((10, 20, 50) if spn == "human" else (1, 2, 5)):
        try:
            s_ex = brentq(lambda sv: life_expectancy(x0, p["a"], p["b"], p["c"], 1.0, sv) - base - add, 1e-3, 1.0)
        except ValueError:
            s_ex = None
        tgt[f"slope_multiplier_to_add_{add}y"] = s_ex
    need[spn] = {"start_age": x0, "baseline_remaining_LE": base,
                 "remaining_LE_if_ageing_frozen_at_start": frozen, **tgt}
out["what_solving_ageing_requires"] = need
print(need)

out["mortality_parameters"] = {k: {kk: float(vv) if not isinstance(vv, list) else vv for kk, vv in v.items()}
                                 for k, v in MORTALITY.items()}
out["example_gains"] = {
    "human_HR0.8_from50_first_order": float(-np.log(0.8) / MORTALITY["human"]["b"]),
    "human_HR0.8_from50_exact_GM": float(st[(st.species == "human") & (st.start_age == 50)]["gain_HR_0.8"].iloc[0]),
    "human_HR0.8_from80_exact_GM": float(st[(st.species == "human") & (st.start_age == 80)]["gain_HR_0.8"].iloc[0]),
    "dog_HR0.8_from7_exact_GM": float(st[(st.species == "dog") & (st.start_age == 7)]["gain_HR_0.8"].iloc[0]),
    "dog_HR0.8_from10_exact_GM": float(st[(st.species == "dog") & (st.start_age == 10)]["gain_HR_0.8"].iloc[0]),
    "human_s0.9_from50_exact_GM": float(st[(st.species == "human") & (st.start_age == 50)]["gain_slow_s_0.9"].iloc[0]),
    "human_s0.9_from20_exact_GM": float(st[(st.species == "human") & (st.start_age == 20)]["gain_slow_s_0.9"].iloc[0]),
    "dog_s0.9_from7_exact_GM": float(st[(st.species == "dog") & (st.start_age == 7)]["gain_slow_s_0.9"].iloc[0]),
}
print(out["example_gains"])
update_results("translation_math", out)
