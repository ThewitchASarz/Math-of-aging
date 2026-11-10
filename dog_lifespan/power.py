"""Part 5: dogs as a model for human aging - how big and how long must a lifespan trial be?

Design (ASSUMPTIONS, stated once here and repeated in findings.md):
  * two arms, 1:1, everyone enrolled at the same age x0 at the same time, followed for T years;
  * no dropout, no competing loss to follow-up; endpoint = death from any cause;
  * treatment multiplies the whole hazard by HR from enrollment onward (proportional hazards);
  * two-sided alpha = 0.05, power = 0.80.
Number of deaths needed (Schoenfeld 1983, via the CASRAI guide):
    D = (z_{1-alpha/2} + z_{power})^2 / (p (1-p) (ln HR)^2),  p = 1/2
Total sample size N = D / P(death within T | alive at x0), averaged over the two arms.
Mortality: Gompertz-Makeham fitted in fit_mortality.py.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, chi2

from common import load_results, update_results, save_csv, cumhaz, life_expectancy, FIG

rng = np.random.default_rng(20261110)
R = load_results()["part1_mortality_fits"]
P = R["parameters_for_reuse"]
cal = {r["size"]: r for r in R["size_class_calibration"]["rows"]}
models = {
    "dog_large_26to45kg_age7": (dict(a=cal["large"]["A_a"], b=cal["large"]["A_b"], c=cal["large"]["c"]), 7.0),
    "dog_all_US_age7": (P["montoya_makeham_1to16"], 7.0),
    "human_US_age65": (P["human_both_makeham_30to94"], 65.0),
    "human_US_age50": (P["human_both_makeham_30to94"], 50.0),
}
ALPHA, POWER = 0.05, 0.80
zsum = norm.ppf(1 - ALPHA / 2) + norm.ppf(POWER)


def events_needed(hr):
    return zsum ** 2 / (0.25 * np.log(hr) ** 2)


def p_event(m, x0, T, hr):
    H = cumhaz(x0, x0 + T, m["a"], m["b"], m["c"])
    return 0.5 * (1 - np.exp(-H)) + 0.5 * (1 - np.exp(-hr * H))


def power_given_N(m, x0, T, hr, N):
    D = N * p_event(m, x0, T, hr)
    return norm.cdf(np.sqrt(D * 0.25) * abs(np.log(hr)) - norm.ppf(1 - ALPHA / 2))


rows = []
for name, (m, x0) in models.items():
    for hr in (0.9, 0.8, 0.7, 0.6):
        for T in (1, 2, 3, 4, 5, 7, 10, 15, 20):
            pe = p_event(m, x0, T, hr)
            D = events_needed(hr)
            rows.append(dict(model=name, x0=x0, HR=hr, T_years=T, P_death=pe, deaths_needed=D, N_total=D / pe))
tab = pd.DataFrame(rows)

# life-expectancy gain implied by each HR (applied from x0 to death)
gains = []
for name, (m, x0) in models.items():
    e_c = life_expectancy(m["a"], m["b"], m["c"], x0)
    for hr in (0.9, 0.8, 0.7, 0.6):
        e_t = life_expectancy(m["a"] * hr, m["b"], m["c"] * hr, x0)
        gains.append(dict(model=name, HR=hr, e_control=e_c, e_treated=e_t, gain_years=e_t - e_c,
                          gain_pct_of_remaining=100 * (e_t - e_c) / e_c,
                          gompertz_shift_lnHR_over_b=-np.log(hr) / m["b"]))
gains = pd.DataFrame(gains)
print(gains.round(3).to_string(index=False))

# Key comparison: same HR = 0.8, 3-year trial in dogs vs humans
def nt(name, hr, T):
    r = tab[(tab.model == name) & (tab.HR == hr) & (tab.T_years == T)]
    return float(r.N_total.iloc[0])
key = {f"{n}_HR{hr}_T{T}": nt(n, hr, T) for n in models for hr in (0.8, 0.7) for T in (3, 5, 10)}
print({k: round(v) for k, v in key.items()})

# duration needed for a fixed budget of participants
def duration_for_N(name, hr, N):
    m, x0 = models[name]
    for T in np.arange(0.25, 60, 0.25):
        if power_given_N(m, x0, T, hr, N) >= POWER:
            return float(T)
    return None
durations = {f"{n}_HR{hr}_N{N}": duration_for_N(n, hr, N) for n in models for hr in (0.8, 0.7) for N in (580, 2000)}
print(durations)

# TRIAD-like design: 580 dogs, >=7 y, >=20 kg, 3 y follow-up (design facts from TAMU press release, Dec 2024)
mL, x0L = models["dog_large_26to45kg_age7"]
triad = {f"HR{hr}": float(power_given_N(mL, x0L, 3, hr, 580)) for hr in (0.9, 0.8, 0.7, 0.6, 0.5)}
print("TRIAD-like analytic power", triad)


# ---------------------------------------------------------------- Monte Carlo check (log-rank)
def sample_times(m, x0, n, hr):
    E = rng.exponential(size=n) / hr
    lo, hi = np.zeros(n), np.full(n, 60.0)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = cumhaz(x0, x0 + mid, m["a"], m["b"], m["c"]) - E
        lo = np.where(f < 0, mid, lo); hi = np.where(f >= 0, mid, hi)
    return 0.5 * (lo + hi)


def logrank_p(t, ev, g):
    order = np.argsort(t); t, ev, g = t[order], ev[order], g[order]
    n = len(t); at_risk = n - np.arange(n); at_risk1 = np.cumsum(g[::-1])[::-1]
    # handle ties by unique event times
    ut, idx = np.unique(t[ev == 1], return_index=False, return_inverse=False), None
    O1 = E1 = V = 0.0
    tt = t[ev == 1]; gg = g[ev == 1]
    for u in np.unique(tt):
        k = np.searchsorted(t, u, side="left")
        nr, nr1 = at_risk[k], at_risk1[k]
        dmask = tt == u; dd = dmask.sum(); d1 = gg[dmask].sum()
        O1 += d1; E1 += dd * nr1 / nr
        if nr > 1: V += dd * (nr1 / nr) * (1 - nr1 / nr) * (nr - dd) / (nr - 1)
    return chi2.sf((O1 - E1) ** 2 / V, 1)


def mc_power(m, x0, T, hr, N, reps=1000):
    hits = 0
    g = np.r_[np.zeros(N // 2), np.ones(N - N // 2)].astype(int)
    for _ in range(reps):
        t = np.r_[sample_times(m, x0, N // 2, 1.0), sample_times(m, x0, N - N // 2, hr)]
        ev = (t <= T).astype(int); t = np.minimum(t, T)
        t = np.round(t, 6)
        hits += logrank_p(t, ev, g) < ALPHA
    return hits / reps


mc = {}
for hr in (0.8, 0.7):
    mc[f"TRIAD_like_HR{hr}"] = dict(analytic=float(power_given_N(mL, x0L, 3, hr, 580)), monte_carlo=mc_power(mL, x0L, 3, hr, 580))
mh, x0h = models["human_US_age65"]
Nh = int(round(nt("human_US_age65", 0.8, 5)))
mc["human65_HR0.8_T5_at_designN"] = dict(N=Nh, analytic=float(power_given_N(mh, x0h, 5, 0.8, Nh)), monte_carlo=mc_power(mh, x0h, 5, 0.8, Nh, reps=300))
print("Monte Carlo check:", mc)

update_results("part5_power", dict(
    assumptions=__doc__.split("Design (ASSUMPTIONS, stated once here and repeated in findings.md):")[1].split("Number of deaths")[0].strip(),
    formula="D = (z_{1-a/2}+z_pow)^2/(0.25 ln(HR)^2); N = D / P(death in T)",
    deaths_needed={f"HR{hr}": float(events_needed(hr)) for hr in (0.9, 0.8, 0.7, 0.6)},
    N_total_key=key, duration_years_for_fixed_N=durations, TRIAD_like_power_analytic=triad,
    monte_carlo_check=mc, LE_gain=gains.to_dict(orient="records"),
    TRIAD_design_source="https://vetmed.tamu.edu/news/press-releases/dog-aging-project-triad/",
    note="Large-dog hazard uses the size-class calibration (scenario A: size acts on Gompertz slope)."))
save_csv(tab, "power_sample_size.csv", ["Series for figures/power_sample_size.png; produced by power.py"])

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for name, ls in (("dog_large_26to45kg_age7", "-"), ("dog_all_US_age7", "--"), ("human_US_age65", "-."), ("human_US_age50", ":")):
    s = tab[(tab.model == name) & (tab.HR == 0.8)]
    ax[0].loglog(s.T_years, s.N_total, ls, marker="o", label=name)
ax[0].axhline(580, color="grey", lw=.8); ax[0].text(1.05, 620, "TRIAD target enrolment (580)", fontsize=7)
ax[0].set_xlabel("follow-up (years)"); ax[0].set_ylabel("total participants for 80% power")
ax[0].set_title("Hazard ratio 0.8 (20% lower death rate)"); ax[0].legend(fontsize=7)
hrs = np.linspace(0.45, 0.95, 50)
ax[1].plot(hrs, [power_given_N(mL, x0L, 3, h, 580) for h in hrs], label="TRIAD-like: 580 dogs ≥7 y, 3 y (analytic)")
ax[1].plot([0.8, 0.7], [mc["TRIAD_like_HR0.8"]["monte_carlo"], mc["TRIAD_like_HR0.7"]["monte_carlo"]], "o", label="Monte Carlo log-rank check")
ax[1].plot(hrs, [power_given_N(mh, x0h, 3, h, 580) for h in hrs], "--", label="same design, 580 humans aged 65")
ax[1].axhline(0.8, color="grey", lw=.8); ax[1].set_xlabel("true hazard ratio"); ax[1].set_ylabel("power")
ax[1].set_title("Power of a 3-year, 580-subject lifespan trial"); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f"{FIG}/power_sample_size.png", dpi=150); plt.close()
print("done")
