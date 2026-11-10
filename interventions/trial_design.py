"""Part 6: why human lifespan trials are nearly impossible, what surrogates can do, why dogs help.

Schoenfeld (1983) events formula for a two-arm 1:1 log-rank test:
   D = (z_{1-alpha/2} + z_{1-beta})^2 / ( pi (1 - pi) (ln HR)^2 ),  pi = 0.5 allocation share
   N = D / P(death during follow-up), with P averaged over both arms.
Mortality uses the fitted Gompertz parameters (common.py); composite-endpoint
event rates for TAME are ILLUSTRATIVE JUDGEMENTS [J].
Run:  python3 trial_design.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from common import MORTALITY, gm_survival, update_results, FIG, DATA

za, zb = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
K = (za + zb) ** 2


def events_needed(hr, pi=0.5):
    return K / (pi * (1 - pi) * np.log(hr) ** 2)


def p_death(species, x0, T, hr=1.0):
    p = MORTALITY[species]
    return 1 - gm_survival(T, x0, p["a"], p["b"], p["c"], hr)


def min_detectable_hr(D, pi=0.5):
    return float(np.exp(-np.sqrt(K / (pi * (1 - pi) * D))))


rows = []
for species, ages, Ts in (("human", [50, 65, 75], [5, 10, 20]), ("dog", [7, 10], [2, 3, 4])):
    for x0 in ages:
        for T in Ts:
            for hr in (0.95, 0.9, 0.8):
                pd_avg = 0.5 * (p_death(species, x0, T) + p_death(species, x0, T, hr))
                D = events_needed(hr)
                rows.append(dict(species=species, start_age=x0, follow_up_years=T, HR=hr,
                                 years_gained_first_order=-np.log(hr) / MORTALITY[species]["b"],
                                 events_needed=D, p_death_in_window=pd_avg, N_total=D / pd_avg))
td = pd.DataFrame(rows)
td.to_csv(f"{DATA}/fig_trial_sample_size.csv", index=False)
print(td.round(3).to_string())

# Figure
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, species in zip(axs, ("human", "dog")):
    d = td[(td.HR == 0.9)]
    d = d[d.species == species]
    for x0, g in d.groupby("start_age"):
        ax.plot(g.follow_up_years, g.N_total, "o-", label=f"start age {x0}")
    ax.set_yscale("log"); ax.set_xlabel("follow-up (years)"); ax.set_ylabel("total participants (80% power)")
    ax.set_title(f"{species}: all-cause mortality trial, HR 0.90 (fitted mortality)")
    ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_trial_sample_size.png", dpi=150); plt.close(fig)

out = {"events_needed": {str(h): float(events_needed(h)) for h in (0.95, 0.9, 0.8, 0.7)}}

# Human: years of follow-up for a half of 50-year-olds to die (what a true lifespan endpoint needs)
p = MORTALITY["human"]
Ts = np.linspace(0, 80, 8001)
S = gm_survival(Ts, 50, p["a"], p["b"], p["c"])
out["human_years_until_half_of_50y_olds_dead"] = float(Ts[np.argmax(S < 0.5)])
pdog = MORTALITY["dog"]
Sd = gm_survival(np.linspace(0, 20, 2001), 7, pdog["a"], pdog["b"], pdog["c"])
out["dog_years_until_half_of_7y_olds_dead"] = float(np.linspace(0, 20, 2001)[np.argmax(Sd < 0.5)])

# TAME-like composite: N=3000 aged 65-79, 6 y follow-up; composite annual event rate is a JUDGEMENT
tame = {}
for rate in (0.02, 0.035, 0.05):
    P = 1 - np.exp(-rate * 6)
    D = 3000 * P
    tame[f"annual_composite_rate_{rate}"] = {"expected_events": float(D), "min_detectable_HR_80pct_power": min_detectable_hr(D)}
out["TAME_like_design_illustrative"] = tame

# Dog trials as designed
dogt = {}
for name, n, x0, T in (("STAY-like (1300 dogs >=10 y, 4 y)", 1300, 10, 4), ("TRIAD-like (580 dogs >=7 y, 3 y)", 580, 7, 3)):
    P = p_death("dog", x0, T)
    dogt[name] = {"p_death_control": float(P), "expected_events": float(n * P),
                  "min_detectable_HR_80pct_power": min_detectable_hr(n * P)}
out["dog_trials_illustrative"] = dogt

# Surrogate endpoint: continuous biomarker with standardised effect d
sur = {}
for d in (0.1, 0.25, 0.4):
    sur[f"d={d}"] = {"n_per_arm": float(2 * K / d ** 2)}
out["biomarker_trial_n_per_arm"] = sur
# If a clock predicts mortality with HR h per SD, a d-SD shift implies HR h**(-d) ONLY IF the
# intervention acts on mortality entirely through the clock (Prentice criterion). Illustrative h [J]:
out["implied_HR_if_clock_fully_causal_illustrative"] = {
    f"h_per_SD={h}, d={d}": float(h ** (-d)) for h in (1.3, 1.6) for d in (0.25, 0.4)}
print(out)
update_results("trial_design", out)
