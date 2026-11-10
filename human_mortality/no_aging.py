"""Part 5: the 'no aging' counterfactual. If the death rate stayed forever at its age-20 or age-30 level
(negligible senescence), what are life expectancy, median lifespan and the survival curve?
Inputs: SSA 2023 period life table (https://www.ssa.gov/oact/STATS/table4c6.html) and the Makeham
term c from results.json['fits'].
Variant A: constant hazard from birth (lifetime = exponential): e = 1/mu, median = ln2/mu.
Variant B: real 2023 hazards up to age x0, then frozen at mu(x0): e0 = int_0^x0 S + S(x0)/mu(x0).
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load_ssa, update_results, RESULTS, FIG, DATA

d = load_ssa()
fits = json.load(open(RESULTS))["fits"]["models"]
out = {"units": "years; hazards per year", "source_url": "https://www.ssa.gov/oact/STATS/table4c6.html", "data_year": 2023}
curves = {"age": np.arange(0, 3001)}


def frozen(sex, x0):
    mu = d[f"mu_{sex}"].values.copy()
    mu_f = np.concatenate([mu[:x0], np.full(3001 - x0, mu[x0])])
    S = np.concatenate([[1.0], np.exp(-np.cumsum(mu_f))])[:3001]
    l = np.concatenate([[1.0], np.exp(-np.cumsum(mu[:x0]))])
    L = l[:-1] * (-np.expm1(-mu[:x0])) / mu[:x0]
    e0 = float(L.sum() + l[-1] / mu[x0])
    median = float(np.interp(0.5, S[::-1], np.arange(3001)[::-1]))
    return e0, median, S


for sex in ["male", "female"]:
    r = {}
    S_real = np.concatenate([[1.0], np.exp(-np.cumsum(d[f"mu_{sex}"].values))])
    curves[f"real_2023_{sex}"] = np.concatenate([S_real, np.zeros(3001 - len(S_real))])
    for x0 in [20, 30]:
        mu0 = float(d.loc[x0, f"mu_{sex}"])
        e0B, medB, S = frozen(sex, x0)
        r[f"age{x0}"] = {
            "hazard_per_year": mu0,
            "A_constant_from_birth": {"life_expectancy": 1 / mu0, "median": np.log(2) / mu0,
                                      "p_survive_to_100": float(np.exp(-mu0 * 100)),
                                      "p_survive_to_500": float(np.exp(-mu0 * 500)),
                                      "p_survive_to_1000": float(np.exp(-mu0 * 1000)),
                                      "age_1pct_still_alive": float(np.log(100) / mu0)},
            "B_real_until_then_frozen": {"e0": e0B, "median": medB,
                                         "p_survive_to_100": float(S[100]), "p_survive_to_500": float(S[500]),
                                         "p_survive_to_1000": float(S[1000])}}
        curves[f"frozen_at_{x0}_{sex}"] = S
    c = fits[sex]["gompertz_makeham"]["params"]["c"]
    cl = fits[sex]["gompertz_makeham"]["ci95"]["c"]
    r["makeham_c_only"] = {"hazard_per_year": c, "life_expectancy": 1 / c, "median": np.log(2) / c,
                           "life_expectancy_ci95": [1 / cl[1], 1 / cl[0]]}
    mmin = float(d[f"mu_{sex}"].min())
    r["minimum_hazard_age_8_9"] = {"hazard_per_year": mmin, "life_expectancy": 1 / mmin, "median": np.log(2) / mmin}
    r["actual_2023"] = {"e0": float(d.loc[0, f"{sex}_ex"])}
    out[sex] = r
out["context_cause_shares_NY2023"] = {
    "heart+cancer share of deaths ages 20-24": (28 + 40) / 762,
    "heart+cancer share of deaths ages 25-44": (653 + 754) / 7844,
    "source": "https://healthweb-back.health.ny.gov/statistics/vital_statistics/2023/table34a.htm"}
update_results("no_aging", out)
pd.DataFrame(curves).to_csv(f"{DATA}/no_aging_survival.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
for ax, xmax in zip(axes, [130, 3000]):
    for sex, col in [("male", "C0"), ("female", "C3")]:
        ax.plot(curves["age"], curves[f"real_2023_{sex}"], color=col, lw=2, label=f"{sex}: actual 2023 period")
        ax.plot(curves["age"], curves[f"frozen_at_20_{sex}"], color=col, ls="--", label=f"{sex}: hazard frozen at age-20 level")
        ax.plot(curves["age"], curves[f"frozen_at_30_{sex}"], color=col, ls=":", label=f"{sex}: hazard frozen at age-30 level")
    ax.set_xlim(0, xmax); ax.set_xlabel("age (years)"); ax.grid(alpha=0.3)
axes[0].set_ylabel("fraction surviving")
axes[0].set_title("Survival with and without aging (to age 130)")
axes[1].set_title("Same curves out to 3000 years")
axes[1].legend(fontsize=7.5)
plt.tight_layout(); plt.savefig(f"{FIG}/no_aging_survival.png", dpi=150)

if __name__ == "__main__":
    print(json.dumps(out, indent=1, default=float))
