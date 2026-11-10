"""Part 4: limits to lifespan. Extreme-value predictions of the largest age at death in a cohort of N births
under (i) Gompertz-Makeham, (ii) gamma-Gompertz (logistic deceleration as extrapolated by SSA),
(iii) Gompertz-Makeham to age 105 then a constant hazard at the level estimated from validated
semi-/supercentenarian data (Barbi et al. 2018: 0.645/yr for 1904 cohort; Belzile et al. 2021: ~0.69/yr),
(iv) the SSA 2023 table itself. Also: chance of a Calment-like record (122.45 y) given n people reaching 110.
Female parameters are used throughout (all top-15 verified record holders are women).
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from common import load_ssa, update_results, H_gm, H_ggm, RESULTS, FIG, DATA

fits = json.load(open(RESULTS))["fits"]["models"]["female"]
gm = fits["gompertz_makeham"]["params"]
gg = fits["gamma_gompertz"]["params"]
d = load_ssa()
mu_ssa = d["mu_female"].values
H_ssa_cum = np.concatenate([[0], np.cumsum(mu_ssa)])
PLATEAU = {"barbi2018_1904cohort": 0.645, "belzile2021_istat": 0.69}
CALMENT_AGE = 122 + 164 / 365.25


def H_model(name, x):
    x = np.asarray(x, float)
    if name == "gompertz_makeham":
        return H_gm(x, gm["a"], gm["b"], gm["c"])
    if name == "gamma_gompertz_from_70":
        # SSA table to age 70, then gamma-Gompertz hazard (ages 70-110 fit) beyond
        h70 = H_ssa_cum[70]
        return np.where(x <= 70, np.interp(x, np.arange(121), H_ssa_cum),
                        h70 + H_ggm(x, gg["a"], gg["b"], 0, gg["s2"]) - H_ggm(70.0, gg["a"], gg["b"], 0, gg["s2"]))
    if name.startswith("gm_then_plateau"):
        lam = PLATEAU[name.split(":")[1]]
        h105 = H_gm(105.0, gm["a"], gm["b"], gm["c"])
        return np.where(x <= 105, H_gm(x, gm["a"], gm["b"], gm["c"]), h105 + lam * (x - 105))
    if name == "ssa_table":
        return np.where(x <= 120, np.interp(x, np.arange(121), H_ssa_cum), H_ssa_cum[120] + mu_ssa[-1] * (x - 120))
    raise ValueError(name)


def quantile_max(name, N, p):
    """Age x such that P(max of N lifetimes <= x) = p, using P = exp(-N S(x)) (Poisson approx, exact as S->0)."""
    target = -np.log(p) / N  # S(x) = target
    f = lambda x: H_model(name, x) - (-np.log(target))
    return brentq(f, 1, 400)


models = ["gompertz_makeham", "gamma_gompertz_from_70", "gm_then_plateau:barbi2018_1904cohort",
          "gm_then_plateau:belzile2021_istat", "ssa_table"]
Ns = [10 ** k for k in range(3, 11)]
rows = []
for mname in models:
    for N in Ns:
        rows.append({"model": mname, "N": N, "median_max_age": quantile_max(mname, N, 0.5),
                     "lo95": quantile_max(mname, N, 0.025), "hi95": quantile_max(mname, N, 0.975)})
df = pd.DataFrame(rows)
df.to_csv(f"{DATA}/max_age_vs_cohort_size.csv", index=False)

# closed-form check (pure Gompertz, c=0 approx): median max = (1/b) ln(1 + (b/a) ln(N/ln2))
a, b, c = gm["a"], gm["b"], gm["c"]
closed = {}
for N in [1e6, 1e9]:
    pure = np.log1p((b / a) * np.log(N / np.log(2))) / b
    f = lambda x: H_gm(x, a, b, 0.0) - np.log(N / np.log(2))
    closed[f"N={N:.0e}"] = {"closed_form": float(pure), "numeric": float(brentq(f, 1, 300))}

# Calment: chance that among n people alive at 110 someone reaches 122.45, with constant hazard lam
calm = {}
for key, lam in PLATEAU.items():
    p1 = float(np.exp(-lam * (CALMENT_AGE - 110)))
    calm[key] = {"hazard_per_year": lam, "p_one_110yo_reaches_122.45": p1,
                 "n_at_110_for_50pct_chance": float(np.log(2) / p1),
                 "p_record_given_n": {str(n): float(-np.expm1(-n * p1)) for n in [100, 1000, 10000]}}
p1gm = float(np.exp(-(H_gm(CALMENT_AGE, a, b, c) - H_gm(110.0, a, b, c))))
calm["gompertz_makeham_female"] = {"p_one_110yo_reaches_122.45": p1gm, "n_at_110_for_50pct_chance": float(np.log(2) / p1gm)}

summary = {
    "note": "N = number of births (female) in a cohort experiencing the stated hazards; the 'n at 110' counts in "
            "the Calment calculation are ILLUSTRATIVE grid values, not data.",
    "median_max_age_table": {m: {str(N): round(v, 2) for N, v in zip(df[df.model == m].N, df[df.model == m].median_max_age)}
                             for m in models},
    "gain_per_10x_cohort_at_N=1e9": {m: float(df[(df.model == m) & (df.N == 10**10)].median_max_age.iloc[0]
                                          - df[(df.model == m) & (df.N == 10**9)].median_max_age.iloc[0]) for m in models},
    "gompertz_closed_form_check": closed,
    "calment_record": calm,
    "record_data": {"Jeanne Calment (France)": "122 y 164 d, died 1997", "Kane Tanaka (Japan)": "119 y 107 d, died 2022",
                    "Sarah Knauss (USA)": "119 y 97 d, died 1999",
                    "oldest living per page (Ethel Caterham, UK)": "117 y 79 d as listed when opened 2026-11-10",
                    "source": "https://en.wikipedia.org/wiki/List_of_the_verified_oldest_people"},
}
update_results("limits", summary)

fig, ax = plt.subplots(figsize=(8.5, 5))
labels = {"gompertz_makeham": "Gompertz-Makeham (no plateau)",
          "gamma_gompertz_from_70": "gamma-Gompertz (SSA-type deceleration, plateau ~2.1/yr)",
          "gm_then_plateau:barbi2018_1904cohort": "GM to 105, then plateau 0.645/yr (Barbi 2018)",
          "gm_then_plateau:belzile2021_istat": "GM to 105, then plateau 0.69/yr (Belzile 2021)",
          "ssa_table": "SSA 2023 female table as published"}
for i, m in enumerate(models):
    sub = df[df.model == m]
    ax.plot(sub.N, sub.median_max_age, "o-", color=f"C{i}", label=labels[m])
    ax.fill_between(sub.N, sub.lo95, sub.hi95, color=f"C{i}", alpha=0.12)
ax.axhline(CALMENT_AGE, color="k", ls=":", lw=1)
ax.text(1.2e3, CALMENT_AGE + 0.4, "Calment record 122.45 y", fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("cohort size N (births, female, US-2023 hazards)")
ax.set_ylabel("age of the longest-lived member (median, 95% band)")
ax.set_title("How old does the oldest person get? Extreme-value prediction by model")
ax.legend(fontsize=7.5); ax.grid(alpha=0.3, which="both")
plt.tight_layout(); plt.savefig(f"{FIG}/max_age_vs_cohort_size.png", dpi=150)

if __name__ == "__main__":
    print(json.dumps(summary, indent=1, default=float))
