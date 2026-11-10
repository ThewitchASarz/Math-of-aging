"""Part 5: transparent Monte Carlo of years of remaining life gained, per intervention and combined.

EVERY INPUT TAGGED  [J]  IS A JUDGEMENT, NOT DATA.  Inputs tagged [D] come from the sources in
data/*.csv. Gompertz-Makeham parameters are fitted (common.py).

Model for one intervention, one simulated world:
  real    ~ Bernoulli(p_real)                       [J]  is the effect real in the target species?
  harm    ~ Bernoulli(p_harm) (only if not real)    [J]  net harm instead (known safety signals)
  theta   = -ln(HR)   (log-hazard reduction applied from the start age)
     from human/dog data:  HR ~ lognormal with 90% interval [lo, hi]           [J, anchored on D]
     from mouse data:      theta = g * K * f
        g = % median-lifespan gain in mice (sex-averaged) ~ U(g_lo, g_hi)      [D]
        K = b_mouse * median_mouse_lifespan ~ U(5, 7)  converts % to log-hazard units
                                                                                [J; mouse Gompertz values unverified]
        f = fraction of the mouse effect that survives translation ~ Beta(2, 6) (mean 0.25)  [J]
  harm theta ~ -U(0, 0.05)                                                     [J]
  years gained = LE(start, HR=exp(-theta)) - LE(start)   exact Gompertz-Makeham integral.

Combination: interventions are grouped by pathway cluster [J]. Within a cluster
  theta_cluster = max_i theta_i + rho * (sum_i theta_i - max_i theta_i),   rho ~ U(0.2, 0.6)  [J]
  (anchor: rapamycin+acarbose in male mice gave 34% vs 23% and 22% alone, i.e. rho ~ 0.5 [D, different cohorts]).
Across clusters, log-hazard reductions add (independent multiplicative HRs).
Run:  python3 monte_carlo.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from common import MORTALITY, life_expectancy, update_results, FIG, DATA

rng = np.random.default_rng(20261110)
N = 200_000

# (name, grade, cluster, p_real[J], p_harm[J], mode, params)
#  mode 'hr': params = (HR_lo, HR_hi)  90% interval if real   [J anchored on data]
#  mode 'mouse': params = (g_lo, g_hi) sex-averaged mouse median gain fraction [D]
HUMAN = [
    ("Blood-pressure control", "A", "vascular", 0.95, 0.0, "hr", (0.88, 0.99)),
    ("LDL lowering", "A", "vascular", 0.95, 0.0, "hr", (0.90, 0.99)),
    ("Exercise (inactive to active)", "B", "vascular", 0.75, 0.0, "hr", (0.80, 0.97)),
    ("Vaccination (zoster, influenza)", "B", "infection", 0.7, 0.0, "hr", (0.97, 0.995)),
    ("Caloric restriction (~12%)", "B", "nutrient", 0.5, 0.05, "hr", (0.80, 0.97)),
    ("GLP-1 receptor agonist", "B", "nutrient", 0.4, 0.05, "hr", (0.85, 0.99)),
    ("SGLT2 inhibitor", "B", "nutrient", 0.3, 0.02, "hr", (0.85, 0.99)),
    ("Rapamycin / rapalogs", "B", "nutrient", 0.35, 0.15, "mouse", (0.11, 0.25)),
    ("Acarbose", "C", "nutrient", 0.3, 0.05, "mouse", (0.04, 0.14)),
    ("Metformin", "C", "nutrient", 0.15, 0.05, "hr", (0.90, 0.99)),
    ("17-alpha-estradiol", "C", "nutrient", 0.1, 0.2, "mouse", (0.05, 0.10)),
    ("Glycine", "C", "nutrient", 0.3, 0.0, "mouse", (0.04, 0.06)),
    ("Taurine", "C", "nutrient", 0.2, 0.0, "mouse", (0.05, 0.12)),
    ("GH/IGF-1 suppression (adult)", "C", "nutrient", 0.1, 0.1, "hr", (0.85, 1.0)),
    ("Senolytics", "C", "damage", 0.15, 0.1, "hr", (0.85, 0.99)),
    ("NAD precursors", "C", "damage", 0.05, 0.0, "hr", (0.95, 0.995)),
    ("Plasma exchange", "C", "damage", 0.1, 0.05, "hr", (0.85, 0.99)),
    ("Telomerase gene therapy", "C", "damage", 0.05, 0.2, "mouse", (0.13, 0.24)),
    ("Klotho", "C", "damage", 0.1, 0.05, "hr", (0.85, 0.99)),
    ("Partial reprogramming (systemic)", "D", "damage", 0.05, 0.2, "hr", (0.50, 0.95)),
    ("Follistatin gene therapy", "D", "damage", 0.05, 0.1, "hr", (0.90, 1.0)),
]
DOG = [
    ("Diet restriction / lean feeding", "A", "nutrient", 0.85, 0.0, "hr", (0.55, 0.85)),
    ("Dental cleaning", "B", "infection", 0.3, 0.02, "hr", (0.85, 0.98)),
    ("Rapamycin (TRIAD pending)", "B", "nutrient", 0.35, 0.1, "mouse", (0.11, 0.25)),
    ("LOY-002 (senior dogs)", "C", "nutrient", 0.35, 0.05, "hr", (0.70, 0.95)),
    ("LOY-001/003 IGF-1 lowering (large dogs)", "C", "nutrient", 0.3, 0.1, "hr", (0.70, 0.95)),
    ("Senolytics", "C", "damage", 0.15, 0.1, "hr", (0.85, 0.99)),
]
START = {"human": 50.0, "dog": 7.0}
START_OVERRIDE = {"LOY-002 (senior dogs)": 10.0}


def gain_interpolator(species, x0):
    p = MORTALITY[species]
    base = life_expectancy(x0, p["a"], p["b"], p["c"])
    th = np.concatenate([np.linspace(-0.2, 3, 161), np.linspace(3.05, 8, 60)])
    g = np.array([life_expectancy(x0, p["a"], p["b"], p["c"], np.exp(-t)) - base for t in th])
    return lambda t: np.interp(t, th, g), base


def lognormal_hr(lo, hi, n):
    # 90% interval [lo, hi] on HR -> theta = -ln HR normal
    t_lo, t_hi = -np.log(hi), -np.log(lo)
    mu = (t_lo + t_hi) / 2; sd = (t_hi - t_lo) / (2 * 1.645)
    return rng.normal(mu, sd, n)


def simulate(entry, n):
    name, grade, cluster, p_real, p_harm, mode, par = entry
    real = rng.random(n) < p_real
    harm = (~real) & (rng.random(n) < p_harm)
    if mode == "hr":
        th = lognormal_hr(par[0], par[1], n)
    else:
        g = rng.uniform(par[0], par[1], n)
        K = rng.uniform(5, 7, n)
        f = rng.beta(2, 6, n)
        th = g * K * f
    theta = np.where(real, th, 0.0)
    theta = np.where(harm, -rng.uniform(0, 0.05, n), theta)
    return theta


def summarise(years):
    return {"mean": float(years.mean()), "p5": float(np.percentile(years, 5)),
            "median": float(np.median(years)), "p95": float(np.percentile(years, 95)),
            "P_gain_gt_0.1y": float((years > 0.1).mean())}


results = {}
rows = []
draws = {}
for species, table in (("human", HUMAN), ("dog", DOG)):
    interp_cache = {}
    thetas = {}
    for e in table:
        x0 = START_OVERRIDE.get(e[0], START[species])
        if x0 not in interp_cache:
            interp_cache[x0] = gain_interpolator(species, x0)
        gfun, base = interp_cache[x0]
        th = simulate(e, N)
        thetas[e[0]] = th
        yrs = gfun(th)
        s = summarise(yrs)
        s.update(dict(species=species, intervention=e[0], grade=e[1], cluster=e[2], p_real=e[3],
                      p_harm=e[4], start_age=x0, baseline_remaining_LE=base))
        rows.append(s)
        draws[(species, e[0])] = yrs
    # ---- combinations (all at species default start age)
    gfun, base = interp_cache[START[species]]
    names = [e[0] for e in table]
    clusters = sorted(set(e[2] for e in table))
    T = np.vstack([thetas[nm] for nm in names])  # (k, N)
    indep = gfun(T.sum(axis=0))
    rho = rng.uniform(0.2, 0.6, N)
    tot = np.zeros(N)
    for cl in clusters:
        idx = [i for i, e in enumerate(table) if e[2] == cl]
        sub = T[idx]
        mx = sub.max(axis=0)
        tot += mx + rho * (sub.sum(axis=0) - mx)
    overl = gfun(tot)
    # established-only (grade A)
    A = [i for i, e in enumerate(table) if e[1] == "A"]
    estab = gfun(T[A].sum(axis=0)) if A else np.zeros(N)
    results[species] = {
        "start_age": START[species], "baseline_remaining_LE": base,
        "combo_all_independent": summarise(indep),
        "combo_all_overlapping_pathways": summarise(overl),
        "combo_grade_A_only": summarise(estab),
        "combined_HR_overlapping_median": float(np.exp(-np.median(tot))),
        "combined_HR_overlapping_p95_best": float(np.exp(-np.percentile(tot, 95))),
    }
    draws[(species, "ALL (overlapping)")] = overl
    draws[(species, "ALL (independent)")] = indep

tab = pd.DataFrame(rows)
order = {"A": 0, "B": 1, "C": 2, "D": 3}
tab = tab.sort_values(["species", "mean"], ascending=[True, False])
cols = ["species", "intervention", "grade", "cluster", "p_real", "p_harm", "start_age",
        "mean", "p5", "median", "p95", "P_gain_gt_0.1y"]
tab[cols].to_csv(f"{DATA}/monte_carlo_ranked.csv", index=False)
print(tab[cols].round(3).to_string())
for k, v in results.items():
    print(k, {kk: (vv if not isinstance(vv, dict) else {a: round(b, 2) for a, b in vv.items()}) for kk, vv in v.items()})

# ---- what would 'solving ageing' need, vs what we have
from common import MORTALITY as P
solve = {}
for species in ("human", "dog"):
    p = P[species]; x0 = START[species]
    mu0 = p["a"] * np.exp(p["b"] * x0)
    solve[species] = {
        "remaining_LE_now": results[species]["baseline_remaining_LE"],
        "remaining_LE_if_ageing_stopped_at_start_age": 1.0 / (mu0 + p["c"]),
        "everything_known_overlapping_mean_gain": results[species]["combo_all_overlapping_pathways"]["mean"],
        "everything_known_overlapping_p95_gain": results[species]["combo_all_overlapping_pathways"]["p95"],
    }
    s = solve[species]
    s["fraction_of_gap_closed_mean"] = s["everything_known_overlapping_mean_gain"] / (
        s["remaining_LE_if_ageing_stopped_at_start_age"] - s["remaining_LE_now"])
results["solving_ageing_comparison"] = solve
print(solve)

# ---- figure: ranked intervals
for species in ("human", "dog"):
    t = tab[tab.species == species].sort_values("mean")
    labels = list(t.intervention) + ["ALL, overlapping pathways", "ALL, independent (optimistic)"]
    means = list(t["mean"]) + [results[species]["combo_all_overlapping_pathways"]["mean"],
                               results[species]["combo_all_independent"]["mean"]]
    lo = list(t.p5) + [results[species]["combo_all_overlapping_pathways"]["p5"], results[species]["combo_all_independent"]["p5"]]
    hi = list(t.p95) + [results[species]["combo_all_overlapping_pathways"]["p95"], results[species]["combo_all_independent"]["p95"]]
    grades = list(t.grade) + ["", ""]
    fig, ax = plt.subplots(figsize=(8, 0.38 * len(labels) + 1.5))
    yy = np.arange(len(labels))
    colors = {"A": "#1b7837", "B": "#5aae61", "C": "#d9a400", "D": "#c0392b", "": "#333"}
    for i in yy:
        ax.plot([lo[i], hi[i]], [i, i], color=colors[grades[i]], lw=2)
        ax.plot(means[i], i, "o", color=colors[grades[i]])
    ax.set_yticks(yy); ax.set_yticklabels([f"{l} [{g}]" if g else l for l, g in zip(labels, grades)], fontsize=7)
    ax.set_xlabel(f"years of remaining life gained from age {START[species]:.0f} (dot = mean, bar = 5-95%)")
    ax.set_title(f"{species}: Monte Carlo (judgement inputs [J]; fitted Gompertz)", fontsize=9)
    ax.axvline(0, color="k", lw=.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_mc_ranked_{species}.png", dpi=150); plt.close(fig)
    pd.DataFrame({"label": labels, "grade": grades, "mean": means, "p5": lo, "p95": hi}).to_csv(
        f"{DATA}/fig_mc_ranked_{species}.csv", index=False)

# ---- figure: total vs 'solving ageing'
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
rowsf = []
for ax, species in zip(axs, ("human", "dog")):
    s = solve[species]
    vals = [s["remaining_LE_now"], s["remaining_LE_now"] + s["everything_known_overlapping_mean_gain"],
            s["remaining_LE_now"] + s["everything_known_overlapping_p95_gain"],
            s["remaining_LE_if_ageing_stopped_at_start_age"]]
    labs = ["today", "+ everything known\n(mean)", "+ everything known\n(95th pct)", "ageing stopped\nat start age"]
    ax.bar(labs, vals, color=["#888", "#5aae61", "#1b7837", "#2a6f97"])
    ax.set_ylabel(f"remaining life expectancy at {START[species]:.0f} (years)")
    ax.set_title(f"{species} (fitted mortality; illustrative comparison)")
    ax.tick_params(axis="x", labelsize=7)
    for l, v in zip(labs, vals):
        rowsf.append({"species": species, "scenario": l.replace("\n", " "), "remaining_LE": v})
fig.tight_layout(); fig.savefig(f"{FIG}/fig_total_vs_solving.png", dpi=150); plt.close(fig)
pd.DataFrame(rowsf).to_csv(f"{DATA}/fig_total_vs_solving.csv", index=False)

results["ranked_table"] = tab[cols].to_dict(orient="records")
results["N_draws"] = N
results["seed"] = 20261110
update_results("monte_carlo", results)

# ---- sensitivity of the combined total to the judgement inputs  [J]
def combo_total(table, species, pmult=1.0, fbeta=(2, 6), rho_rng=(0.2, 0.6), n=50_000):
    global rng
    rng = np.random.default_rng(7)
    gfun, _ = gain_interpolator(species, START[species])
    T = []
    for e in table:
        e2 = list(e); e2[3] = min(1.0, e[3] * pmult)
        name, grade, cluster, p_real, p_harm, mode, par = e2
        real = rng.random(n) < p_real
        harm = (~real) & (rng.random(n) < p_harm)
        if mode == "hr":
            th = lognormal_hr(par[0], par[1], n)
        else:
            th = rng.uniform(par[0], par[1], n) * rng.uniform(5, 7, n) * rng.beta(*fbeta, n)
        th = np.where(real, th, 0.0); th = np.where(harm, -rng.uniform(0, 0.05, n), th)
        T.append(th)
    T = np.vstack(T)
    rho = rng.uniform(*rho_rng, n); tot = np.zeros(n)
    for cl in sorted(set(e[2] for e in table)):
        idx = [i for i, e in enumerate(table) if e[2] == cl]
        mx = T[idx].max(axis=0); tot += mx + rho * (T[idx].sum(axis=0) - mx)
    y = gfun(tot)
    return {"mean": float(y.mean()), "p95": float(np.percentile(y, 95))}

sens = {}
for species, table in (("human", HUMAN), ("dog", DOG)):
    sens[species] = {
        "base": combo_total(table, species),
        "p_real_x0.5": combo_total(table, species, pmult=0.5),
        "p_real_x1.5": combo_total(table, species, pmult=1.5),
        "all_p_real_1_(everything_works)": combo_total(table, species, pmult=100),
        "translation_f_mean_0.5": combo_total(table, species, fbeta=(4, 4)),
        "rho_0_(pure_redundancy)": combo_total(table, species, rho_rng=(0.0, 1e-9)),
        "rho_1_(fully_additive)": combo_total(table, species, rho_rng=(1.0, 1.0 + 1e-9)),
    }
print("sensitivity:", sens)
pd.DataFrame([{"species": s, "scenario": k, **v} for s, d in sens.items() for k, v in d.items()]).to_csv(
    f"{DATA}/monte_carlo_sensitivity.csv", index=False)
update_results("monte_carlo_sensitivity", sens)
