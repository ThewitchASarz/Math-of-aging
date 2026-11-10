"""Part 2: what moves life expectancy.
(a) lowering baseline (age-independent) mortality, (b) slowing the rate of aging (Gompertz slope),
(c) delaying the onset of aging by T years; proportional reductions and Keyfitz entropy;
Taeuber paradox via a competing-risks cause-elimination calculation.

All empirical inputs: SSA 2023 period life table (https://www.ssa.gov/oact/STATS/table4c6.html),
GM parameters from results.json['fits'] (run fit_mortality_laws.py first),
cause-of-death shares by age: New York State 2023 vital statistics Table 34a
(https://healthweb-back.health.ny.gov/statistics/vital_statistics/2023/table34a.htm) and US totals from
NCHS Data Brief 521 (https://cdc.gov/nchs/products/databriefs/db521.htm).
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load_ssa, update_results, keyfitz_entropy, life_expectancy_from_mu, mu_gm, RESULTS, FIG, DATA

d = load_ssa()
fits = json.load(open(RESULTS))["fits"]["models"]
AMAX = 250
ages = np.arange(AMAX + 1)  # table ends at 119; beyond that the senescent part is extended log-linearly
X0 = 30  # age at which "aging" interventions start to act on the empirical table


def e0_from_mu(mu):
    """e0 from single-year hazards (piecewise-constant within year, exact integral) with open tail 1/mu."""
    mu = np.asarray(mu, float)
    l = np.concatenate([[1.0], np.exp(-np.cumsum(mu))])
    # integral within year of exp(-mu t) = (1-e^-mu)/mu
    L = l[:-1] * (-np.expm1(-mu)) / mu
    return float(L.sum() + l[-1] / mu[-1])


def mu_ext(sex):
    """Observed hazard 0-119, extended to AMAX by growing the non-Makeham part at the GM slope b
    (only matters for scenarios that push many survivors past 119)."""
    mu = d[f"mu_{sex}"].values.copy()
    c = fits[sex]["gompertz_makeham"]["params"]["c"]
    b = fits[sex]["gompertz_makeham"]["params"]["b"]
    extra = np.arange(1, AMAX - 119 + 1)
    return np.concatenate([mu, c + (mu[-1] - c) * np.exp(b * extra)])


def split(sex):
    """Observed hazard = c (age-independent, GM Makeham term) + s_x (everything else, 'senescent+other')."""
    mu = mu_ext(sex)
    c = fits[sex]["gompertz_makeham"]["params"]["c"]
    s = mu - c
    s[s <= 0] = 1e-9  # childhood: total hazard is below c; keep tiny positive
    return mu, c, s


def baseline_change(sex, frac):
    mu, c, s = split(sex)
    return np.where(ages >= 0, np.maximum(mu - frac * c, 1e-9), mu)


def slow_aging(sex, eps):
    """For x >= X0, log of senescent part grows at (1-eps) times its observed rate."""
    mu, c, s = split(sex)
    out = mu.copy()
    sel = ages >= X0
    s_new = s[X0] * (s[sel] / s[X0]) ** (1 - eps)
    out[sel] = c + s_new
    return out


def delay_aging(sex, T):
    """Aging pauses for T years at age X0: senescent part stays at s(X0) for T years, then resumes."""
    mu, c, s = split(sex)
    out = mu.copy()
    T = int(T)
    for i, x in enumerate(ages):
        if x >= X0:
            src = max(X0, x - T)
            out[i] = c + s[src]
    return out


def proportional(sex, delta):
    mu = mu_ext(sex)
    return mu * (1 - delta)


res = {"inputs": "SSA 2023 period table; GM Makeham c from fits; interventions act from age 30 unless noted",
       "units": "years of period life expectancy at birth (e0)"}
rows = []
for sex in ["male", "female"]:
    mu, c, s = split(sex)
    e0 = e0_from_mu(mu)
    H = keyfitz_entropy(d[f"{sex}_qx"].values)
    r = {"e0_baseline": e0, "keyfitz_entropy_H": H}
    # (a) baseline
    r["remove_all_makeham_c"] = e0_from_mu(baseline_change(sex, 1.0)) - e0
    r["halve_makeham_c"] = e0_from_mu(baseline_change(sex, 0.5)) - e0
    # proportional all-age reductions vs Keyfitz first-order prediction
    r["proportional"] = {}
    for dl in [0.01, 0.1, 0.2, 0.5]:
        exact = e0_from_mu(proportional(sex, dl)) - e0
        r["proportional"][str(dl)] = {"exact_gain": exact, "keyfitz_first_order_delta*H*e0": dl * H * e0}
    # (b) slow aging
    r["slow_aging"] = {str(eps): e0_from_mu(slow_aging(sex, eps)) - e0 for eps in [0.1, 0.2, 0.5, 1.0]}
    # (c) delay
    r["delay_aging"] = {str(T): e0_from_mu(delay_aging(sex, T)) - e0 for T in [1, 5, 10, 20]}
    res[sex] = r
    for eps in np.linspace(0, 0.6, 31):
        rows.append({"sex": sex, "kind": "slow_aging_eps", "x": eps, "gain": e0_from_mu(slow_aging(sex, eps)) - e0})
    for T in range(0, 31):
        rows.append({"sex": sex, "kind": "delay_T_years", "x": T, "gain": e0_from_mu(delay_aging(sex, T)) - e0})
    for dl in np.linspace(0, 0.9, 31):
        rows.append({"sex": sex, "kind": "proportional_delta", "x": dl, "gain": e0_from_mu(proportional(sex, dl)) - e0})

# ---- analytic checks on pure Gompertz (female fit) ----
a, b = fits["female"]["gompertz"]["params"]["a"], fits["female"]["gompertz"]["params"]["b"]
eG = lambda a_, b_: life_expectancy_from_mu(lambda t: mu_gm(t, a_, b_))[0]
e_base = eG(a, b)
an = {"model": "pure Gompertz, female fit (a, b) applied from birth",
      "e0": e_base,
      "delay_10y_numeric": eG(a * np.exp(-b * 10), b) - e_base,
      "delay_10y_theory": 10.0,
      "halve_hazard_numeric": eG(a / 2, b) - e_base,
      "halve_hazard_theory_ln2/b": np.log(2) / b,
      "entropy_numeric": None,
      "entropy_theory_1/(b e0)": 1 / (b * e_base)}
t = np.linspace(0, 300, 300001)
S = np.exp(-(a / b) * np.expm1(b * t))
Hn = -np.trapezoid(S * np.log(np.maximum(S, 1e-300)), t) / np.trapezoid(S, t)
an["entropy_numeric"] = float(Hn)
# slope: anchor hazard at age 30, reduce b by 10%
mu30 = a * np.exp(b * 30)
def e_anchor(bb):
    return life_expectancy_from_mu(lambda tt: np.where(tt < 30, a * np.exp(b * tt), mu30 * np.exp(bb * (tt - 30))))[0]
an["slope_minus10pct_anchor30_numeric"] = e_anchor(0.9 * b) - e_base
an["slope_minus10pct_first_order_0.1*((e0-30)-1/b)"] = 0.1 * ((e_base - 30) - 1 / b)
res["gompertz_analytic_checks"] = an

# ---- Taeuber paradox: competing risks with cause shares by age ----
ny = pd.DataFrame({  # NY State 2023, Table 34a (Total row per age group). Heart not in top 5 under 10 -> 0
    "lo": [0, 1, 10, 20, 25, 45, 65, 75, 85], "hi": [0, 9, 19, 24, 44, 64, 74, 84, 250],
    "all": [825, 303, 560, 762, 7844, 27197, 30749, 40700, 53342],
    "heart": [0, 0, 20, 28, 653, 5798, 7420, 10445, 16805],
    "cancer": [0, 42, 56, 40, 754, 7071, 9009, 9307, 5979]})
with open(f"{DATA}/ny2023_cause_shares_by_age.csv", "w") as fh_:
    fh_.write("# Source: New York State Dept of Health, Vital Statistics 2023, Table 34a (Total rows), transcribed\n"
              "# URL: https://healthweb-back.health.ny.gov/statistics/vital_statistics/2023/table34a.htm\n"
              "# heart = Diseases of the Heart deaths, cancer = Malignant Neoplasms deaths; 0 where not in the age group's top 5\n")
    ny.to_csv(fh_, index=False)
ny_heart_share = 41187 / 162287
ny_cancer_share = 32260 / 162287
us_heart_share = 680981 / 3090964
us_cancer_share = 613352 / 3090964
fr = {"heart": np.zeros(len(ages)), "cancer": np.zeros(len(ages))}
for _, rw in ny.iterrows():
    sel = (ages >= rw.lo) & (ages <= rw.hi)
    fr["heart"][sel] = rw.heart / rw["all"]
    fr["cancer"][sel] = rw.cancer / rw["all"]
taeuber = {"method": "independent competing risks: mu_x(-k) = mu_x (1 - f_k(age group)); "
                     "f from NY State 2023 deaths by age (both sexes) applied to US SSA 2023 hazards by sex; "
                     "'US-scaled' multiplies f by US/NY all-age share ratio",
           "us_2023_share_heart": us_heart_share, "us_2023_share_cancer": us_cancer_share,
           "ny_2023_share_heart": ny_heart_share, "ny_2023_share_cancer": ny_cancer_share}
for sex in ["male", "female"]:
    mu = mu_ext(sex)
    e0 = e0_from_mu(mu)
    out = {}
    for scale_name, sh, sc in [("ny_fractions", 1, 1),
                               ("us_scaled", us_heart_share / ny_heart_share, us_cancer_share / ny_cancer_share)]:
        fh, fc = np.clip(fr["heart"] * sh, 0, 1), np.clip(fr["cancer"] * sc, 0, 1)
        out[scale_name] = {
            "eliminate_heart": e0_from_mu(mu * (1 - fh)) - e0,
            "eliminate_cancer": e0_from_mu(mu * (1 - fc)) - e0,
            "eliminate_both": e0_from_mu(mu * (1 - fh - fc)) - e0,
        }
    # what share of deaths does that represent? (life-table deaths weighted)
    l = np.concatenate([[1.0], np.exp(-np.cumsum(mu))])
    dx = l[:-1] - l[1:]
    out["lifetable_prob_die_of_cancer_ny_fractions"] = float(np.sum(dx * fr["cancer"]))
    out["lifetable_prob_die_of_heart_ny_fractions"] = float(np.sum(dx * fr["heart"]))
    taeuber[sex] = out
taeuber["published_comparisons"] = {
    "NCHS 1999-2001 cause-elimination tables (NVSR 61(9), 2013)": "heart disease: almost 4 years; cancer: more than 3 years (abstract)",
    "Keyfitz 1977, US males 1964": "eliminating cancer: about 2.27 years",
    "1989 US tables (Wikipedia, Taeuber paradox page)": "cancer 3.36 years; major cardiovascular disease 6.73 years"}
res["taeuber"] = taeuber
update_results("le_sensitivity", res)
df = pd.DataFrame(rows)
df.to_csv(f"{DATA}/le_sensitivity.csv", index=False)

# figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, kind, xl in zip(axes, ["proportional_delta", "slow_aging_eps", "delay_T_years"],
                        ["fraction cut from every age's death rate (delta)",
                         "fractional slowing of the aging slope after age 30 (epsilon)",
                         "years aging is paused at age 30 (T)"]):
    for sex, col in [("male", "C0"), ("female", "C3")]:
        sub = df[(df.sex == sex) & (df.kind == kind)]
        ax.plot(sub.x, sub.gain, color=col, label=sex)
    ax.set_xlabel(xl); ax.grid(alpha=0.3)
axes[0].set_ylabel("gain in period life expectancy at birth (years)")
axes[0].set_title("(a) cut all death rates proportionally")
axes[1].set_title("(b) slow the rate of aging")
axes[2].set_title("(c) delay the onset of aging")
axes[2].plot([0, 30], [0, 30], "k:", lw=1, label="gain = T")
axes[2].legend()
plt.suptitle("What moves life expectancy: US 2023 period table (SSA), interventions applied to the hazard")
plt.tight_layout(); plt.savefig(f"{FIG}/le_sensitivity.png", dpi=150)

if __name__ == "__main__":
    print(json.dumps(res, indent=1, default=float))
