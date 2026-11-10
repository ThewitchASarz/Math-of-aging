"""Trial fingerprint: did a treatment delay death (level shift) or slow aging (slope change)?

Model (Gompertz, per-animal entry age x0, time on trial t):
    control hazard   h(t) = a * exp(b*x0) * exp(b*t)
    treated hazard   h(t) = r * a * exp(b*x0) * exp(s*b*t)
r is the hazard ratio at the moment treatment starts (r < 1: risk cut, "delayed death").
s is the multiplier on the rate of aging from then on (s < 1: "slowed aging").

Four nested models are fitted by maximum likelihood, with left truncation at entry and
right censoring at the end of follow-up:
    M0  r = 1, s = 1      no effect
    M1  r free, s = 1     level shift only (ordinary proportional hazards)
    M2  r = 1, s free     slope change only
    M3  r free, s free    both
Likelihood-ratio tests:
    any effect      M3 vs M0   (2 degrees of freedom)
    touched aging   M3 vs M1   (1 degree of freedom): is s different from 1?
    level shift     M3 vs M2   (1 degree of freedom): is r different from 1?

Use on real data:
    from fingerprint import fingerprint
    result = fingerprint(entry_age, years_followed, died, treated)
or from the command line:
    python3 fingerprint.py trial.csv      # columns: entry_age, years_followed, died, treated
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import exp1
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))

# Fitted dog mortality (UK VetCompass, Gompertz, ages 2-19); see dog_lifespan/results.json
_FALLBACK_DOG = dict(a=0.004967182098875431, b=0.31954873909049325)


def dog_parameters():
    path = os.path.join(HERE, "..", "dog_lifespan", "results.json")
    try:
        with open(path) as fh:
            d = json.load(fh)["part1_mortality_fits"]["RECOMMENDED_dog_parameters"]
        return dict(a=d["a_per_year"], b=d["b_per_year"])
    except (OSError, KeyError, ValueError):
        return dict(_FALLBACK_DOG)


# ----------------------------------------------------------------------------- likelihood
def _negloglik(theta, x0, t, d, arm, free_r, free_s, tau=None):
    """theta = [ln a, ln b, (ln r), (ln s)] with r, s present only when free.

    tau = years of treatment. None means treated for the whole follow-up. If treatment stops at tau,
    the risk cut r ends with it, while the aging "saved" during treatment is kept: afterwards the
    animal ages at the normal rate from a hazard that is exp((s-1)*b*tau) times the control's.
    """
    la, lb = theta[0], theta[1]
    k = 2
    lr = theta[k] if free_r else 0.0
    k += 1 if free_r else 0
    ls = theta[k] if free_s else 0.0
    b = np.exp(lb)
    log_zc = la + b * x0                            # log control hazard at entry
    lr_i = lr * arm
    bs = b * np.exp(ls * arm)                       # rate of aging while on treatment
    t1 = t if tau is None else np.minimum(t, tau)   # years on treatment
    t2 = t - t1                                     # years after treatment stopped
    with np.errstate(over="ignore", invalid="ignore"):
        cum = (np.exp(log_zc + lr_i) * np.expm1(bs * t1) / bs
               + np.exp(log_zc + bs * t1) * np.expm1(b * t2) / b)
        logh = np.where(t2 > 0, log_zc + bs * t1 + b * t2, log_zc + lr_i + bs * t1)
        val = np.sum(cum) - np.sum(d * logh)
    return float(val) if np.isfinite(val) else 1e300


def _fit(x0, t, d, arm, free_r, free_s, start, tau=None):
    best = None
    starts = [start]
    if free_r or free_s:                            # a second start guards against a poor local optimum
        alt = np.array(start, dtype=float).copy()
        alt[2:] = -0.3
        starts.append(alt)
    for s0 in starts:
        res = minimize(_negloglik, s0, args=(x0, t, d, arm, free_r, free_s, tau), method="L-BFGS-B",
                       bounds=[(-15, 2), (-5, 2)] + [(-4, 4)] * (len(s0) - 2),
                       options=dict(ftol=1e-13, gtol=1e-9, maxiter=2000, maxfun=20000))
        if best is None or res.fun < best.fun:
            best = res
    return best


def fingerprint(entry_age, years_followed, died, treated, alpha=0.05, treatment_years=None):
    """Fit M0-M3 and return estimates, likelihood-ratio tests and a plain-language verdict.

    treatment_years: how long treatment lasted, if shorter than follow-up (None = treated throughout)."""
    x0 = np.asarray(entry_age, float)
    t = np.maximum(np.asarray(years_followed, float), 1e-9)
    d = np.asarray(died, float)
    arm = np.asarray(treated, float)
    m0 = _fit(x0, t, d, arm, False, False, np.array([np.log(0.005), np.log(0.3)]), tau=treatment_years)
    base = m0.x
    m1 = _fit(x0, t, d, arm, True, False, np.r_[base, 0.0], tau=treatment_years)
    m2 = _fit(x0, t, d, arm, False, True, np.r_[base, 0.0], tau=treatment_years)
    m3 = _fit(x0, t, d, arm, True, True, np.r_[m1.x, 0.0], tau=treatment_years)
    m1b = _fit(x0, t, d, arm, True, False, m3.x[:3], tau=treatment_years)
    m2b = _fit(x0, t, d, arm, False, True, np.r_[m3.x[:2], m3.x[3]], tau=treatment_years)
    m1 = m1b if m1b.fun < m1.fun else m1
    m2 = m2b if m2b.fun < m2.fun else m2
    ll = {k: -m.fun for k, m in dict(M0=m0, M1=m1, M2=m2, M3=m3).items()}
    ll["M3"] = max(ll["M3"], ll["M1"], ll["M2"])   # nested: M3 can never fit worse
    p_any = float(chi2.sf(2 * (ll["M3"] - ll["M0"]), 2))
    p_slope = float(chi2.sf(max(2 * (ll["M3"] - ll["M1"]), 0.0), 1))
    p_level = float(chi2.sf(max(2 * (ll["M3"] - ll["M2"]), 0.0), 1))
    out = dict(
        n=int(len(x0)), deaths=int(d.sum()),
        a=float(np.exp(m3.x[0])), b=float(np.exp(m3.x[1])),
        hazard_ratio_at_start_r=float(np.exp(m3.x[2])), rate_of_aging_multiplier_s=float(np.exp(m3.x[3])),
        level_only_hazard_ratio=float(np.exp(m1.x[2])), slope_only_multiplier=float(np.exp(m2.x[2])),
        loglik=ll, p_any_effect=p_any, p_touched_aging=p_slope, p_level_shift=p_level,
    )
    if p_any >= alpha:
        verdict = "No effect detected."
    elif p_slope < alpha and out["rate_of_aging_multiplier_s"] < 1:
        verdict = "Effect detected, and the rate of aging itself is slower (slope change)."
    elif p_slope < alpha:
        verdict = "Effect detected, with a faster rate of aging on treatment (slope change, harmful direction)."
    else:
        verdict = ("Effect detected, but the data cannot tell slowed aging from delayed death "
                   "(a level shift alone explains it).")
    out["verdict"] = verdict
    return out


# ----------------------------------------------------------------------------- simulation
def remaining_life(x0, a, b, r=1.0, s=1.0):
    """Mean remaining life at age x0 under the trial model (closed form, Theorem 1)."""
    bs = b * s
    m = r * a * np.exp(b * x0) / bs
    return float(np.exp(m) * exp1(m) / bs)


def effect_for_gain(gain, x_ref, a, b, kind):
    """The r (kind='level') or s (kind='slope') that adds `gain` years of mean remaining life at x_ref."""
    base = remaining_life(x_ref, a, b)
    if kind == "level":
        return brentq(lambda r: remaining_life(x_ref, a, b, r=r) - base - gain, 1e-6, 1.0)
    return brentq(lambda s: remaining_life(x_ref, a, b, s=s) - base - gain, 1e-3, 1.0)


def draw_entry_ages(rng, n, lo, hi, a, b):
    """Entry ages of living dogs between lo and hi: density proportional to Gompertz survival."""
    grid = np.linspace(lo, hi, 2001)
    w = np.exp(-(a / b) * np.expm1(b * grid))
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, grid)


def simulate_trial(rng, n, lo, hi, follow_up, a, b, r=1.0, s=1.0, treatment_years=None):
    x0 = draw_entry_ages(rng, n, lo, hi, a, b)
    arm = (np.arange(n) % 2).astype(float)          # 1:1 allocation
    rng.shuffle(arm)
    zc = a * np.exp(b * x0)
    r_i = np.where(arm == 1, r, 1.0)
    bs = b * np.where(arm == 1, s, 1.0)
    need = -np.log(rng.random(n))                   # cumulative hazard at which each animal dies
    tau = np.inf if treatment_years is None else treatment_years
    with np.errstate(over="ignore", invalid="ignore"):
        cum_on = r_i * zc * np.expm1(bs * min(tau, 1e3)) / bs if np.isfinite(tau) else np.full(n, np.inf)
        t_on = np.log1p(bs * need / (r_i * zc)) / bs
        t_off = tau + np.log1p(b * (need - cum_on) / (zc * np.exp(bs * tau))) / b if np.isfinite(tau) else t_on
    time_to_death = np.where(need <= cum_on, t_on, t_off)
    d = (time_to_death <= follow_up).astype(float)
    t = np.minimum(time_to_death, follow_up)
    return x0, t, d, arm


def power(n, lo, hi, follow_up, a, b, r=1.0, s=1.0, sims=400, seed=0, alpha=0.05,
          treatment_years=None, analysed_as_treated_for=None):
    """analysed_as_treated_for lets the analysis assume a different treatment length from the truth
    ("whole" = analyse as if treated throughout), to show what a mis-specified analysis does."""
    rng = np.random.default_rng(seed)
    analysis_tau = treatment_years if analysed_as_treated_for is None else (
        None if analysed_as_treated_for == "whole" else analysed_as_treated_for)
    any_eff = touched = level = picked_slope = slower = 0
    deaths = []
    for _ in range(sims):
        f = fingerprint(*simulate_trial(rng, n, lo, hi, follow_up, a, b, r, s, treatment_years),
                        alpha=alpha, treatment_years=analysis_tau)
        any_eff += f["p_any_effect"] < alpha
        touched += f["p_touched_aging"] < alpha
        slower += (f["p_touched_aging"] < alpha) and (f["rate_of_aging_multiplier_s"] < 1)
        level += f["p_level_shift"] < alpha
        picked_slope += f["loglik"]["M2"] > f["loglik"]["M1"]
        deaths.append(f["deaths"])
    return dict(sims=sims, mean_deaths=float(np.mean(deaths)),
                p_detect_any_effect=any_eff / sims,
                p_conclude_rate_of_aging_changed=touched / sims,
                p_conclude_rate_of_aging_slower=slower / sims,
                p_conclude_level_shifted=level / sims,
                share_where_slope_model_fits_better=picked_slope / sims)


def triad_as_designed(a, b, sims=400):
    """The rapamycin dog trial as actually designed: 580 dogs aged 7+, ONE year of treatment,
    then two more years of follow-up. Entry ages and mortality are still simplified."""
    n, lo, hi, fu, tau = 580, 7.0, 13.0, 3.0, 1.0
    rows = {
        "no effect": power(n, lo, hi, fu, a, b, sims=sims, seed=61, treatment_years=tau),
        "risk cut 40% during the treatment year only": power(n, lo, hi, fu, a, b, r=0.6, sims=sims, seed=62, treatment_years=tau),
        "aging halved during the treatment year (kept afterwards)": power(n, lo, hi, fu, a, b, s=0.5, sims=sims, seed=63, treatment_years=tau),
        "WRONG ANALYSIS: risk cut 40% for one year, analysed as if treated throughout":
            power(n, lo, hi, fu, a, b, r=0.6, sims=sims, seed=62, treatment_years=tau, analysed_as_treated_for="whole"),
    }
    return dict(design="580 dogs, entry ages 7-13, 1 year of treatment, 3 years of follow-up", scenarios=rows,
                age_shift_from_halving_aging_for_one_year_years=0.5)


DESIGNS = {
    # name: (dogs, youngest entry age, oldest entry age, years of follow-up)
    "STAY-like (1,300 dogs aged 10+, 4 years)": (1300, 10.0, 16.0, 4.0),
    "TRIAD-like (580 dogs aged 7+, 3 years)": (580, 7.0, 13.0, 3.0),
}


def main():
    p = dog_parameters()
    a, b = p["a"], p["b"]
    out = dict(
        note=("Simulated trials under the fitted UK dog Gompertz law. Entry ages, allocation and follow-up "
              "are simplified versions of the real designs: 1:1 allocation, no dropout, treatment for the "
              "whole follow-up, all-cause death the only endpoint. Not a statement about the real trials' "
              "analysis plans."),
        dog_parameters=p, alpha=0.05)

    # effects matched on the same gain in mean remaining life at age 10
    x_ref = 10.0
    effects = {}
    for gain in (0.5, 1.0):
        r = effect_for_gain(gain, x_ref, a, b, "level")
        s = effect_for_gain(gain, x_ref, a, b, "slope")
        effects[f"+{gain} y"] = dict(level_r=r, slope_s=s)
    out["effects_matched_at_age_10"] = effects
    out["baseline_remaining_life_at_10"] = remaining_life(x_ref, a, b)

    table = {}
    for name, (n, lo, hi, fu) in DESIGNS.items():
        rows = {"no effect": power(n, lo, hi, fu, a, b, sims=600, seed=11)}
        for label, e in effects.items():
            rows[f"delayed death, {label}"] = power(n, lo, hi, fu, a, b, r=e["level_r"], sims=400, seed=21)
            rows[f"slowed aging, {label}"] = power(n, lo, hi, fu, a, b, s=e["slope_s"], sims=400, seed=31)
        table[name] = rows
        print(name)
        for k, v in rows.items():
            print(f"  {k:28s} deaths {v['mean_deaths']:6.0f}  any {v['p_detect_any_effect']:.2f}  "
                  f"aging-changed {v['p_conclude_rate_of_aging_changed']:.2f}")
    out["power_by_design"] = table

    # how big would a trial have to be to tell slowed aging from delayed death?
    sweep = []
    for n in (580, 1300, 2600, 5200):
        for lo, hi, fu in ((10.0, 16.0, 4.0), (7.0, 13.0, 4.0)):
            row = dict(dogs=n, entry_ages=f"{lo:.0f}-{hi:.0f}", follow_up_years=fu)
            for label, e in effects.items():
                res = power(n, lo, hi, fu, a, b, s=e["slope_s"], sims=300, seed=41)
                row["mean_deaths"] = res["mean_deaths"]
                row[f"power_to_show_aging_slowed_gain_{label[1:4]}y"] = res["p_conclude_rate_of_aging_changed"]
            fa = power(n, lo, hi, fu, a, b, r=effects["+1.0 y"]["level_r"], sims=300, seed=51)
            row["false_alarm_when_only_delayed"] = fa["p_conclude_rate_of_aging_changed"]
            sweep.append(row)
            print(row)
    out["size_needed_to_show_slowed_aging"] = sweep

    out["triad_as_designed"] = triad_as_designed(a, b)

    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", "power_sweep.csv"), "w") as fh:
        fh.write("# Simulated power to conclude that the rate of aging changed, when the true effect is slowed aging "
                 "worth +0.5 or +1 year of remaining life at age 10\n")
        cols = list(sweep[0].keys())
        fh.write(",".join(cols) + "\n")
        for r_ in sweep:
            fh.write(",".join(str(r_[k]) for k in cols) + "\n")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--triad-only":
        p_ = dog_parameters()
        with open(os.path.join(HERE, "results.json")) as fh:
            saved = json.load(fh)
        saved["triad_as_designed"] = triad_as_designed(p_["a"], p_["b"])
        with open(os.path.join(HERE, "results.json"), "w") as fh:
            json.dump(saved, fh, indent=1, default=float)
        print(json.dumps(saved["triad_as_designed"], indent=1))
    elif len(sys.argv) > 1:
        import pandas as pd
        df = pd.read_csv(sys.argv[1], comment="#")
        print(json.dumps(fingerprint(df["entry_age"], df["years_followed"], df["died"], df["treated"]), indent=1))
    else:
        main()
