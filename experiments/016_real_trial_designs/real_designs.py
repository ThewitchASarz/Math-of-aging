"""Experiment 016: the two dog trials now running, simulated with their own published designs.

Part A  The rapamycin trial (TRIAD) as published: 580 dogs, Gompertz level 0.0066 and slope 0.29 per year,
        entry ages 7-8 (half), 8-9 (quarter), 9-10 (quarter), 12 months of dosing, 36 months of follow-up.
        Source: Coleman et al. 2025, GeroScience 47:2851 (PMC12181551), as recorded in
        interventions/triad_published_inputs.py and papers/literature-check.md. The paper itself could not
        be opened in the run that wrote this script.
Part B  A trial like Loyal's STAY study: 1,300 dogs aged 10 or older, 4 years, treatment throughout
        (published: DVM360, 13 January 2026 and Loyal, 9 October 2024, as recorded in papers/literature-check.md).
        ASSUMPTIONS: UK death rates (dog_lifespan/results.json), no dropout, and the entry ages (two guesses).
Part C  Experiment 010 repeated with Part A's dogs: treated throughout, followed for 10 years.

Model for a treated dog (t = years since enrolment, tau = years of dosing, h0 = a * exp(b * entry age)):
    t <  tau :  hazard = h0 * r1 * exp(b * s1 * t)
    t >= tau :  hazard = h0 * exp(b * s1 * tau) * r2 * exp(b * (t - tau))
Placebo dogs have r1 = s1 = r2 = 1. This is the model of experiment 003.

Analyses of each simulated trial:
    any-difference test : two-sided log-rank test at the 5% level (ASSUMPTION: the paper names a Cox model)
    slope test          : experiment 001 (likelihood-ratio test of a slope change, level effect allowed)
    two-period test     : experiment 006 (is the risk ratio after dosing different from 1?), Part A only

Run:  python3 real_designs.py            (about 2 minutes on 2 cores)
      python3 real_designs.py --quick
"""
import importlib.util
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from multiprocessing import Pool

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "experiments", "001_trial_fingerprint", "fingerprint.py")
two = _load("two_period", "experiments", "006_stop_and_watch", "two_period.py")
triad = _load("triad_published_inputs", "interventions", "triad_published_inputs.py")

SEED = 20261117
TRIAD_A, TRIAD_B = triad.A, triad.B              # 0.0066 and 0.29 per year, published
Z = 1.959964
LOOKS = (200, 400, 580)                          # early looks after 200 and 400 dogs: published. Rule: ASSUMPTION
LONG = 1e9                                       # "tau" for a treatment that never stops


def ages_triad(n, rng):
    """Published mix: half enter at 7-8, a quarter at 8-9, a quarter at 9-10."""
    return triad.entry_ages(n, rng)


def ages_even_10_14(n, rng):
    """ASSUMPTION (as in experiment 001): entry ages evenly spread from 10 to 14."""
    return rng.uniform(10.0, 14.0, n)


def ages_population_10_16(n, rng, a, b):
    """ASSUMPTION: entry ages shaped like the dog population aged 10 to 16, that is, in step with the share
    of dogs still alive at each age under the UK fit."""
    grid = np.linspace(10.0, 16.0, 6001)
    surv = np.exp(-a / b * (np.exp(b * grid) - np.exp(b * 10.0)))
    cdf = np.cumsum(surv)
    cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    return np.interp(rng.random(n), cdf, grid)


def simulate(x0, follow, tau, r1, s1, r2, a, b, rng):
    """Two-phase death times by inverse transform (the formulas of experiment 003, any entry ages, any tau)."""
    n = x0.size
    arm = (np.arange(n) % 2).astype(float)
    treated = arm == 1
    R1, S1, R2 = (np.where(treated, v, 1.0) for v in (r1, s1, r2))
    h0 = a * np.exp(b * x0)
    e = -np.log(rng.random(n))
    bs1 = b * S1
    if tau >= LONG:
        t = np.log1p(e * bs1 / (h0 * R1)) / bs1
    else:
        H1 = h0 * R1 * np.expm1(bs1 * tau) / bs1
        t_first = np.log1p(np.minimum(e, H1) * bs1 / (h0 * R1)) / bs1
        h2 = h0 * np.exp(bs1 * tau) * R2
        t_second = tau + np.log1p(np.maximum(e - H1, 0.0) * b / h2) / b
        t = np.where(e < H1, t_first, t_second)
    died = (t <= follow).astype(float)
    return x0, arm, np.minimum(t, follow), died


def obf_constant():
    """Constant c of an O'Brien-Fleming-type rule: stop at look k if |z_k| > c / sqrt(share of dogs seen).
    Chosen so that a trial with no effect stops or finishes 'significant' 5% of the time over the three looks."""
    rng = np.random.default_rng(SEED)
    frac = np.array(LOOKS) / LOOKS[-1]
    steps = rng.standard_normal((400000, 3)) * np.sqrt(np.diff(np.concatenate([[0.0], frac])))
    worst = np.max(np.abs(np.cumsum(steps, axis=1)), axis=1)     # |z_k| * sqrt(frac_k) = |B(frac_k)|
    return float(np.quantile(worst, 0.95))


# Cases: (part, label, n dogs, entry ages, years of follow-up, years of dosing, (r1, s1, r2), analyses)
def build_cases():
    A = [("no effect", (1.0, 1.0, 1.0)),
         ("risk cut 20% in the dosing year, then fades", (0.8, 1.0, 1.0)),
         ("risk cut 50% in the dosing year, then fades", (0.5, 1.0, 1.0)),
         ("risk cut 20%, lasts all 3 years", (0.8, 1.0, 0.8)),
         ("risk cut 38%, lasts all 3 years", (0.62, 1.0, 0.62)),
         ("aging slowed 20% in the dosing year, gain kept", (1.0, 0.8, 1.0)),
         ("aging stopped in the dosing year, gain kept", (1.0, 1e-6, 1.0))]
    B = [("no effect", (1.0, 1.0, 1.0)), ("risk cut 20%", (0.8, 1.0, 1.0)), ("risk cut 30%", (0.7, 1.0, 1.0)),
         ("aging slowed 20%", (1.0, 0.8, 1.0)), ("aging slowed 30%", (1.0, 0.7, 1.0))]
    C = [("no effect", (1.0, 1.0, 1.0)), ("risk cut 50%", (0.5, 1.0, 1.0)), ("aging slowed 20%", (1.0, 0.8, 1.0))]
    cases = [("A", "rapamycin trial as published", 580, "triad", 3.0, 1.0, name, eff) for name, eff in A]
    for ages in ("even 10-14", "population 10-16"):
        cases += [("B", f"STAY-like, entry ages {ages}", 1300, ages, 4.0, LONG, name, eff) for name, eff in B]
    for n in (580, 1300):
        cases += [("C", f"rapamycin trial's dogs, treated throughout, 10 years, {n} dogs", n, "triad", 10.0, LONG,
                   name, eff) for name, eff in C]
    return cases


def one_case(args):
    (part, design, n, ages, follow, tau, effect, (r1, s1, r2)), nsim, seed, c_obf = args
    uk_a, uk_b = fp.dog_parameters()
    a, b = (uk_a, uk_b) if part == "B" else (TRIAD_A, TRIAD_B)
    rng = np.random.default_rng(seed)
    c1 = chi2.ppf(0.95, 1)
    logrank = slope = slower = lasting = seq = seq_early = 0
    s_hat, placebo_dead = [], []
    for _ in range(nsim):
        if ages == "triad":
            x0 = ages_triad(n, rng)
        elif ages == "even 10-14":
            x0 = ages_even_10_14(n, rng)
        else:
            x0 = ages_population_10_16(n, rng, a, b)
        data = simulate(x0, follow, tau, r1, s1, r2, a, b, rng)
        _, arm, t, died = data
        placebo_dead.append(float(died[arm == 0].mean()))
        logrank += abs(triad.logrank_z(t, died, arm)) > Z
        ll = fp.fit_all(data, a, b)
        flag = 2 * (ll["full"][0] - ll["level"][0]) > c1
        sh = float(np.exp(ll["full"][1][3]))
        slope += flag
        slower += flag and sh < 1.0
        s_hat.append(sh)
        if part == "A":
            on = two.fit("on_only", data, a, b)
            both = two.fit("both", data, a, b, extra=[on[1][2], 0.0])
            lasting += 2 * (both[0] - on[0]) > c1
            # Early looks. ASSUMPTION: dogs are counted in order of enrolment and each look uses complete
            # 3-year records of the dogs seen so far.
            for k, m in enumerate(LOOKS):
                z = triad.logrank_z(t[:m], died[:m], arm[:m])
                if abs(z) > c_obf / np.sqrt(m / LOOKS[-1]):
                    seq += 1
                    seq_early += k < len(LOOKS) - 1
                    break
    row = dict(part=part, design=design, effect=effect, n_dogs=n, follow_up_years=follow,
               dosing_years=None if tau >= LONG else tau, n_trials=nsim,
               placebo_share_dead=float(np.mean(placebo_dead)),
               any_difference_detected=logrank / nsim, slope_change_flagged=slope / nsim,
               flagged_as_slower_aging=slower / nsim, slope_estimate_median=float(np.median(s_hat)))
    if part == "A":
        row.update(lasting_benefit_reported=lasting / nsim, any_difference_with_early_looks=seq / nsim,
                   stopped_at_an_early_look=seq_early / nsim)
    return row


def risk_cut_for_one_year(a, b, ages):
    """The lasting risk ratio r that adds one year of life on average over the given entry ages."""
    gain = lambda r: np.mean([fp.years_gained(x, r, 1.0, a, b) for x in ages]) - 1.0
    return float(brentq(gain, 0.05, 0.99))


def main():
    quick = "--quick" in sys.argv
    nsim = 30 if quick else 1000
    c_obf = obf_constant()
    cases = [(c, nsim, SEED + 100 * i, c_obf) for i, c in enumerate(build_cases())]
    with Pool(2) as pool:
        rows = pool.map(one_case, cases, chunksize=1)
    uk_a, uk_b = fp.dog_parameters()
    rng = np.random.default_rng(SEED)
    pop_ages = ages_population_10_16(4000, rng, uk_a, uk_b)
    side = {
        "obf_constant": c_obf,
        "obf_thresholds_for_z_at_each_look": [c_obf / float(np.sqrt(m / LOOKS[-1])) for m in LOOKS],
        "stay_like_mean_entry_age": {"even 10-14": 12.0, "population 10-16": float(pop_ages.mean())},
        "stay_like_risk_ratio_that_adds_one_year": {
            "even 10-14": risk_cut_for_one_year(uk_a, uk_b, np.linspace(10.0, 14.0, 9)),
            "population 10-16": risk_cut_for_one_year(uk_a, uk_b, np.quantile(pop_ages, np.linspace(0.05, 0.95, 10)))},
        "stay_like_years_gained_by_20_percent_cut": {
            "even 10-14": float(np.mean([fp.years_gained(x, 0.8, 1.0, uk_a, uk_b) for x in np.linspace(10, 14, 9)]))},
    }
    out = dict(
        note="Simulated trials. Part A and C use the rapamycin trial's published level, slope and entry ages. "
             "Part B uses the UK fit and assumed entry ages. 'any_difference_detected' is a two-sided 5% log-rank "
             "test; 'slope_change_flagged' is experiment 001's test; 'lasting_benefit_reported' is experiment 006's.",
        triad_level_per_year=TRIAD_A, triad_slope_per_year=TRIAD_B, uk_level_per_year=uk_a, uk_slope_per_year=uk_b,
        seed=SEED, trials_per_case=nsim, side_calculations=side, cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        extra = (f" | lasting {r['lasting_benefit_reported']:.2f} | early looks {r['any_difference_with_early_looks']:.2f}"
                 if r["part"] == "A" else "")
        print(f"{r['part']} {r['design'][:44]:44s} | {r['effect'][:44]:44s} | placebo dead {r['placebo_share_dead']:.2f}"
              f" | any {r['any_difference_detected']:.3f} | slope {r['slope_change_flagged']:.3f}"
              f" (slower {r['flagged_as_slower_aging']:.3f}) | s-hat {r['slope_estimate_median']:.2f}{extra}")
    print(json.dumps(side, indent=1))


if __name__ == "__main__":
    main()
