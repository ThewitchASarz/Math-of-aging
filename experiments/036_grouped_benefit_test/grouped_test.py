"""Experiment 036: is a grouped test of "does the benefit differ by frailty score?" sound, and how sharp is it?

Model that makes the data (experiments 032, 034 and 035, unchanged):
    hazard(t) = z * a * exp(b * x0) * r * exp(b * t)
  z     hidden personal risk factor: gamma, mean 1 and variance V = 0.3 among dogs alive at age 2 (ASSUMPTION).
  r     treatment's cut in risk for that dog; 1 in the control group. No dog's rate of aging changes.
  a, b  UK fit, dog_lifespan/results.json (read through experiment 001).
Trial: 1,300 dogs, half treated, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

Frailty scores (ASSUMPTION), experiment 034's make_score, imported:
  rank score   w = RHO * q + sqrt(1 - RHO^2) * noise, q = bell-curve number with the dog's rank in z among
               dogs of its entry age
  risk score   the same with q replaced by ln z, put on a standard scale across the dogs in the trial

True effects (sizes are ASSUMPTIONS). p = the dog's rank (0 to 1) in hidden risk among dogs of its entry age.
  halves   r = r_sturdy if p <= 0.5, else r_frail
  middle   r = 0.3 if 0.25 < p <= 0.75, else 0.7

Tests of "does the benefit differ by score?". Every model is
    ln hazard = ln a + c * (x0 - 9.5) + X . beta + b * t          (no change in the rate of aging anywhere)
and each test is a 5% likelihood-ratio test of a smaller X against a larger one.  G_k = 1 if the score is in
fifth k of all dogs' scores, T_k the same for thirds.
  straight-line        [w, arm]                  against  [w, arm, arm * w]                     1 extra setting
  five, treated only   [w, arm]                  against  [w, arm * G_1, ..., arm * G_5]        4 extra
  five, both arms      [w, G_2..G_5, arm]        against  [w, G_2..G_5, arm * G_1..arm * G_5]   4 extra
  three, both arms     [w, T_2, T_3, arm]        against  [w, T_2, T_3, arm * T_1..arm * T_3]   2 extra
Log-likelihood (experiment 034's, imported): sum of died * (A + b t) - exp(A) * (exp(b t) - 1) / b.

Run:  python3 grouped_test.py            (about 5 minutes on 2 cores)
      python3 grouped_test.py --quick
"""
import importlib.util
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from multiprocessing import Pool

import numpy as np
from scipy.stats import chi2, norm
from scipy.stats import gamma as gamma_dist

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, *parts):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", *parts))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = _load("fingerprint", "001_trial_fingerprint", "fingerprint.py")
e32 = _load("frailty_score", "032_frailty_score", "frailty_score.py")
e34 = _load("entry_age", "034_entry_age_setting", "entry_age.py")

SEED = 20261136
V, AGE0, MID = e32.V, e32.AGE0, e34.MID     # spread 0.3: ASSUMPTION
N, LO, HI, FOLLOW = e32.N, e32.LO, e32.HI, e32.FOLLOW   # entry ages 7 to 12: ASSUMPTION
SCORES = (("rank", 0.5), ("rank", 0.7), ("rank", 0.9), ("rank", 1.0), ("risk", 0.7), ("risk", 1.0))  # ASSUMPTION
# name: (pattern, r for sturdier dogs or middle dogs, r for the others)   all sizes: ASSUMPTION
CASES = {
    "no effect": ("halves", 1.0, 1.0),
    "risk cut 50% in every dog": ("halves", 0.5, 0.5),
    "sturdy helped more (70% and 30% cuts)": ("halves", 0.3, 0.7),
    "sturdy helped a little more (60% and 40% cuts)": ("halves", 0.4, 0.6),
    "middle helped more (70% middle half, 30% both ends)": ("middle", 0.3, 0.7),
}
TESTS = ("straight-line", "five groups, treated only", "five groups, both arms", "three groups, both arms")
EXTRA = dict(zip(TESTS, (1, 4, 4, 2)))


def simulate(pattern, r_in, r_out, a, b, rng):
    """Experiment 032's simulate() with the same random draws in the same order, and a choice of pattern."""
    x0 = rng.uniform(LO, HI, N)
    arm = (np.arange(N) % 2).astype(float)
    H0 = (a / b) * (np.exp(b * x0) - np.exp(b * AGE0))
    u = rng.gamma(1.0 / V, 1.0, N)
    z = u * V / (1.0 + V * H0)
    p = np.clip(gamma_dist.cdf(u, 1.0 / V), 1e-12, 1 - 1e-12)
    q = norm.ppf(p)
    if pattern == "halves":                 # r_in for the sturdier half
        r_dog = np.where(p > 0.5, r_out, r_in)
    else:                                   # r_in for the middle half of the ranking
        r_dog = np.where((p > 0.25) & (p <= 0.75), r_in, r_out)
    r = np.where(arm == 1, r_dog, 1.0)
    t = np.log1p(-np.log(rng.random(N)) * b / (z * a * np.exp(b * x0) * r)) / b
    died = (t <= FOLLOW).astype(float)
    noise = rng.standard_normal(N)
    return (x0, arm, np.minimum(t, FOLLOW), died), q, noise


def fit_level(data, X, start_beta, a, b):
    """Best log-likelihood of the model with columns X and no change in the rate of aging."""
    start = [np.log(a) + b * MID, np.log(b), *start_beta, 0.0]
    return e34.fit(data, X, False, start, False, a, b)


def groups(w, k):
    cuts = np.quantile(w, np.arange(1, k) / k)
    idx = np.searchsorted(cuts, w)
    return [(idx == j).astype(float) for j in range(k)]


def analyse(data, w, a, b):
    """The four likelihood-ratio statistics for 'the benefit differs by score'."""
    x0, arm, t, died = data
    xc = x0 - MID
    G, T = groups(w, 5), groups(w, 3)
    out = {}
    same, th = fit_level(data, np.column_stack([xc, w, arm]), [b, 0.0, 0.0], a, b)
    c, g, lr = th[2], th[3], th[4]
    line, _ = fit_level(data, np.column_stack([xc, w, arm, arm * w]), [c, g, lr, 0.0], a, b)
    out[TESTS[0]] = 2 * (line - same)
    five_t, _ = fit_level(data, np.column_stack([xc, w] + [arm * Gk for Gk in G]), [c, g] + [lr] * 5, a, b)
    out[TESTS[1]] = 2 * (five_t - same)
    for name, K in ((TESTS[2], G), (TESTS[3], T)):
        k = len(K)
        base, th = fit_level(data, np.column_stack([xc, w] + K[1:] + [arm]), [c, g] + [0.0] * (k - 1) + [lr], a, b)
        full, _ = fit_level(data, np.column_stack([xc, w] + K[1:] + [arm * Kj for Kj in K]),
                            list(th[2:3 + k]) + [th[3 + k]] * k, a, b)
        out[name] = 2 * (full - base)
    return out


def check_same_draws(a, b):
    """With the 'halves' pattern, simulate() here must give exactly experiment 032's simulated trial."""
    mine = simulate("halves", 0.3, 0.7, a, b, np.random.default_rng(7))
    theirs = e32.simulate(0.3, 0.7, 1.0, a, b, np.random.default_rng(7))
    return float(max(np.max(np.abs(m - t)) for m, t in zip((*mine[0], mine[1], mine[2]),
                                                           (*theirs[0], theirs[1], theirs[2]))))


def one_chunk(args):
    name, nsim, a, b, seed = args
    pattern, r_in, r_out = CASES[name]
    rng = np.random.default_rng(seed)
    crit = {tname: chi2.ppf(0.95, EXTRA[tname]) for tname in TESTS}
    fired = {(sc, tname): 0 for sc in SCORES for tname in TESTS}
    negative = 0
    for _ in range(nsim):
        data, q, noise = simulate(pattern, r_in, r_out, a, b, rng)
        for sc in SCORES:
            w = e34.make_score(sc[0], sc[1], data, q, noise, a, b)
            for tname, stat in analyse(data, w, a, b).items():
                fired[(sc, tname)] += int(stat > crit[tname])
                negative += int(stat < -1e-4)
    return name, nsim, fired, negative


def main():
    quick = "--quick" in sys.argv
    a, b = fp.dog_parameters()
    same_draws = check_same_draws(a, b)
    assert same_draws < 1e-12, same_draws
    chunks, per = (1, 6) if quick else (4, 150)
    jobs = [(name, per, a, b, SEED + 1000 * i + k) for i, name in enumerate(CASES) for k in range(chunks)]
    with Pool(2) as pool:
        parts = pool.map(one_chunk, jobs, chunksize=1)
    rows = []
    for name, (pattern, r_in, r_out) in CASES.items():
        mine = [p for p in parts if p[0] == name]
        n = sum(p[1] for p in mine)
        for sc in SCORES:
            row = dict(case=name, pattern=pattern, risk_multiplier_sturdier_or_middle=r_in,
                       risk_multiplier_others=r_out, score_kind=sc[0], score_quality=sc[1], n_trials=n)
            for tname in TESTS:
                row["fired: " + tname] = sum(p[2][(sc, tname)] for p in mine) / n
            rows.append(row)
    out = dict(note="Simulated 10-year trials of 1,300 dogs, half treated from enrolment. Each 'fired' figure is the "
                    "share of trials in which a 5% likelihood-ratio test says the benefit differs by frailty "
                    "score. No dog's rate of aging changes in any case. ASSUMPTIONS: spread 0.3, entry ages 7 to "
                    "12, the form and quality of the scores, and every size and pattern of benefit.",
               seed=SEED, uk_a_per_year=a, uk_b_per_year=b, spread=V, n_dogs=N, entry_ages=[LO, HI], years=FOLLOW,
               extra_settings=EXTRA, max_difference_from_experiment_032_simulation=same_draws,
               fits_where_larger_model_fitted_worse=int(sum(p[3] for p in parts)), cases=rows)
    with open(os.path.join(HERE, "results_quick.json" if quick else "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['case'][:28]:28s} | {r['score_kind']} {r['score_quality']:.1f} | "
              + " | ".join(f"{r['fired: ' + tname]:.3f}" for tname in TESTS))
    print("larger model fitted worse:", out["fits_where_larger_model_fitted_worse"])


if __name__ == "__main__":
    main()
