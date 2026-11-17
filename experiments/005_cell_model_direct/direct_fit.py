"""Experiment 005: fit the worn-out-cell model to dog life tables by solving the model directly.

Experiment 002 used a shortcut: the death rate at each age was computed as if production of worn-out cells
were frozen at that age's value. Here the model itself is solved, with production rising continuously.

Model (Karin et al. 2019), time t in days, age in years:
    dX = [eta*age - beta*X/(kappa+X)] dt + sqrt(2*eps) dW,   X(0) = 0, X cannot go below 0, death when X = Xc.

Direct solution: the spread of X across dogs, P(x, t), obeys the Fokker-Planck equation
    dP/dt = -d/dx[(eta*age - beta*x/(kappa+x)) P] + eps * d2P/dx2
with no flow through x = 0 and P = 0 at x = Xc. It is solved on a grid of N cells (finite volumes with
exponentially fitted fluxes, implicit time steps). The share alive is the total of P. The death rate in the
year of age a is -ln(S(a+1)/S(a)), the same definition used for the life tables.

Inputs: the published human settings and the two dog life tables, exactly as in experiment 002. There are no
assumed inputs in this experiment other than the model itself and kappa held at the human value (ASSUMPTION,
carried over from experiment 002: death records cannot pin down more than about two settings).

Run:  python3 direct_fit.py            (about 3 minutes on 2 cores)
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from multiprocessing import Pool
from scipy.linalg import solve_banded
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(HERE, "..", "002_dog_cell_model"))
from dog_cell_model import HUMAN, BOUNDS, interval_hazard, load_table, observed   # noqa: E402
from sde_check import simulate_death_ages                                          # noqa: E402

SEED = 20261117
N_CELLS = 240          # grid cells between 0 and Xc
DT_DAYS = 1.0          # time step
FITS = {"A": ("eta",), "B": ("eta", "beta"), "C": ("eta", "eps"), "D": ("eta", "Xc")}
TABLES = {"UK (Teng 2022), ages 2-17": ("teng2022_uk_dog_life_table.csv", 2, 17),
          "US (Montoya 2023), ages 2-16": ("montoya2023_us_dog_life_table.csv", 2, 16)}
# Starting points for the search, as multiples of the human settings. The first of each list is close to
# experiment 002's answer; the others are spread out to catch a different best fit.
STARTS = {"A": [(10.0,), (4.0,), (20.0,)],
          "B": [(3.6, 0.52), (6.0, 1.0), (2.0, 0.3)],
          "C": [(6.0, 1.8), (3.0, 1.0), (9.0, 3.0)],
          "D": [(5.0, 0.6), (8.0, 1.0), (3.0, 0.4)]}


def bernoulli(z):
    """z / (exp(z) - 1), safe near 0 and for large |z|."""
    z = np.clip(z, -500.0, 500.0)
    small = np.abs(z) < 1e-8
    return np.where(small, 1.0 - 0.5 * z, z / np.expm1(np.where(small, 1.0, z)))


def survival(eta, beta, kappa, eps, Xc, max_age, n=N_CELLS, dt=DT_DAYS, return_density=False):
    """Share of dogs alive at ages 0, 1, ..., max_age (whole years), production rising continuously."""
    h = Xc / n
    xi = h * np.arange(1, n)                    # interfaces between cells
    removal = beta * xi / (kappa + xi)
    removal_wall = beta * Xc / (kappa + Xc)
    m = np.zeros(n)
    m[0] = 1.0                                  # every dog starts with no worn-out cells
    steps_per_year = int(round(365.0 / dt))
    surv = [1.0]
    ab = np.zeros((3, n))
    for year in range(int(max_age)):
        for k in range(steps_per_year):
            age = year + (k + 0.5) / steps_per_year          # production at the middle of the step
            p = eta * age
            v = p - removal
            pe = v * h / eps
            out_right = (eps / h) * bernoulli(-pe)           # flow from cell i to i+1, per unit in cell i
            out_left = (eps / h) * bernoulli(pe)             # flow from cell i+1 to i, per unit in cell i+1
            vw = p - removal_wall
            wall = (eps / (h / 2)) * bernoulli(-vw * (h / 2) / eps)   # flow out through x = Xc (death)
            # mass balance, backward Euler: (I - dt/h * L) m_new = m_old
            c = dt / h
            diag = np.ones(n)
            diag[:-1] += c * out_right
            diag[1:] += c * out_left
            diag[-1] += c * wall
            ab[0, 1:] = -c * out_left
            ab[1] = diag
            ab[2, :-1] = -c * out_right
            m = solve_banded((1, 1), ab, m, overwrite_b=True, check_finite=False)
        surv.append(float(m.sum()))
    s = np.array(surv)
    return (s, m) if return_density else s


def direct_hazard(age_starts, **par):
    """Death rate in each year of age [a, a+1) from the direct solution."""
    age_starts = np.asarray(age_starts, int)
    s = survival(max_age=int(age_starts.max()) + 1, **par)
    with np.errstate(all="ignore"):
        return -np.log(s[age_starts + 1] / s[age_starts])


def params_from(theta, free):
    par = dict(HUMAN)
    for name, value in zip(free, theta):
        par[name] = HUMAN[name] * np.exp(value)
    return par


def residuals(theta, free, ages, log_obs):
    with np.errstate(all="ignore"):
        r = np.log(direct_hazard(ages, **params_from(theta, free))) - log_obs
    return np.where(np.isfinite(r), r, 50.0)


def fit_one(job):
    label, name = job[:2]
    fname, lo, hi = TABLES[label]
    if len(job) == 4:              # extra fit on a narrower age range, one starting point
        lo, hi = job[2], job[3]
    ages, haz, _ = observed(load_table(fname), lo, hi)
    free = FITS[name]
    starts = STARTS[name][:1] if len(job) == 4 else STARTS[name]
    log_obs = np.log(haz)
    lob = [BOUNDS[f][0] for f in free]
    hib = [BOUNDS[f][1] for f in free]
    best, all_costs = None, []
    for s0 in starts:
        res = least_squares(residuals, np.log(s0), args=(free, ages, log_obs), bounds=(lob, hib),
                            xtol=1e-7, ftol=1e-9, diff_step=1e-4)
        all_costs.append(float(np.sqrt(np.mean(res.fun ** 2))))
        if best is None or res.cost < best.cost:
            best = res
    par = params_from(best.x, free)
    model = direct_hazard(ages, **par)
    shortcut = interval_hazard(ages.astype(float), **par)
    gap = np.log(model / shortcut)
    out = dict(free_settings=list(free), rmse_log_death_rate=float(np.sqrt(np.mean(best.fun ** 2))),
               rmse_from_each_start=all_costs,
               max_abs_log_error=float(np.max(np.abs(best.fun))),
               model_death_rate=[float(v) for v in model],
               shortcut_at_these_settings=dict(
                   rmse_log_gap_direct_vs_shortcut=float(np.sqrt(np.mean(gap ** 2))),
                   max_abs_log_gap=float(np.max(np.abs(gap))),
                   shortcut_rmse_against_table=float(np.sqrt(np.mean((np.log(shortcut) - log_obs) ** 2)))))
    for f, v in zip(free, best.x):
        out[f"{f}_ratio_to_human"] = float(np.exp(v))
    return label, name, out


def simulate_with_crossing_check(eta, beta, kappa, eps, Xc, n, dt, max_age, seed):
    """Brute-force simulation, written separately from the solver. Between two time points a path can touch
    the lethal load and come back; the chance of that, exp(-(Xc-x0)(Xc-x1)/(eps*dt)), is used to catch those
    deaths (the standard Brownian-bridge correction)."""
    rng = np.random.default_rng(seed)
    x = np.zeros(n)
    death_age = np.full(n, np.nan)
    idx = np.arange(n)
    sd = np.sqrt(2.0 * eps * dt)
    for k in range(int(max_age * 365.0 / dt)):
        age = (k + 0.5) * dt / 365.0
        x0 = x[idx]
        x1 = np.abs(x0 + (eta * age - beta * x0 / (kappa + x0)) * dt + sd * rng.standard_normal(idx.size))
        with np.errstate(over="ignore"):
            touched = rng.random(idx.size) < np.exp(-(Xc - x0) * (Xc - x1) / (eps * dt))
        dead = (x1 >= Xc) | touched
        death_age[idx[dead]] = (k + 1) * dt / 365.0
        x[idx] = x1
        idx = idx[~dead]
        if idx.size == 0:
            break
    return death_age


def yearly_rates(death_age, ages):
    rows = []
    for a in ages:
        at_risk = int(np.sum(~(death_age < a)))
        died = int(np.sum((death_age >= a) & (death_age < a + 1)))
        if at_risk < 200 or died < 20:
            continue
        rows.append((int(a), -np.log(1.0 - died / at_risk), died, np.sqrt((1.0 - died / at_risk) / died)))
    return rows


def solver_guard(exp2):
    """Direct solution against brute-force simulation, experiment 002's fit-B settings for the US table."""
    fitb = exp2["tables"]["US (Montoya 2023), ages 2-16"]["fits"]["B"]
    par = dict(HUMAN)
    par["eta"] = HUMAN["eta"] * fitb["eta_ratio_to_human"]
    par["beta"] = HUMAN["beta"] * fitb["beta_ratio_to_human"]
    ages = np.arange(2, 17)
    direct = dict(zip(ages.tolist(), direct_hazard(ages, **par)))
    fine = -np.diff(np.log(survival(max_age=17, n=2 * N_CELLS, dt=DT_DAYS / 4, **par)))[2:]
    grid_gap = float(np.max(np.abs(np.log(np.array(list(direct.values())) / fine))))
    sim = yearly_rates(simulate_with_crossing_check(n=40000, dt=0.25, max_age=17.0, seed=SEED, **par), ages)
    rows = [dict(age=a, direct=float(direct[a]), simulated=float(h), simulated_over_direct=float(h / direct[a]),
                 deaths=d, z=float(np.log(h / direct[a]) / se)) for a, h, d, se in sim]
    lg = np.log([r["simulated_over_direct"] for r in rows])
    # the simulation of experiment 002's sde_check.py, unchanged: it looks for death only at the end of each
    # half-day step, so it misses some deaths
    old = yearly_rates(simulate_death_ages(n=40000, dt=0.5, seed=SEED, **par), ages)
    lg_old = np.log([h / direct[a] for a, h, d, se in old])
    return dict(note="Direct solution against a brute-force simulation of 40,000 dogs (time step 0.25 day, with "
                     "a check for crossings between time points), experiment 002 fit-B settings, US table. "
                     "z = log(simulated/direct) / its rough standard error; values within about 2 to 3 are "
                     "simulation noise across 15 ages.",
                settings=par, by_age=rows,
                mean_log_ratio=float(lg.mean()), rmse_log_ratio=float(np.sqrt(np.mean(lg ** 2))),
                max_abs_z=float(np.max(np.abs([r["z"] for r in rows]))),
                mean_z=float(np.mean([r["z"] for r in rows])),
                max_abs_log_gap_grid_refinement=grid_gap,
                experiment_002_simulation=dict(
                    note="The simulation in experiment 002's sde_check.py (half-day steps, death checked only "
                         "at the end of each step) against the direct solution. A negative mean means that "
                         "simulation reads low.",
                    mean_log_ratio=float(lg_old.mean()), rmse_log_ratio=float(np.sqrt(np.mean(lg_old ** 2)))))


def human_check():
    """Extra, not part of the prediction: how good is the shortcut at the published human settings?"""
    ages = np.arange(40, 90, 5)
    s = survival(max_age=90, n=N_CELLS, dt=DT_DAYS, **HUMAN)
    with np.errstate(all="ignore"):
        direct = -np.log(s[ages + 1] / s[ages])
    shortcut = interval_hazard(ages.astype(float), **HUMAN)
    gap = np.log(direct / shortcut)
    return dict(note="Extra, not a test of the predictions. Published human settings, ages 40 to 85: "
                     "direct solution against the frozen-production shortcut.",
                ages=[int(a) for a in ages], direct=[float(v) for v in direct],
                shortcut=[float(v) for v in shortcut],
                direct_over_shortcut=[float(v) for v in np.exp(gap)],
                max_abs_log_gap=float(np.max(np.abs(gap))))


def main():
    t0 = time.time()
    with open(os.path.join(HERE, "..", "002_dog_cell_model", "results.json")) as fh:
        exp2 = json.load(fh)
    jobs = [(label, name) for label in TABLES for name in FITS]
    # Extra, added after the independent check and so not a test of the predictions: the two-setting fits
    # again on ages 5 to 15 only, leaving out the youngest ages, where the model misses most.
    extra_jobs = [(label, name, 5, 15) for label in TABLES for name in ("B", "C", "D")]
    with Pool(2) as pool:
        both = pool.map(fit_one, jobs + extra_jobs, chunksize=1)
    fitted, extra = both[:len(jobs)], both[len(jobs):]
    out = dict(note="Fits of the saturating-removal model to dog life tables with the model solved directly "
                    "(production rising continuously). Least squares on the log of the yearly death rate. kappa "
                    "is fixed at the human value. Ratios are dog setting / published human setting. "
                    "'experiment_002' repeats the shortcut-based values for comparison.",
               human_settings=HUMAN, grid_cells=N_CELLS, time_step_days=DT_DAYS, seed=SEED, tables={})
    for label in TABLES:
        fname, lo, hi = TABLES[label]
        ages, haz, _ = observed(load_table(fname), lo, hi)
        out["tables"][label] = dict(ages=[float(a) for a in ages], observed_hazard=[float(h) for h in haz], fits={})
    for label, name, res in fitted:
        old = exp2["tables"][label]["fits"][name]
        res["experiment_002"] = {k: old[k] for k in old if k.endswith("_ratio_to_human") or k == "rmse_log_death_rate"}
        res["rmse_change_from_experiment_002"] = res["rmse_log_death_rate"] - old["rmse_log_death_rate"]
        out["tables"][label]["fits"][name] = res
        print(f"{label[:14]:14s} fit {name} rmse {res['rmse_log_death_rate']:.3f} (002: {old['rmse_log_death_rate']:.3f})  "
              + "  ".join(f"{f}x{res[f + '_ratio_to_human']:.2f}" for f in FITS[name])
              + f"   starts {['%.3f' % c for c in res['rmse_from_each_start']]}")
    for label in TABLES:
        f = out["tables"][label]["fits"]
        errs = {k: f[k]["rmse_log_death_rate"] for k in ("B", "C", "D")}
        out["tables"][label]["order_best_to_worst_BCD"] = sorted(errs, key=errs.get)
        out["tables"][label]["spread_of_BCD_errors"] = max(errs.values()) - min(errs.values())
        out["tables"][label]["residuals_fit_B_log_model_minus_log_observed"] = [
            float(np.log(m / o)) for m, o in zip(f["B"]["model_death_rate"], out["tables"][label]["observed_hazard"])]
        ex = {name: res for lab, name, res in extra if lab == label}
        out["tables"][label]["extra_ages_5_to_15"] = dict(
            note="Added after the independent check; not a test of the predictions. One starting point.",
            fits={name: {k: v for k, v in r.items() if k.endswith("_ratio_to_human") or k == "rmse_log_death_rate"}
                  for name, r in ex.items()},
            order_best_to_worst_BCD=sorted(ex, key=lambda k: ex[k]["rmse_log_death_rate"]))
        print(label[:14], "ages 5-15:", {k: round(v["rmse_log_death_rate"], 3) for k, v in ex.items()})
    out["solver_guard"] = solver_guard(exp2)
    out["human_check"] = human_check()
    print("guard: rmse log ratio", round(out["solver_guard"]["rmse_log_ratio"], 4), "max |z|",
          round(out["solver_guard"]["max_abs_z"], 2), "grid gap",
          round(out["solver_guard"]["max_abs_log_gap_grid_refinement"], 4))
    print("human: max gap", round(out["human_check"]["max_abs_log_gap"], 4))
    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("seconds", round(time.time() - t0))


if __name__ == "__main__":
    main()
