"""Experiment 002: fit the saturating-removal ("worn-out-cell") model of Karin et al. 2019 to dog life tables.

Model: senescent-cell load X follows  dX = [eta*t - beta*X/(kappa+X)] dt + sqrt(2*eps) dW, death when X
reaches Xc. X settles in days to weeks while production eta*t changes over years, so the hazard at age t is
taken as 1 / (mean first-passage time from 0 to Xc with production p = eta*t held fixed).

Human reference settings (published, Karin et al. 2019, Fig. 4d):
    eta = 0.00135 per day per year, beta = 0.15 per day, kappa = 0.5, eps = 0.142 per day, Xc = 17.

Question: which setting must differ for dogs? Fits (kappa held at the human value throughout):
    A  eta only            B  eta and beta           C  eta and eps           D  eta and Xc
    E  eta, beta and eps   (extra; not part of the written prediction)
    B_bg, C_bg, D_bg       (extra, added after the first results: the same fits with a constant background
                            death rate of 0.005 per year, the value fitted to both dog tables)

Run:  python3 dog_cell_model.py            (about 2 minutes)
      python3 dog_cell_model.py --quick    (no bootstrap; for tests)
"""
import json
import os
import sys

# One thread per process: the linear-algebra library otherwise starts a thread pool for tiny arrays,
# which made this script about 30 times slower on a 4-core test machine.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
SEED = 20261110
HUMAN = dict(eta=0.00135, beta=0.15, kappa=0.5, eps=0.142, Xc=17.0)
BACKGROUND = 0.005     # per year: age-independent dog death rate from the Gompertz-Makeham fits (both tables)
# name: (settings allowed to differ from the human values, background death rate added to the model)
FITS = {"A": (("eta",), 0.0), "B": (("eta", "beta"), 0.0), "C": (("eta", "eps"), 0.0),
        "D": (("eta", "Xc"), 0.0), "E": (("eta", "beta", "eps"), 0.0),
        # extra fits, added after seeing the first results; not part of the written prediction
        "B_bg": (("eta", "beta"), BACKGROUND), "C_bg": (("eta", "eps"), BACKGROUND),
        "D_bg": (("eta", "Xc"), BACKGROUND)}
OFFSETS = np.array([0.1, 0.3, 0.5, 0.7, 0.9])     # points inside each one-year age interval


def potential(x, p, beta, kappa):
    """U with U'(x) = beta*x/(kappa+x) - p."""
    return -p * x + beta * (x - kappa * np.log1p(x / kappa))


def hazard_per_year(ages, eta, beta, kappa, eps, Xc, n=2001):
    """365 / (mean first-passage time in days), production frozen at p = eta * age."""
    p = (eta * np.asarray(ages, float))[:, None]
    x = np.linspace(0.0, Xc, n)[None, :]
    u = potential(x, p, beta, kappa) / eps
    um = u.max(axis=1, keepdims=True)
    inner = cumulative_trapezoid(np.exp(-(u - um)), x[0], axis=1, initial=0.0)
    mfpt = np.trapezoid(np.exp(u - um) * inner, x[0], axis=1) / eps
    return 365.0 / mfpt


def interval_hazard(age_starts, **par):
    """Mean model hazard over each interval [age, age+1)."""
    pts = (np.asarray(age_starts, float)[:, None] + OFFSETS[None, :]).ravel()
    return hazard_per_year(pts, **par).reshape(len(age_starts), len(OFFSETS)).mean(axis=1)


def load_table(name):
    path = os.path.join(ROOT, "dog_lifespan", "data", name)
    rows = [l.strip().split(",") for l in open(path) if not l.startswith("#")]
    cols = rows[0]
    data = {c: np.array([float(r[i]) for r in rows[1:]]) for i, c in enumerate(cols)}
    return data


def observed(table, lo, hi):
    keep = (table["age_start"] >= lo) & (table["age_start"] <= hi)
    ages = table["age_start"][keep]
    haz = -np.log(1.0 - table["qx"][keep])          # exact mean hazard over the year
    return ages, haz, keep


def params_from(theta, free):
    par = dict(HUMAN)
    for name, value in zip(free, theta):
        par[name] = HUMAN[name] * np.exp(value)
    return par


def residuals(theta, free, ages, log_obs, lam=0.0):
    par = params_from(theta, free)
    with np.errstate(all="ignore"):
        h = interval_hazard(ages, **par) + lam
    r = np.log(h) - log_obs
    return np.where(np.isfinite(r), r, 50.0)


STARTS = {"eta": np.log([2.0, 4.0, 8.0, 12.0]), "beta": np.log([0.4, 0.7, 1.0, 2.0]),
          "eps": np.log([0.6, 1.0, 2.0, 3.5]), "Xc": np.log([0.4, 0.7, 1.0, 1.5])}
BOUNDS = {"eta": np.log([0.5, 60.0]), "beta": np.log([0.05, 20.0]), "eps": np.log([0.05, 20.0]),
          "Xc": np.log([0.1, 5.0])}


def fit(free, ages, haz, starts=None, lam=0.0):
    log_obs = np.log(haz)
    lo = [BOUNDS[f][0] for f in free]
    hi = [BOUNDS[f][1] for f in free]
    if starts is None:
        grids = np.meshgrid(*[STARTS[f] for f in free], indexing="ij")
        starts = np.stack([g.ravel() for g in grids], axis=1)
    best = None
    for s0 in starts:
        res = least_squares(residuals, s0, args=(free, ages, log_obs, lam), bounds=(lo, hi), xtol=1e-10, ftol=1e-10)
        if best is None or res.cost < best.cost:
            best = res
    rmse = float(np.sqrt(np.mean(best.fun ** 2)))
    return best.x, rmse


def summarise(theta, free, ages, haz, lam=0.0):
    par = params_from(theta, free)
    model = interval_hazard(ages, **par) + lam
    slope = np.polyfit(ages + 0.5, np.log(model), 1)[0]
    out = {f"{f}_ratio_to_human": float(np.exp(v)) for f, v in zip(free, theta)}
    out["model_doubling_time_years"] = float(np.log(2) / slope)
    out["max_abs_log_error"] = float(np.max(np.abs(np.log(model) - np.log(haz))))
    # what a lab could measure: typical load and clearance half-life at a young-adult age (dog age 3)
    p = par["eta"] * 3.0
    if p < par["beta"]:
        xs = par["kappa"] * p / (par["beta"] - p)
        out["young_adult_load_at_age_3"] = float(xs)
        out["clearance_half_life_days_at_age_3"] = float(np.log(2) * (par["kappa"] + xs) ** 2
                                                         / (par["beta"] * par["kappa"]))
        out["production_as_share_of_clearance_capacity_at_age_3"] = float(p / par["beta"])
    return out


def human_reference():
    """Same quantities for the published human settings at age 30."""
    p = HUMAN["eta"] * 30.0
    xs = HUMAN["kappa"] * p / (HUMAN["beta"] - p)
    return dict(young_adult_load_at_age_30=float(xs),
                clearance_half_life_days_at_age_30=float(np.log(2) * (HUMAN["kappa"] + xs) ** 2
                                                         / (HUMAN["beta"] * HUMAN["kappa"])),
                production_as_share_of_clearance_capacity_at_age_30=float(p / HUMAN["beta"]))


def main():
    quick = "--quick" in sys.argv
    tables = {"UK (Teng 2022), ages 2-17": ("teng2022_uk_dog_life_table.csv", 2, 17),
              "US (Montoya 2023), ages 2-16": ("montoya2023_us_dog_life_table.csv", 2, 16)}
    out = dict(note="Least-squares fits of the saturating-removal model to dog life tables, on the log of the "
                    "yearly death rate. kappa is fixed at the human value. Ratios are dog setting / published "
                    "human setting.", human_settings=HUMAN, human_reference=human_reference(), seed=SEED,
               tables={})
    for label, (fname, lo, hi) in tables.items():
        tab = load_table(fname)
        ages, haz, keep = observed(tab, lo, hi)
        res = {}
        for name, (free, lam) in FITS.items():
            theta, rmse = fit(free, ages, haz, lam=lam)
            res[name] = dict(free_settings=list(free), background_death_rate_per_year=lam,
                             rmse_log_death_rate=rmse, **summarise(theta, free, ages, haz, lam))
            res[name]["_theta"] = [float(v) for v in theta]
            print(f"{label[:14]:14s} fit {name} {str(free):24s} rmse {rmse:.3f}  "
                  + "  ".join(f"{f}x{np.exp(v):.2f}" for f, v in zip(free, theta)))
        out["tables"][label] = dict(ages=[float(a) for a in ages], observed_hazard=[float(h) for h in haz], fits=res)

        # sampling uncertainty for the UK table, which publishes death counts
        if "dx" in tab and not quick:
            rng = np.random.default_rng(SEED)
            lx, dx = tab["lx"][keep], tab["dx"][keep]
            boot = {name: [] for name in ("B", "C")}
            for _ in range(200):
                d_star = rng.binomial(lx.astype(int), dx / lx)
                h_star = -np.log(1.0 - d_star / lx)
                for name in boot:
                    free = FITS[name][0]
                    th, rm = fit(free, ages, h_star, starts=[res[name]["_theta"]])
                    boot[name].append(list(np.exp(th)) + [rm])
            for name, vals in boot.items():
                arr = np.array(vals)
                ci = np.percentile(arr, [2.5, 97.5], axis=0)
                res[name]["bootstrap_95ci"] = {f"{f}_ratio_to_human": [float(ci[0, i]), float(ci[1, i])]
                                               for i, f in enumerate(FITS[name][0])}
                res[name]["bootstrap_95ci"]["rmse_log_death_rate"] = [float(ci[0, -1]), float(ci[1, -1])]
                res[name]["bootstrap_note"] = ("Sampling error only (binomial resampling of the published death "
                                               "counts, 200 replicates). It does not include model error.")
    for tab in out["tables"].values():
        for r in tab["fits"].values():
            r.pop("_theta", None)
    name = "results_quick.json" if quick else "results.json"
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
