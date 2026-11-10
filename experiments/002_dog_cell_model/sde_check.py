"""Check of the shortcut used in dog_cell_model.py.

dog_cell_model.py computes the death rate at each age as if production were frozen at that age's value
(the "quasi-static" shortcut). This script simulates the model itself, with production rising continuously,
and compares the simulated death rates with the shortcut. Added after the independent check found the
shortcut is off by up to about 20 percent for dogs.

Run:  python3 sde_check.py      (about 1 minute)
"""
import json
import os

import numpy as np

from dog_cell_model import HUMAN, interval_hazard

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261110


def simulate_death_ages(eta, beta, kappa, eps, Xc, n=40000, dt=0.5, max_age=18.0, seed=SEED):
    """Euler-Maruyama for dX = [eta*t - beta*X/(kappa+X)]dt + sqrt(2 eps) dW, reflecting at 0, death at Xc.
    Time step dt in days; age t in years."""
    rng = np.random.default_rng(seed)
    x = np.zeros(n)
    alive = np.ones(n, bool)
    death_age = np.full(n, np.nan)
    sd = np.sqrt(2.0 * eps * dt)
    steps = int(max_age * 365.0 / dt)
    for k in range(steps):
        age = k * dt / 365.0
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        xa = x[idx]
        xa = xa + (eta * age - beta * xa / (kappa + xa)) * dt + sd * rng.standard_normal(idx.size)
        xa = np.abs(xa)                                  # reflecting boundary at 0
        dead = xa >= Xc
        death_age[idx[dead]] = age
        alive[idx[dead]] = False
        x[idx] = xa
    return death_age


def main():
    with open(os.path.join(HERE, "results.json")) as fh:
        res = json.load(fh)
    fit = res["tables"]["US (Montoya 2023), ages 2-16"]["fits"]["B"]
    par = dict(HUMAN)
    par["eta"] = HUMAN["eta"] * fit["eta_ratio_to_human"]
    par["beta"] = HUMAN["beta"] * fit["beta_ratio_to_human"]
    death_age = simulate_death_ages(**par)
    ages = np.arange(2, 17)
    rows = []
    for a in ages:
        at_risk = np.sum(~(death_age < a))               # alive at the start of the interval
        died = np.sum((death_age >= a) & (death_age < a + 1))
        if at_risk < 200 or died < 20:
            continue
        simulated = -np.log(1.0 - died / at_risk)
        shortcut = float(interval_hazard(np.array([float(a)]), **par)[0])
        rows.append(dict(age=int(a), simulated_death_rate=float(simulated), shortcut_death_rate=shortcut,
                         ratio=float(simulated / shortcut), deaths=int(died)))
    logs = np.log([r["ratio"] for r in rows])
    out = dict(note="Simulated model against the frozen-production shortcut, fit-B settings for the US dog table. "
                    "ratio = simulated / shortcut.",
               settings=par, dogs_simulated=int(death_age.size), time_step_days=0.5, by_age=rows,
               max_abs_log_error=float(np.max(np.abs(logs))), rmse_log_error=float(np.sqrt(np.mean(logs ** 2))))
    with open(os.path.join(HERE, "sde_check.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"age {r['age']:2d}  simulated {r['simulated_death_rate']:.4f}  shortcut {r['shortcut_death_rate']:.4f}"
              f"  ratio {r['ratio']:.2f}")
    print("max abs log error", round(out["max_abs_log_error"], 3), " rmse", round(out["rmse_log_error"], 3))


if __name__ == "__main__":
    main()
