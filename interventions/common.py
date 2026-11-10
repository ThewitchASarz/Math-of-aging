"""Shared helpers for the interventions stream.

Gompertz-Makeham hazard:  mu(x) = c + a * exp(b * x)
  x : age in years
  a : hazard at age 0 of the ageing (Gompertz) component, per year
  b : Gompertz slope, per year (how fast the ageing hazard grows)
  c : Makeham constant, age-independent background hazard, per year

PARAMETERS ARE FITTED, not assumed. They are read from the other parts of this project:
  human: Gompertz fit to the US Social Security 2023 period life table, both sexes, ages 40-94
         (dog_lifespan/results.json -> part1_mortality_fits.human_both_gompertz_40to94)
  dog:   Gompertz fit to the UK VetCompass dog life table, ages 2-19 (Teng et al. 2022)
         (dog_lifespan/results.json -> part1_mortality_fits.RECOMMENDED_dog_parameters)
If that file is missing, the same values are used from the constants below.
Earlier versions of this code used provisional values (human b = 0.087, dog b = 0.25 per year).
"""
import json
import os

import numpy as np
from scipy import integrate, special

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

_FALLBACK = {
    "human": {"a": 5.576363367663802e-05, "b": 0.08511920061726733, "c": 0.0},
    "dog": {"a": 0.004967182098875431, "b": 0.31954873909049325, "c": 0.0},
}


def _load_fitted():
    out = {k: dict(v) for k, v in _FALLBACK.items()}
    path = os.path.join(HERE, "..", "dog_lifespan", "results.json")
    try:
        with open(path) as fh:
            fits = json.load(fh)["part1_mortality_fits"]
        h = fits["human_both_gompertz_40to94"]
        d = fits["RECOMMENDED_dog_parameters"]
        out["human"].update(a=h["a_per_year"], b=h["b_per_year"])
        out["dog"].update(a=d["a_per_year"], b=d["b_per_year"])
    except (OSError, KeyError, ValueError):
        pass
    out["human"]["b_range"] = [0.0825, 0.1014]   # fitted range across sexes and model forms
    out["dog"]["b_range"] = [0.3054, 0.3439]     # fitted range across UK/US tables and model forms
    for sp in out.values():
        sp["mode"] = float(np.log(sp["b"] / sp["a"]) / sp["b"])   # modal age at death of the Gompertz part
    return out


MORTALITY = _load_fitted()


def gm_survival(t, x0, a, b, c=0.0, hr=1.0, s=1.0):
    """Survival from age x0 to x0+t when, from x0 on, the Gompertz part is multiplied by hr
    and its slope is multiplied by s (s<1 = slower ageing). Makeham c is untouched."""
    m0 = a * np.exp(b * x0)
    bs = b * s
    if abs(bs) < 1e-12:
        H = hr * m0 * t
    else:
        H = hr * m0 * (np.exp(bs * t) - 1.0) / bs
    return np.exp(-c * t - H)


def life_expectancy(x0, a, b, c=0.0, hr=1.0, s=1.0, tmax=None):
    """Remaining life expectancy at age x0 (numerical integral of survival)."""
    if tmax is None:
        tmax = 2000.0 if c == 0 else min(2000.0, 60.0 / max(c, 1e-6))
    f = lambda t: gm_survival(t, x0, a, b, c, hr, s)
    # split integral for accuracy
    pts = np.unique(np.concatenate([[0], np.geomspace(1, tmax, 30)]))
    tot = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        tot += integrate.quad(f, lo, hi, limit=200)[0]
    return tot


def le_gompertz_closed(x0, a, b, hr=1.0, s=1.0):
    """Closed form for pure Gompertz (c=0): e = exp(M) E1(M) / (b s), M = hr a e^{b x0} / (b s)."""
    bs = b * s
    M = hr * a * np.exp(b * x0) / bs
    return np.exp(M) * special.exp1(M) / bs


def update_results(key, value):
    path = os.path.join(HERE, "results.json")
    res = {}
    if os.path.exists(path):
        with open(path) as fh:
            res = json.load(fh)
    res[key] = value
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=float)
