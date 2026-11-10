"""Shared helpers for the dog_lifespan stream.

Mortality models (x = age in years):
  Gompertz            mu(x) = a * exp(b x)
  Gompertz-Makeham    mu(x) = c + a * exp(b x)
Cumulative hazard     H(x0,x1) = c (x1-x0) + (a/b) (exp(b x1) - exp(b x0))
Survival              S(x) = exp(-H(0,x))
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "figures")
RESULTS = os.path.join(HERE, "results.json")
os.makedirs(DATA, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

LN2 = np.log(2.0)


def read(name):
    return pd.read_csv(os.path.join(DATA, name), comment="#")


def cumhaz(x0, x1, a, b, c=0.0):
    return c * (x1 - x0) + (a / b) * (np.exp(b * x1) - np.exp(b * x0))


def hazard(x, a, b, c=0.0):
    return c + a * np.exp(b * np.asarray(x, float))


def surv(x, a, b, c=0.0, x0=0.0):
    return np.exp(-cumhaz(x0, np.asarray(x, float), a, b, c))


def life_expectancy(a, b, c=0.0, x0=0.0, upper=200.0):
    """Remaining life expectancy at x0: integral of S(x)/S(x0) from x0 to infinity."""
    f = lambda x: np.exp(-cumhaz(x0, x, a, b, c))
    val, _ = integrate.quad(f, x0, x0 + upper, limit=400)
    return val


def median_age(a, b, c=0.0, x0=0.0):
    from scipy.optimize import brentq
    return brentq(lambda x: cumhaz(x0, x, a, b, c) - LN2, x0, x0 + 300)


def update_results(section, payload):
    res = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as fh:
            res = json.load(fh)
    res[section] = payload
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=2, default=float)


def load_results():
    with open(RESULTS) as fh:
        return json.load(fh)


def save_csv(df, name, header_lines=()):
    path = os.path.join(DATA, name)
    with open(path, "w") as fh:
        for h in header_lines:
            fh.write("# " + h + "\n")
        df.to_csv(fh, index=False)
    return path
