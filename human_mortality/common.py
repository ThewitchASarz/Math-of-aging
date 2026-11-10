"""Shared helpers for the human_mortality stream.

Data: SSA 2023 period life table (2026 Trustees Report),
https://www.ssa.gov/oact/STATS/table4c6.html (transcribed into data/ssa_2023_period_life_table.csv).
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "figures")
RESULTS = os.path.join(HERE, "results.json")
SSA_URL = "https://www.ssa.gov/oact/STATS/table4c6.html"
SSA_DESC = "US Social Security Administration period life table, data year 2023 (2026 Trustees Report)"
os.makedirs(FIG, exist_ok=True)
os.makedirs(DATA, exist_ok=True)


def load_ssa():
    """Return DataFrame indexed by integer age 0..119 with qx, lx, ex for male/female,
    plus mu_<sex> = -ln(1-qx): the average force of mortality over [x, x+1) (per year)."""
    d = pd.read_csv(os.path.join(DATA, "ssa_2023_period_life_table.csv"), comment="#")
    d = d.set_index("age")
    for s in ("male", "female"):
        d[f"mu_{s}"] = -np.log1p(-d[f"{s}_qx"])
    return d


def update_results(section, payload):
    """Merge payload into results.json under key `section`."""
    res = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            res = json.load(f)
    res[section] = payload
    with open(RESULTS, "w") as f:
        json.dump(res, f, indent=2, default=float)


# ---------------- hazard models (x = age in years, rates per year) ----------------
def H_gm(x, a, b, c=0.0):
    """Cumulative hazard of Gompertz-Makeham, from age 0: c x + (a/b)(e^{bx}-1)."""
    return c * x + (a / b) * np.expm1(b * x)


def mu_gm(x, a, b, c=0.0):
    return c + a * np.exp(b * x)


def H_ggm(x, a, b, c, s2):
    """Cumulative hazard of gamma-Gompertz-Makeham (population level, from age 0):
    c x + (1/s2) ln(1 + s2 (a/b)(e^{bx}-1))."""
    return c * x + np.log1p(s2 * (a / b) * np.expm1(b * x)) / s2


def mu_ggm(x, a, b, c, s2):
    return c + a * np.exp(b * x) / (1.0 + s2 * (a / b) * np.expm1(b * x))


def H_kannisto(x, a, b, x0=80.0):
    """Kannisto logistic: mu = a e^{b(x-x0)} / (1 + a e^{b(x-x0)}); H from x0 = (1/b) ln((1+a e^{b(x-x0)})/(1+a))."""
    return np.log1p(a * np.exp(b * (x - x0))) / b


def mu_kannisto(x, a, b, x0=80.0):
    e = a * np.exp(b * (x - x0))
    return e / (1 + e)


def life_expectancy_from_mu(mu_func, x_from=0.0, x_max=400.0, n=400001):
    """Numerical e(x_from) = integral of S(t)/S(x_from) for a hazard function, trapezoid on fine grid."""
    t = np.linspace(x_from, x_max, n)
    m = mu_func(t)
    H = np.concatenate([[0], np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(t))])
    S = np.exp(-H)
    return float(np.trapezoid(S, t)), t, S


def lifetable_e0(qx, start=0):
    """Life expectancy at `start` from single-year qx (uniform deaths within year), closing the table
    at the last age with e = 1/mu for the open interval."""
    q = np.asarray(qx[start:], float)
    l = np.concatenate([[1.0], np.cumprod(1 - q)])
    L = 0.5 * (l[:-1] + l[1:])
    # open-ended tail after last age: remaining survivors l[-1] live 1/mu_last more years
    mu_last = -np.log1p(-q[-1])
    return float(L.sum() + l[-1] / mu_last)


def keyfitz_entropy(qx):
    """Keyfitz life-table entropy H = -int l ln l / int l, from single-year qx (piecewise-exponential)."""
    q = np.asarray(qx, float)
    mu = -np.log1p(-q)
    # fine grid within each year
    sub = 200
    t = np.arange(len(q) * sub) / sub
    mu_t = np.repeat(mu, sub)
    H = np.concatenate([[0], np.cumsum(mu_t / sub)])[:-1]
    l = np.exp(-H)
    return float(np.sum(l * H) / np.sum(l))
