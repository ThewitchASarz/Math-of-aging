"""Shared parameters and helpers for the escape_velocity stream.

Rewritten in the verified run (the earlier interrupted version was discarded).
Every empirical number below carries the URL of a source opened in this session.
Values tagged PROVISIONAL, ILLUSTRATIVE or JUDGEMENT are not measurements.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

# Categorical palette (dataviz default instance), fixed order, never cycled.
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "font.size": 10, "axes.titlesize": 11, "axes.titlecolor": INK, "legend.frameon": False,
})

SRC = {
    "cdc2024": "https://cdc.gov/nchs/products/databriefs/db548.htm",
    "hus2016": "https://www.cdc.gov/nchs/data/hus/2016/015.pdf",
    "oeppen_vaupel": "https://u.demog.berkeley.edu/~jrw/Biblio/Eprints/%20M-O/oeppen.vaupel.2002.pdf",
    "japan2024": "https://www.mhlw.go.jp/english/database/db-hw/lifetb24/dl/lifetb24-01.pdf",
    "japan2024_hist": "https://www.mhlw.go.jp/english/database/db-hw/lifetb24/dl/lifetb24-02.pdf",
    "japan2017": "https://www.mhlw.go.jp/english/database/db-hw/lifetb17/dl/lifetb17-06.pdf",
    "olshansky_scimex": "https://www.scimex.org/newsfeed/our-lifespans-arent-getting-much-longer",
    "olshansky_news": "https://www.news-medical.net/news/20241008/Study-suggests-humans-are-approaching-a-biological-limit-to-life.aspx",
    "ho_hendi_bmj": "https://www.bmj.com/lookup/doi/10.1136/bmj.k2562",
    "scott2021": "https://www.london.edu/faculty-and-research/academic-research/a/alls-well-that-ages-well-the-economic-value-of-targeting-aging",
    "mcmillan2024": "https://phys.org/news/2024-02-small-nosed-dogs-longest.amp",
}

# ------------------------------------------------------------ human data
# US 2024 age-specific death rates per 100,000 (CDC NCHS data brief 548).
US2024_RATES = {  # (lower, upper) : rate
    (25, 35): 124.5, (35, 45): 213.9, (45, 55): 386.9,
    (55, 65): 859.8, (65, 75): 1772.9, (75, 85): 4264.0,
}
US2024_RATE_5_14 = 14.5e-5      # all causes, ages 5-14 (lowest age group)
US2024_RATE_15_24 = 66.9e-5     # all causes, ages 15-24
US2024_INJURY_AADR = 53.3e-5    # unintentional injury, age-adjusted
B_HUMAN = (0.085, 0.09)    # 0.085 = fitted both-sexes Gompertz slope (0.0851, US 2023, ages 40-94); 0.09 lies inside the fitted range 0.083-0.101
B_DOG = (0.30, 0.32, 0.34)   # FITTED dog Gompertz slopes: 0.305 (US), 0.320 (UK), 0.344 (US, with background term); see dog_lifespan/results.json
DOG_MEDIAN = 12.5   # McMillan 2024: medians 12.7 (pure breeds), 12.0 (crossbreeds)
C_DOG = 0.005       # FITTED dog background hazard (Gompertz-Makeham c = 0.0050 in both the UK and US tables)


def gm_survival(x, a, b, c):
    """Gompertz-Makeham survival from age 0 to age x: hazard mu(x)=a e^{bx}+c."""
    return np.exp(-c * x - (a / b) * np.expm1(b * x))


def gm_hazard(x, a, b, c):
    return a * np.exp(b * x) + c


def period_e(x, a, b, c):
    """Period remaining life expectancy at age x (closed integrand, numeric quad)."""
    z = a * np.exp(b * x)  # senescent hazard at age x
    f = lambda s: np.exp(-c * s - (z / b) * np.expm1(b * s))
    upper = 50.0 / b + 1.0
    # find s where integrand < 1e-12 to bound the integral
    while f(upper) > 1e-12 and upper < 1e5:
        upper *= 2
    val, _ = quad(f, 0, upper, limit=500)
    return val


def fit_us2024():
    """Fit Gompertz-Makeham to US 2024 group rates (midpoint approximation)."""
    mids = np.array([(lo + hi) / 2 for lo, hi in US2024_RATES])
    m = np.array(list(US2024_RATES.values())) / 1e5

    def res(p):
        la, b, lc = p
        return np.log(np.exp(la) * np.exp(b * mids) + np.exp(lc)) - np.log(m)

    sol = least_squares(res, x0=[np.log(2e-5), 0.085, np.log(5e-4)])
    la, b, lc = sol.x
    return dict(a=float(np.exp(la)), b=float(b), c=float(np.exp(lc)),
                rmse_log=float(np.sqrt(np.mean(sol.fun ** 2))))


def human_params(b, c):
    """Recalibrate a so that mu(70) equals the US 2024 65-74 rate."""
    target = US2024_RATES[(65, 75)] / 1e5
    return dict(a=(target - c) / np.exp(b * 70.0), b=b, c=c)


def dog_params(b, c=C_DOG, median=DOG_MEDIAN):
    a = brentq(lambda a: gm_survival(median, a, b, c) - 0.5, 1e-9, 1.0)
    return dict(a=a, b=b, c=c)


def save_json(name, obj):
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(obj, fh, indent=2, default=float)


def finish(fig, name, note=None):
    if note:
        fig.text(0.01, 0.005, note, fontsize=7, color=MUTED, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.03 if note else 0, 1, 1))
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=150)
    plt.close(fig)


def cohort_e(x0, p, rho):
    """Expected remaining life for someone aged x0 if a(t) falls at rate rho forever.

    Hazard along the person's own life: mu(s) = z0*exp((b-rho)*s) + c, z0 = a e^{b x0}.
    """
    z0 = p["a"] * np.exp(p["b"] * x0)
    k = p["b"] - rho

    def H(s):
        g = z0 * (np.expm1(k * s) / k if abs(k) > 1e-12 else s)
        return p["c"] * s + g
    hi = 1.0
    while H(hi) < 60 and hi < 1e6:
        hi *= 2
    hi = brentq(lambda s: H(s) - 60, 0, hi) if H(hi) >= 60 else hi
    pts = np.linspace(0, hi, 20001)
    return float(np.trapezoid(np.exp(-np.array([H(u) for u in pts])), pts))
