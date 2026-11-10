"""Shared helpers for the damage_repair stream: results.json merging, CSV saving, figure style."""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DAT = os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)
os.makedirs(DAT, exist_ok=True)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def save_results(section, payload):
    """Merge `payload` into results.json under key `section`."""
    path = os.path.join(HERE, "results.json")
    data = {}
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
    data[section] = _clean(payload)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def save_csv(name, df, header_comment):
    """Save a figure's plotted series to data/<name>.csv with a comment header."""
    path = os.path.join(DAT, name + ".csv")
    with open(path, "w") as f:
        for line in header_comment.strip().splitlines():
            f.write("# " + line + "\n")
        df.to_csv(f, index=False)


def gompertz_slope(ages, hazard, lo, hi):
    """Least-squares slope of log hazard on age over [lo, hi]; returns (b, doubling time)."""
    ages = np.asarray(ages)
    hazard = np.asarray(hazard)
    m = (ages >= lo) & (ages <= hi) & (hazard > 0) & np.isfinite(hazard)
    if m.sum() < 3:
        return np.nan, np.nan
    b = np.polyfit(ages[m], np.log(hazard[m]), 1)[0]
    return b, np.log(2) / b
