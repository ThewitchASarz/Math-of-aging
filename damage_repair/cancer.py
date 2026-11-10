"""
Cancer and the multistage model.
 1. Armitage-Doll: k sequential rare hits -> incidence ~ t^(k-1). Exact Markov-chain check.
 2. Non-identifiability: a 2-hit model with clonal expansion of the intermediate stage mimics slope 6.
 3. Peto's paradox: required suppression when cell number N and lifespan T grow; effect of Cagan scaling.
 4. Dogs: Nunney (2024, R Soc Open Sci) fit  p = A4 * W * T^4,  T = 13.3 - 0.08 W.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from math import factorial
import sympy as sp

from common import FIG, save_results, save_csv


def exact_k_stage(k, u, ages):
    """One cell, stages 0..k, rate u for each step. Returns P(cell reached stage k by t)."""
    Q = np.zeros((k + 1, k + 1))
    for i in range(k):
        Q[i, i] = -u
        Q[i, i + 1] = u
    out = []
    for t in ages:
        out.append(expm(Q * t)[0, k])
    return np.array(out)


def main():
    res = {}
    ages = np.linspace(20, 90, 71)

    # 1. Armitage-Doll exact vs power law
    N, u, k = 1e8, 5e-4, 7         # ILLUSTRATIVE: N cells at risk, per-step rate per year
    Pk = exact_k_stage(k, u, ages)
    # tissue hazard = N * dP/dt / (1 - N P)  approx N dP/dt (P tiny)
    dP = np.gradient(Pk, ages)
    h_exact = N * dP
    h_power = N * u ** k * ages ** (k - 1) / factorial(k - 1)
    slope = np.polyfit(np.log(ages), np.log(h_exact), 1)[0]
    t = sp.symbols("t", positive=True)
    uu = sp.symbols("u", positive=True)
    # leading term of P_k(t) for small u: Erlang CDF expansion
    erlang = 1 - sp.exp(-uu * t) * sum((uu * t) ** i / sp.factorial(i) for i in range(k))
    lead = sp.series(erlang, uu, 0, k + 1).removeO()
    res["armitage_doll"] = dict(
        N=N, u=u, k=k, note="ILLUSTRATIVE N,u",
        loglog_slope_exact=float(slope), expected_slope=k - 1,
        max_rel_err_power_law=float(np.max(np.abs(h_power / h_exact - 1))),
        leading_term_symbolic=str(sp.simplify(lead)),
        observed_slope_human_cancers="about 6 (incidence ~ sixth power of age; BJC 2004 editorial on Armitage & Doll 1954)",
        implied_hits_naive=7)

    # 2. Two-stage clonal expansion mimicking slope 6
    g = 0.1   # ILLUSTRATIVE net growth rate of intermediate clones per year
    h_tsce = (np.exp(g * ages) - 1) / g
    loc = g * ages * np.exp(g * ages) / (np.exp(g * ages) - 1)
    res["clonal_expansion"] = dict(
        model="h(t) = N u1 u2 (e^{g t} - 1)/g",
        local_loglog_slope_at={str(a): float(np.interp(a, ages, loc)) for a in (40, 60, 80)},
        fitted_slope_40_80=float(np.polyfit(np.log(ages[(ages >= 40) & (ages <= 80)]),
                                            np.log(h_tsce[(ages >= 40) & (ages <= 80)]), 1)[0]),
        message="2 hits + clonal growth 0.1/yr gives log-log slope ~6 over ages 40-80: slope alone does not identify the number of hits")

    # 3. Peto's paradox: lifetime risk ~ 1 - exp(-N (uT)^k / k!)
    k4 = 4
    ref_risk = 0.3                 # ILLUSTRATIVE reference lifetime risk
    S_grid = np.array([1e-3, 1e-2, 1e-1, 1, 10, 100, 1000])
    tab = []
    for S in S_grid:
        for tau in (0.05, 0.5, 1.0, 2.0):
            x = -np.log(1 - ref_risk) * S * tau ** k4
            risk_naive = 1 - np.exp(-x)
            suppression = (S * tau ** k4) ** (1 / k4)       # fold reduction in EVERY per-step rate to keep risk = ref
            # if per-year mutation rate scales as 1/T (Cagan), u*T is constant -> only S matters
            risk_cagan = 1 - (1 - ref_risk) ** S             # u*T held constant: risk depends on S only
            tab.append(dict(cell_ratio=S, lifespan_ratio=tau, risk_naive=risk_naive,
                            fold_rate_reduction_needed=suppression, risk_if_uT_constant=risk_cagan,
                            fold_reduction_needed_if_uT_constant=S ** (1 / k4)))
    df = pd.DataFrame(tab)
    res["peto"] = dict(k=k4, ref_risk=ref_risk, note="ILLUSTRATIVE reference risk; ratios relative to a human-like reference",
                       example_100x_cells_same_lifespan=df[(df.cell_ratio == 100) & (df.lifespan_ratio == 1)].iloc[0].to_dict(),
                       example_1000x_cells_2x_lifespan=df[(df.cell_ratio == 1000) & (df.lifespan_ratio == 2)].iloc[0].to_dict())

    # 4. Dogs (Nunney 2024)
    A4 = 7.84e-7                       # ALL MEANS fit
    W = np.linspace(2, 80, 157)
    T = 13.3 - 0.08 * W
    p = A4 * W * T ** 4
    Ws = sp.symbols("W", positive=True)
    peak = sp.solve(sp.diff(Ws * (sp.Rational(133, 10) - sp.Rational(8, 100) * Ws) ** 4, Ws), Ws)
    p_fixed_T = A4 * W * 11.0 ** 4
    res["dogs_nunney"] = dict(
        A4=A4, lifespan_fit="T = 13.3 - 0.08 W", peak_weight_kg=[float(x) for x in peak if 0 < float(x) < 166],
        p_at={str(w): float(A4 * w * (13.3 - 0.08 * w) ** 4) for w in (5, 10, 20, 33.25, 50, 70)},
        p_if_lifespan_fixed_11y={str(w): float(A4 * w * 11 ** 4) for w in (5, 33.25, 70)},
        source_reported="breeds 29-42 kg were over-represented among highest-risk; >60 kg breeds lower")
    save_results("cancer", res)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].loglog(ages, h_exact, label="exact 7-stage chain")
    ax[0].loglog(ages, h_power, "--", label="N u^k t^(k-1)/(k-1)!")
    ax[0].loglog(ages, h_tsce * h_exact[40] / np.interp(60, ages, h_tsce), ":", label="2-hit + clonal growth (rescaled)")
    ax[0].set(xlabel="Age (log)", ylabel="Incidence (log, arb.)", title="Multistage incidence (ILLUSTRATIVE)")
    ax[0].legend(fontsize=7)
    for tau in (0.5, 1.0, 2.0):
        d = df[df.lifespan_ratio == tau]
        ax[1].semilogx(d.cell_ratio, d.risk_naive, "o-", label=f"naive, lifespan x{tau}")
    d = df[df.lifespan_ratio == 1.0]
    ax[1].semilogx(d.cell_ratio, d.risk_if_uT_constant, "k--", label="mutation rate ~ 1/lifespan")
    ax[1].set(xlabel="Cells at risk relative to reference", ylabel="Lifetime cancer risk",
              title="Peto's paradox: naive multistage (ILLUSTRATIVE)")
    ax[1].legend(fontsize=7)
    ax[2].plot(W, p, label="Nunney 2024 fit: A4 W T(W)^4")
    ax[2].plot(W, p_fixed_T, ":", label="same, lifespan held at 11 yr")
    ax[2].axvline(33.25, c="k", ls=":")
    ax[2].set(xlabel="Breed weight (kg)", ylabel="Lifetime cancer mortality", title="Dogs: size vs cancer")
    ax[2].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/cancer.png", dpi=150); plt.close(fig)
    save_csv("cancer_multistage", pd.DataFrame(dict(age=ages, h_exact=h_exact, h_power=h_power, h_tsce=h_tsce)),
             "ILLUSTRATIVE multistage incidence curves.")
    save_csv("cancer_peto", df, "ILLUSTRATIVE Peto's paradox table (k=4, reference risk 0.3).")
    save_csv("cancer_dogs", pd.DataFrame(dict(weight_kg=W, lifespan=T, p_cancer=p, p_fixed_lifespan=p_fixed_T)),
             "Nunney 2024 R Soc Open Sci fit (https://escholarship.org/content/qt8w0923ff/qt8w0923ff.pdf)")
    print("cancer done")


if __name__ == "__main__":
    main()
