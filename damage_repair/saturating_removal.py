"""
Saturating-removal (SR) model of senescent cells.
Karin, Agrawal, Porat, Krizhanovsky & Alon (2019) Nat Commun 10:5495
  parameters: https://pmc.ncbi.nlm.nih.gov/articles/PMC6889273
  equation and noise convention (sqrt(2 eps)): Alon lecture notes
  https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/lecture10_srmodel_june28_2020.pdf

    dX = [ eta*t - beta*X/(kappa+X) ] dt + sqrt(2*eps) dW,   reflecting at X=0, death when X >= Xc

X  senescent-cell (SnC) abundance, dimensionless.  t  age in years; time step in days.
eta  production slope (day^-1 year^-1): production per day at age t (years) is p = eta*t.
beta maximal removal rate (day^-1); kappa half-saturation; eps noise (day^-1); Xc death threshold.
lam  age-independent extrinsic hazard (year^-1).

Sections:
 1. Exact quasi-static hazard (mean first-passage time, MFPT) and a closed-form Laplace approximation.
 2. Euler-Maruyama SDE check.
 3. Refits: approximate human target (x2 per 8 yr) and approximate dog-like target (x2 per 3 yr, e0 11.2 yr).
 4. Interventions: more removal (beta), less production (eta), periodic removal pulses.
 5. Condition for a non-rising hazard.
Runtime target < 3 min.
"""
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
import sympy as sp

from common import FIG, save_results, save_csv, gompertz_slope

T0 = time.time()
HUMAN_PUB = dict(eta=0.00135, beta=0.15, kappa=0.5, eps=0.142, Xc=17.0, lam=0.4e-3)   # Fig. 4d
MOUSE_PUB = dict(eta=0.084, beta=0.15, kappa=0.5, eps=0.16, Xc=17.0, lam=0.0)         # Fig. 4b


# ------------------------------------------------------------------ 1. quasi-static hazard
def U(x, p, beta, kappa):
    """Potential with U'(x) = -(drift) = beta x/(kappa+x) - p (constant dropped)."""
    return -p * x + beta * (x - kappa * np.log1p(x / kappa))


def mfpt_days(p, beta, kappa, eps, Xc, n=2501):
    """Exact MFPT (days) from X=0 (reflecting) to Xc (absorbing) with production p held fixed.
    T = (1/eps) int_0^Xc dy e^{U(y)/eps} int_0^y e^{-U(z)/eps} dz.  Vectorised over arrays p, beta."""
    p = np.atleast_1d(np.asarray(p, float))[:, None]
    beta = np.broadcast_to(np.atleast_1d(np.asarray(beta, float)), (p.shape[0],))[:, None]
    x = np.linspace(0.0, Xc, n)[None, :]
    u = U(x, p, beta, kappa) / eps
    # e^{u(y)} int_0^y e^{-u(z)} dz = int_0^y e^{u(y)-u(z)} dz ; shift by row max for stability
    um = u.max(axis=1, keepdims=True)
    inner = cumulative_trapezoid(np.exp(-(u - um)), x[0], axis=1, initial=0.0)
    integrand = np.exp(u - um) * inner
    return np.trapezoid(integrand, x[0], axis=1) / eps


def qs_hazard(ages, eta, beta, kappa, eps, Xc, lam=0.0, beta_fn=None, prod_fn=None):
    ages = np.asarray(ages, float)
    p = eta * ages if prod_fn is None else np.array([prod_fn(a) for a in ages])
    b = beta * np.ones_like(ages) if beta_fn is None else np.array([beta_fn(a) for a in ages])
    return 365.0 / mfpt_days(p, b, kappa, eps, Xc) + lam


def laplace_hazard(ages, eta, beta, kappa, eps, Xc):
    """Closed form (per year): h = 365 * U'(Xc) * sqrt(U''(X*)/(2 pi eps)) * exp(-[U(Xc)-U(X*)]/eps),
    X* = kappa p/(beta-p). Valid when the barrier dU >> eps and U'(Xc) > 0 (i.e. p < beta Xc/(kappa+Xc))."""
    p = eta * np.asarray(ages, float)
    p = np.where(p < beta, p, np.nan)          # formula undefined once production exceeds beta
    xs = kappa * p / (beta - p)
    dU = U(Xc, p, beta, kappa) - U(xs, p, beta, kappa)
    U1c = beta * Xc / (kappa + Xc) - p
    U2s = beta * kappa / (kappa + xs) ** 2
    return 365.0 * U1c * np.sqrt(U2s / (2 * np.pi * eps)) * np.exp(-dU / eps)


def laplace_slope(age, eta, beta, kappa, eps, Xc):
    """d ln h/dt (per year) from the closed form:
    eta (Xc - X*)/eps  -  eta/U'(Xc)  -  eta/(beta - p)      [barrier term, boundary-slope term, curvature term]"""
    p = eta * age
    xs = kappa * p / (beta - p)
    return eta * (Xc - xs) / eps - eta / (beta * Xc / (kappa + Xc) - p) - eta / (beta - p)


def life_stats(ages, h):
    H = cumulative_trapezoid(h, ages, initial=0.0)
    S = np.exp(-H)
    e0 = np.trapezoid(S, ages)                      # restricted to the age grid
    med = np.interp(0.5, S[::-1], ages[::-1]) if S[-1] < 0.5 else np.nan
    return S, e0, med


def symbolic_checks():
    """(i) envelope identity d(dU)/dp = -(Xc - X*); (ii) X* solves p = beta X/(kappa+X);
    (iii) slope of ln h from the closed form, differentiated symbolically, equals laplace_slope()."""
    x, p, b, k, Xc, eta, t, e = sp.symbols("x p beta kappa X_c eta t epsilon", positive=True)
    Ux = -p * x + b * (x - k * sp.log(1 + x / k))
    xs = sp.solve(sp.Eq(sp.diff(Ux, x), 0), x)[0]
    dU = Ux.subs(x, Xc) - Ux.subs(x, xs)
    env_ok = sp.simplify(sp.diff(dU, p) + (Xc - xs)) == 0
    U1 = sp.diff(Ux, x).subs(x, Xc)
    U2 = sp.diff(Ux, x, 2).subs(x, xs)
    lnh = sp.log(U1) + sp.Rational(1, 2) * sp.log(U2) - dU / e
    slope = sp.diff(lnh.subs(p, eta * t), t)
    target = eta * (Xc - xs.subs(p, eta * t)) / e - eta / U1.subs(p, eta * t) - eta / (b - eta * t)
    vals = {eta: 0.00135, b: 0.15, k: 0.5, e: 0.142, Xc: 17, t: 60}
    slope_ok = abs(float(slope.subs(vals)) - float(target.subs(vals))) < 1e-10
    return dict(xstar=str(xs), envelope_identity=bool(env_ok), slope_formula_matches_symbolic=bool(slope_ok))


# ------------------------------------------------------------------ 2. SDE
def simulate_sde(eta, beta, kappa, eps, Xc, lam, n, dt=1.0, max_age=150.0, seed=0,
                 pulse=None, beta_fn=None, prod_fn=None):
    """Euler-Maruyama (dt in days) with reflection at 0.  Returns age at death; np.inf = alive at max_age.
    pulse = (start_age, period_days, fraction_removed)."""
    rng = np.random.default_rng(seed)
    X = np.zeros(n)
    idx = np.arange(n)
    death = np.full(n, np.inf)
    sq = np.sqrt(2 * eps * dt)
    p_ext = lam * dt / 365.0
    next_pulse = None if pulse is None else pulse[0] * 365.0
    for k in range(1, int(max_age * 365 / dt) + 1):
        td = k * dt
        a = td / 365.0
        b = beta if beta_fn is None else beta_fn(a)
        prod = eta * a if prod_fn is None else prod_fn(a)
        X += (prod - b * X / (kappa + X)) * dt + sq * rng.standard_normal(X.size)
        np.abs(X, out=X)
        if next_pulse is not None and td >= next_pulse:
            X *= (1.0 - pulse[2])
            next_pulse += pulse[1]
        dead = X >= Xc
        if p_ext > 0:
            dead |= rng.random(X.size) < p_ext
        if dead.any():
            death[idx[dead]] = a
            X, idx = X[~dead], idx[~dead]
            if X.size == 0:
                break
    return death


def sde_stats(d, max_age=150.0):
    """Censoring-safe: restricted mean lifespan to max_age (survivors counted as max_age)."""
    dd = np.minimum(d, max_age)
    return dict(restricted_mean=float(dd.mean()), median=float(np.median(dd)),
                frac_alive_at_max=float(np.mean(~np.isfinite(d))))


def lifetable_hazard(d, bin_w=1.0, max_age=150):
    edges = np.arange(0, max_age + bin_w, bin_w)
    deaths, _ = np.histogram(d[np.isfinite(d)], edges)
    at_risk = len(d) - np.concatenate([[0], np.cumsum(deaths)[:-1]])
    py = (at_risk - deaths) * bin_w + deaths * bin_w / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.where(py > 0, deaths / py, np.nan)
    return edges[:-1] + bin_w / 2, h, deaths


# ------------------------------------------------------------------ main
def main():
    res = {}
    P = HUMAN_PUB
    ages = np.linspace(0.0, 150.0, 1501)
    res["symbolic"] = symbolic_checks()

    # --- 1. quasi-static, published human set
    h_pub = qs_hazard(ages, P["eta"], P["beta"], P["kappa"], P["eps"], P["Xc"], P["lam"])
    h_int = h_pub - P["lam"]
    S, e0, med = life_stats(ages, h_pub)
    h_lap = laplace_hazard(ages, P["eta"], P["beta"], P["kappa"], P["eps"], P["Xc"])
    loc_b = np.gradient(np.log(h_int), ages)
    res["human_published_quasistatic"] = dict(
        params=P, e0=e0, median=med,
        doubling_40_80=gompertz_slope(ages, h_int, 40, 80)[1],
        doubling_80_100=gompertz_slope(ages, h_int, 80, 100)[1],
        doubling_100_110=gompertz_slope(ages, h_int, 100, 110)[1],
        local_gompertz_b={str(a): float(np.interp(a, ages, loc_b)) for a in [40, 60, 80, 100, 110]},
        laplace_slope={str(a): float(laplace_slope(a, P["eta"], P["beta"], P["kappa"], P["eps"], P["Xc"]))
                       for a in [40, 60, 80, 100]},
        leading_order_slope_eta_Xc_over_eps=P["eta"] * P["Xc"] / P["eps"],
        laplace_over_exact_hazard={str(a): float(np.interp(a, ages, h_lap / h_int)) for a in [40, 60, 80, 100]},
        hazard_per_year={str(a): float(np.interp(a, ages, h_pub)) for a in [30, 50, 70, 80, 90, 100, 110, 120]},
        age_production_equals_removal_at_threshold=P["beta"] * P["Xc"] / (P["kappa"] + P["Xc"]) / P["eta"],
        age_production_equals_beta=P["beta"] / P["eta"])

    # --- 2. SDE check
    d_pub = simulate_sde(P["eta"], P["beta"], P["kappa"], P["eps"], P["Xc"], P["lam"], n=20000, seed=1)
    mids, h_sde, deaths = lifetable_hazard(d_pub)
    d_fine = simulate_sde(P["eta"], P["beta"], P["kappa"], P["eps"], P["Xc"], P["lam"], n=4000, dt=0.25,
                          seed=2, max_age=130)
    S_grid_mean = life_stats(ages, h_pub)[1]
    res["human_published_sde"] = dict(n=20000, dt_days=1.0, **sde_stats(d_pub),
                                      doubling_40_80=gompertz_slope(mids, h_sde, 40, 80)[1],
                                      quasistatic_e0=S_grid_mean,
                                      dt_check_0p25day_n4000=sde_stats(d_fine, 130))
    print(f"[{time.time()-T0:.0f}s] SDE baseline done")

    # --- 3a. refit to approximate human target
    tgt_a = np.arange(40, 91, 5.0)
    tgt_h = (1 / 3000) * 2 ** ((tgt_a - 25) / 8.0)

    def rh(th):
        eta, eps, beta = np.exp(th)
        return np.log(qs_hazard(tgt_a, eta, beta, P["kappa"], eps, P["Xc"], P["lam"])) - np.log(tgt_h)

    fh = least_squares(rh, np.log([P["eta"], P["eps"], P["beta"]]))
    eh, sh, bh = np.exp(fh.x)
    h_fit_h = qs_hazard(ages, eh, bh, P["kappa"], sh, P["Xc"], P["lam"])
    _, e0h, medh = life_stats(ages, h_fit_h)
    res["human_refit_8yr"] = dict(
        target="hazard=(1/3000)*2^((age-25)/8), ages 40-90; APPROXIMATE target (Smithsonian summary of Gompertz law)",
        eta=eh, eps=sh, beta=bh, kappa=P["kappa"], Xc=P["Xc"],
        doubling_40_90=gompertz_slope(ages, h_fit_h - P["lam"], 40, 90)[1], e0=e0h, median=medh,
        rms_log_error=float(np.sqrt(np.mean(fh.fun ** 2))))

    # --- 3b. dog-like target: Gompertz, doubling 3 yr (approximate, from project brief, UNVERIFIED here),
    # level set so that e0 = 11.2 yr (VetCompass UK companion dogs, rvc.ac.uk news page).
    b_d = np.log(2) / 3.0
    gg = np.linspace(0, 40, 4001)

    def e0_gomp(a):
        return np.trapezoid(np.exp(-(a / b_d) * (np.exp(b_d * gg) - 1)), gg)
    lo, hi = 1e-5, 1.0
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if e0_gomp(mid) > 11.2 else (lo, mid)
    a_d = np.sqrt(lo * hi)
    tgt_ad = np.arange(4, 17, 1.0)
    tgt_hd = a_d * np.exp(b_d * tgt_ad)

    def rd(th, free):
        v = dict(eta=P["eta"], eps=P["eps"], beta=P["beta"])
        v.update(zip(free, np.exp(th)))
        return np.log(qs_hazard(tgt_ad, v["eta"], v["beta"], P["kappa"], v["eps"], P["Xc"])) - np.log(tgt_hd)

    dog = {}
    aa = ages[ages <= 40]
    for key, free, x0 in [("eta_only", ["eta"], [P["eta"] * 7]),
                          ("eta_eps", ["eta", "eps"], [P["eta"] * 5, P["eps"] * 2]),
                          ("eta_beta", ["eta", "beta"], [P["eta"] * 3, P["beta"] / 2])]:
        f = least_squares(lambda th: rd(th, free), np.log(x0))
        v = dict(eta=P["eta"], eps=P["eps"], beta=P["beta"])
        v.update(zip(free, np.exp(f.x)))
        hh = qs_hazard(aa, v["eta"], v["beta"], P["kappa"], v["eps"], P["Xc"])
        _, e0d, medd = life_stats(aa, hh)
        dog[key] = dict(params=v, eta_fold_vs_human=v["eta"] / P["eta"], beta_fold=v["beta"] / P["beta"],
                        eps_fold=v["eps"] / P["eps"], doubling_4_16=gompertz_slope(aa, hh, 4, 16)[1],
                        e0=e0d, median=medd, rms_log_error=float(np.sqrt(np.mean(f.fun ** 2))))
    # dog-like interventions on both successful fits (model experiments)
    for key in ["eta_eps", "eta_beta"]:
        v = dog[key]["params"]
        base = life_stats(aa, qs_hazard(aa, v["eta"], v["beta"], P["kappa"], v["eps"], P["Xc"]))[1]
        iv = {}
        for nm, bf, ef in [("removal_x1.25", 1.25, 1.0), ("removal_x1.5", 1.5, 1.0),
                           ("production_x0.8", 1.0, 0.8), ("production_x0.5", 1.0, 0.5)]:
            hh = qs_hazard(aa, v["eta"] * ef, v["beta"] * bf, P["kappa"], v["eps"], P["Xc"])
            iv[nm] = dict(e0=life_stats(aa, hh)[1], delta_e0=life_stats(aa, hh)[1] - base,
                          pct_gain=100 * (life_stats(aa, hh)[1] / base - 1),
                          doubling_4_16=gompertz_slope(aa, hh, 4, 16)[1])
        dog[key]["interventions"] = iv
    res["dog_like_fits"] = dict(target=dict(gompertz_a=a_d, gompertz_b=b_d, doubling=3.0, e0=11.2,
                                            note="ILLUSTRATIVE dog-like target: doubling ~3 yr is an approximate "
                                                 "target (unverified here); e0 11.2 yr from VetCompass UK."),
                                fits=dog)

    # time-rescaling invariance: h(t; s*eta) = h(s t; eta) => lifespan x Gompertz-b is fixed if only eta changes
    s = 7.0
    hs = qs_hazard(ages / s, P["eta"] * s, P["beta"], P["kappa"], P["eps"], P["Xc"])
    res["eta_rescaling"] = dict(
        claim="Quasi-static hazard depends on age only through p=eta*t, so scaling eta by s compresses the "
              "whole curve by s: median/s and doubling/s. Ratio median/doubling is invariant.",
        max_abs_log_diff_check=float(np.nanmax(np.abs(np.log(hs[1:]) - np.log(h_int[1:])))),
        human_median_over_doubling=med / gompertz_slope(ages, h_int, 40, 80)[1],
        dog_target_median_over_doubling=float(np.interp(0.5, *(lambda S: (S[::-1], aa[::-1]))(
            np.exp(-(a_d / b_d) * (np.exp(b_d * aa) - 1))))) / 3.0)
    print(f"[{time.time()-T0:.0f}s] fits done")

    # --- 4a. interventions (quasi-static)
    B, E = P["beta"], P["eta"]
    scen = {
        "baseline": {},
        "removal_x1.25_lifelong": dict(beta_fn=lambda a: B * 1.25),
        "removal_x1.5_lifelong": dict(beta_fn=lambda a: B * 1.5),
        "removal_x1.5_from_60": dict(beta_fn=lambda a: B * (1.5 if a >= 60 else 1.0)),
        "production_x0.8_lifelong": dict(prod_fn=lambda a: 0.8 * E * a),
        "production_x0.5_lifelong": dict(prod_fn=lambda a: 0.5 * E * a),
        "production_slope_x0.5_from_60": dict(prod_fn=lambda a: E * (a if a < 60 else 60 + 0.5 * (a - 60))),
        "production_frozen_at_60": dict(prod_fn=lambda a: E * min(a, 60.0)),
        "noise_x0.8_lifelong": "eps",
    }
    out, curves = {}, {"age": ages}
    for name, kw in scen.items():
        if kw == "eps":
            hh = qs_hazard(ages, E, B, P["kappa"], 0.8 * P["eps"], P["Xc"], P["lam"])
        else:
            hh = qs_hazard(ages, E, B, P["kappa"], P["eps"], P["Xc"], P["lam"], **kw)
        _, e0s, meds = life_stats(ages, hh)
        out[name] = dict(e0_restricted_150=e0s, median=meds,
                         doubling_60_90=gompertz_slope(ages, hh - P["lam"], 60, 90)[1],
                         hazard_80=float(np.interp(80, ages, hh)), hazard_100=float(np.interp(100, ages, hh)),
                         hazard_120=float(np.interp(120, ages, hh)))
        curves[name] = hh
    for name in out:
        out[name]["delta_e0"] = out[name]["e0_restricted_150"] - out["baseline"]["e0_restricted_150"]
    res["interventions_quasistatic"] = out

    # --- 4b. SDE: check removal x1.5 and periodic pulses (pulses cannot be done quasi-statically)
    sde = {}
    sde["none"] = sde_stats(d_pub)
    d_b = simulate_sde(E, B, P["kappa"], P["eps"], P["Xc"], P["lam"], n=8000, seed=3, beta_fn=lambda a: B * 1.5)
    sde["removal_x1.5_lifelong"] = sde_stats(d_b)
    for name, pulse in [("remove50pct_yearly_from60", (60, 365.0, 0.5)),
                        ("remove50pct_monthly_from60", (60, 30.0, 0.5)),
                        ("remove50pct_monthly_from75", (75, 30.0, 0.5))]:
        d = simulate_sde(E, B, P["kappa"], P["eps"], P["Xc"], P["lam"], n=8000, seed=4, pulse=pulse)
        m2, h2, dd = lifetable_hazard(d)
        sde[name] = sde_stats(d)
        ok = dd >= 20
        sde[name]["doubling_after_start"] = gompertz_slope(m2[ok], h2[ok], pulse[0] + 5, pulse[0] + 40)[1]
        curves_p = (m2, h2, dd)
        if name == "remove50pct_monthly_from60":
            pulse_curve = curves_p
    for name in sde:
        sde[name]["delta_restricted_mean"] = sde[name]["restricted_mean"] - sde["none"]["restricted_mean"]
    res["interventions_sde"] = sde
    print(f"[{time.time()-T0:.0f}s] SDE interventions done")

    # --- 5. stationary condition: production frozen at p -> constant hazard H(p)
    p_grid = np.linspace(0.0, 0.16, 33)
    Hp = 365.0 / mfpt_days(p_grid, np.full_like(p_grid, B), P["kappa"], P["eps"], P["Xc"])
    res["stationary_condition"] = dict(
        statement="With production held at constant p the hazard is constant H(p): aging stops in this model. "
                  "Hazard stays low only while p is well below the removal ceiling beta (and below "
                  "beta*Xc/(kappa+Xc) at the threshold).",
        H_per_year={f"{p:.3f}": float(h) for p, h in zip(p_grid, Hp)},
        production_at_age={str(a): E * a for a in [20, 40, 60, 80, 100]})
    res["runtime_s"] = time.time() - T0
    save_results("saturating_removal", res)

    # ---------------- figures
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ok = np.isfinite(h_sde) & (deaths >= 20)
    ax[0].semilogy(ages, h_pub, label="Exact quasi-static (published human params)")
    ax[0].semilogy(ages[ages < 105], (h_lap + P["lam"])[ages < 105], ":", label="Closed-form Laplace approx.")
    ax[0].semilogy(mids[ok], h_sde[ok], "o", ms=3, label="SDE simulation (n=20,000)")
    ax[0].semilogy(ages, h_fit_h, "--", label="Refit to approx. x2 per 8 yr")
    ax[0].semilogy(tgt_a, tgt_h, "k+", label="Approx. target 1/3000 at 25, x2 / 8 yr")
    ax[0].set(xlim=(20, 125), ylim=(1e-4, 3), xlabel="Age (years)", ylabel="Hazard (per year)",
              title="Saturating-removal model: human")
    ax[0].legend(fontsize=7)
    for key, sty in [("eta_only", ":"), ("eta_eps", "-"), ("eta_beta", "--")]:
        v = dog[key]["params"]
        ax[1].semilogy(aa[aa <= 25], qs_hazard(aa[aa <= 25], v["eta"], v["beta"], P["kappa"], v["eps"], P["Xc"]),
                       sty, label=f"fit varying {key.replace('_', ' ')}")
    ax[1].semilogy(tgt_ad, tgt_hd, "k+", label="ILLUSTRATIVE dog-like target (x2/3 yr, e0 11.2)")
    ax[1].set(xlabel="Age (years)", title="Dog-like target (illustrative)", ylim=(1e-3, 5))
    ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/saturating_removal_hazard.png", dpi=150); plt.close(fig)
    save_csv("saturating_removal_hazard",
             pd.DataFrame(dict(age=ages, hazard_exact_qs_human_pub=h_pub, hazard_laplace_intrinsic=h_lap,
                               hazard_refit_8yr=h_fit_h)),
             "Karin et al. 2019 SR model, published human parameters (https://pmc.ncbi.nlm.nih.gov/articles/PMC6889273); "
             "refit to an APPROXIMATE x2-per-8-yr target. Hazard per year.")
    save_csv("saturating_removal_hazard_sde", pd.DataFrame(dict(age_mid=mids, hazard_sde=h_sde, deaths=deaths)),
             "SDE simulation, published human parameters, n=20000, dt=1 day, seed 1.")
    dogdf = pd.DataFrame({"age": aa})
    for key in dog:
        v = dog[key]["params"]
        dogdf[f"hazard_fit_{key}"] = qs_hazard(aa, v["eta"], v["beta"], P["kappa"], v["eps"], P["Xc"])
    dogdf["target_ILLUSTRATIVE"] = a_d * np.exp(b_d * aa)
    save_csv("saturating_removal_dog_fits", dogdf,
             "ILLUSTRATIVE dog-like target (Gompertz doubling 3 yr, e0 11.2 yr from VetCompass) and SR-model fits.")

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for name in scen:
        ax.semilogy(ages, curves[name], label=f"{name} (e0 {out[name]['e0_restricted_150']:.0f})")
    okp = pulse_curve[2] >= 20
    ax.semilogy(pulse_curve[0][okp], pulse_curve[1][okp], "k.", ms=4, label="SDE: remove 50% monthly from 60")
    ax.set(xlim=(40, 150), ylim=(1e-3, 5), xlabel="Age (years)", ylabel="Hazard (per year)",
           title="Model 'senolytic' experiments (human params) - model output, not trial data")
    ax.legend(fontsize=6.3); fig.tight_layout()
    fig.savefig(f"{FIG}/saturating_removal_interventions.png", dpi=150); plt.close(fig)
    save_csv("saturating_removal_interventions", pd.DataFrame(curves),
             "Quasi-static hazard per year under model interventions; Karin et al. 2019 human parameters. "
             "Model experiments, not trial results.")
    print(f"saturating_removal done in {time.time()-T0:.0f}s")


if __name__ == "__main__":
    main()
