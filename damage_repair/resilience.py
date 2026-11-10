"""
Whole-organism dynamics.
 (a) Loss of resilience: Ornstein-Uhlenbeck state with a recovery rate that declines linearly with age
     (Pyrkov et al. 2021 structure).  Estimation of the critical age; hazard from barrier crossing.
 (b) Frailty network: simplified version of the Farrell/Taneja/Mitnitski/Rockwood/Rutenberg generic network
     model.  Scale-free network, damage/repair rates exponential in local damaged fraction.  ILLUSTRATIVE.
 (c) Reliability theory (Gavrilov & Gavrilova 2001) with and without repair: exact Markov chain.
"""
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import expm
import sympy as sp

from common import FIG, save_results, save_csv, gompertz_slope

T0 = time.time()


# ======================================================= (a) critical slowing down
def recovery_rate_line(tau40_wk, tau90_wk):
    """Linear decline of recovery rate r = 1/tau through ages 40 and 90; returns (r0 per week at age 0, t_c)."""
    r40, r90 = 1 / tau40_wk, 1 / tau90_wk
    slope = (r90 - r40) / 50.0
    t_c = 40 - r40 / slope
    return r40 - slope * 40, t_c


def ou_mfpt_hazard(eps_per_day, D, Xc, n=2001):
    """Exact quasi-static hazard (per year) for dx = -eps x dt + sqrt(2D) dW, absorbing at |x| = Xc.
    Symmetric: reflecting at 0, absorbing at Xc.  T = (1/D) int_0^Xc e^{U(y)/D} int_0^y e^{-U(z)/D} dz dy."""
    eps = np.atleast_1d(eps_per_day)[:, None]
    x = np.linspace(0, Xc, n)[None, :]
    u = 0.5 * eps * x ** 2 / D
    um = u.max(axis=1, keepdims=True)
    inner = cumulative_trapezoid(np.exp(-(u - um)), x[0], axis=1, initial=0)
    T = np.trapezoid(np.exp(u - um) * inner, x[0], axis=1) / D
    return 365.0 / T


def resilience():
    out = {}
    for label, (a, b) in {"results_2to8wk": (2.0, 8.0), "discussion_3to6wk": (3.0, 6.0)}.items():
        r0, tc = recovery_rate_line(a, b)
        out[label] = dict(tau40_wk=a, tau90_wk=b, extrapolated_critical_age=tc)
    # simulation of the estimation problem (synthetic individuals, weekly samples, 1 year per person)
    rng = np.random.default_rng(5)
    r0, tc = recovery_rate_line(3.0, 6.0)
    ages = rng.uniform(40, 90, 3000)
    rate_wk = r0 * (1 - ages / tc)
    n_obs = 52
    a1 = np.exp(-rate_wk)                     # lag-1 (1 week) autocorrelation
    x = np.zeros((ages.size, n_obs))
    x[:, 0] = rng.normal(0, 1, ages.size)     # stationary variance normalised to 1 at each age? no: B/eps
    var_stat = 1.0 / rate_wk                  # sigma^2 = B/eps with B = 1 (ILLUSTRATIVE units)
    x[:, 0] *= np.sqrt(var_stat)
    for j in range(1, n_obs):
        x[:, j] = a1 * x[:, j - 1] + np.sqrt(var_stat * (1 - a1 ** 2)) * rng.normal(size=ages.size)
    xc = x - x.mean(1, keepdims=True)
    ac = (xc[:, 1:] * xc[:, :-1]).sum(1) / (xc ** 2).sum(1)
    bins = np.arange(40, 91, 5)
    mids, rate_hat, var_hat = [], [], []
    for lo in bins[:-1]:
        m = (ages >= lo) & (ages < lo + 5)
        mids.append(lo + 2.5)
        rate_hat.append(-np.log(np.clip(ac[m].mean(), 1e-3, 0.999)))
        var_hat.append(x[m].var())
    mids, rate_hat, var_hat = map(np.array, (mids, rate_hat, var_hat))
    p = np.polyfit(mids, rate_hat, 1)
    pv = np.polyfit(mids, 1 / var_hat, 1)
    out["simulated_estimation"] = dict(
        true_t_c=tc, n_people=3000, weeks_each=n_obs,
        t_c_from_recovery_rate=-p[1] / p[0], t_c_from_inverse_variance=-pv[1] / pv[0],
        note="Synthetic data; short series bias the autocorrelation estimate (known small-sample bias).")
    # hazard from barrier crossing with declining recovery rate
    #   ln h ~ -eps(t) Xc^2/(2D) => d ln h/dt = eps0 Xc^2/(2 D t_c)  (Gompertz)
    t = sp.symbols("t", positive=True)
    e0, Xc_, D_, tc_ = sp.symbols("epsilon_0 X_c D t_c", positive=True)
    lnh = -e0 * (1 - t / tc_) * Xc_ ** 2 / (2 * D_)
    out["gompertz_from_resilience_symbolic"] = str(sp.diff(lnh, t))
    # Barrier chosen NUMERICALLY so the exact first-passage hazard has Gompertz b = 0.085/yr over ages 50-90
    # (rate quoted by Pyrkov et al. 2021).  The leading-order formula overstates the slope because the
    # prefactor (~ eps^{3/2}) falls with age.
    from scipy.optimize import brentq, least_squares
    eps0_day = r0 / 7.0
    D = 1.0
    agegrid = np.linspace(0, tc + 5, 800)

    def haz(barrier0, tcx, e0d=eps0_day, Dx=D, ages=agegrid):
        Xc_ = np.sqrt(2 * Dx * barrier0 / e0d)
        ee = np.maximum(e0d * (1 - ages / tcx), 0.0)
        return ou_mfpt_hazard(ee, Dx, Xc_), Xc_
    b_of = lambda B0: np.log(2) / gompertz_slope(agegrid, haz(B0, tc)[0], 50, 90)[1]
    barrier0 = brentq(lambda B0: b_of(B0) - 0.085, 3.0, 60.0)
    h, Xc = haz(barrier0, tc)
    out["hazard_model"] = dict(t_c=tc, barrier_at_birth_in_noise_units=barrier0,
                               leading_order_barrier_would_be=0.085 * tc,
                               gompertz_doubling_50_90=gompertz_slope(agegrid, h, 50, 90)[1],
                               hazard_at={str(a): float(np.interp(a, agegrid, h)) for a in (40, 80, 100, 120, tc)},
                               approx_human_target_at_40=(1 / 3000) * 2 ** (15 / 8))
    # Joint fit of (barrier, t_c) to the approximate human target hazard (1/3000 at 25, x2 per 8 yr), ages 40-90
    tg_a = np.arange(40, 91, 5.0)
    tg_h = (1 / 3000) * 2 ** ((tg_a - 25) / 8.0)
    def resid(th):
        B0, tcx = np.exp(th)
        return np.log(haz(B0, tcx, ages=tg_a)[0]) - np.log(tg_h)
    fj = least_squares(resid, np.log([barrier0 * 1.5, 160.0]), bounds=(np.log([1, 95]), np.log([500, 2000])))
    out["joint_fit_to_human_target"] = dict(barrier0=float(np.exp(fj.x[0])), t_c=float(np.exp(fj.x[1])),
                                            rms_log_error=float(np.sqrt(np.mean(fj.fun ** 2))),
                                            note="recovery rate at age 0 fixed from the 3-wk/6-wk line; target APPROXIMATE")
    # interventions: (i) recovery rate x1.2 at all ages, (ii) noise D x0.8, (iii) slower decline (t_c +20)
    ivs = {"baseline": (eps0_day, tc, D), "recovery_x1.2": (1.2 * eps0_day, tc, D),
           "noise_x0.8": (eps0_day, tc, 0.8 * D), "decline_slower_tc+20": (eps0_day, tc + 20, D)}
    curves = {"age": agegrid}
    for name, (e_, tc2, D2) in ivs.items():
        ee = np.maximum(e_ * (1 - agegrid / tc2), 0)
        hh = ou_mfpt_hazard(ee, D2, Xc)       # same threshold Xc in all scenarios
        H = cumulative_trapezoid(hh, agegrid, initial=0)
        S = np.exp(-H)
        out.setdefault("interventions", {})[name] = dict(
            e0=float(np.trapezoid(S, agegrid)), doubling_50_90=gompertz_slope(agegrid, hh, 50, 90)[1],
            age_hazard_reaches_1=float(agegrid[np.argmax(hh > 1)]) if (hh > 1).any() else None)
        curves[name] = hh
    return out, (mids, rate_hat, p), curves


# ======================================================= (b) frailty network
def scale_free(N, m, rng):
    """Barabasi-Albert preferential attachment."""
    A = np.zeros((N, N), np.float32)
    targets = list(range(m))
    repeated = []
    for v in range(m, N):
        for t in set(targets):
            A[v, t] = A[t, v] = 1
        repeated += list(set(targets)) + [v] * m
        targets = list(rng.choice(repeated, m))
    return A


def network_model(A, n_ind, G0, gp, gm, R, dt, max_age, rng, scale_plus=1.0, scale_minus=1.0,
                  pulse=None):
    """Damage rate G+ = G0 exp(gp f), repair rate G- = (G0/R) exp(-gm f), f = damaged fraction of neighbours.
    Death when both of the two most-connected nodes are damaged. pulse = (start, period, fraction repaired)."""
    N = A.shape[0]
    k = A.sum(1)
    order = np.argsort(-k)
    mort = order[:2]
    fi_nodes = order[2:32]
    Ainv = A / k[None, :]                      # column-normalised so (d @ A/k)_i = sum_j d_j A_ji / k_i
    d = np.zeros((n_ind, N), np.float32)
    alive = np.ones(n_ind, bool)
    death = np.full(n_ind, np.inf)
    steps = int(max_age / dt)
    fi_t = np.zeros(steps)
    next_pulse = None if pulse is None else pulse[0]
    for s in range(steps):
        age = (s + 1) * dt
        da = d[alive]
        f = da @ Ainv
        pp = 1 - np.exp(-scale_plus * G0 * np.exp(gp * f) * dt)
        pm = 1 - np.exp(-scale_minus * (G0 / R) * np.exp(-gm * f) * dt)
        u = rng.random(da.shape, dtype=np.float32)
        da = np.where(da == 0, (u < pp).astype(np.float32), (u >= pm).astype(np.float32))
        if next_pulse is not None and age >= next_pulse:
            fix = rng.random(da.shape) < pulse[2]
            da[fix] = 0
            next_pulse += pulse[1]
        d[alive] = da
        fi_t[s] = da[:, fi_nodes].mean() if da.shape[0] else np.nan
        dead_now = (da[:, mort[0]] == 1) & (da[:, mort[1]] == 1)
        idx = np.where(alive)[0][dead_now]
        death[idx] = age
        alive[idx] = False
        if not alive.any():
            break
    return death, np.arange(1, steps + 1) * dt, fi_t


def lt_hazard(death, max_age, w=2.0):
    edges = np.arange(0, max_age + w, w)
    dd, _ = np.histogram(death[np.isfinite(death)], edges)
    at = len(death) - np.concatenate([[0], np.cumsum(dd)[:-1]])
    py = (at - dd) * w + dd * w / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        return edges[:-1] + w / 2, np.where(py > 0, dd / py, np.nan), dd


def frailty():
    rng = np.random.default_rng(11)
    A = scale_free(200, 5, rng)
    # ILLUSTRATIVE parameters; G0 tuned by hand for a human-like median death age
    G0, gp, gm, R = 0.0016, 7.5, 6.5, 1.5
    out = dict(params=dict(N=200, m=5, mean_degree=float(A.sum() / 200), G0=G0, gamma_plus=gp, gamma_minus=gm,
                           R=R, note="ILLUSTRATIVE; simplified from Farrell et al. (arXiv 1601.07063) which used N=800,"
                                     " Kramers-type rates and Gillespie updates"))
    curves, fis = {}, {}
    scen = {"baseline": {}, "damage_rate_x0.8": dict(scale_plus=0.8), "repair_rate_x1.5": dict(scale_minus=1.5),
            "repair_rate_x10": dict(scale_minus=10.0), "repair_rate_x100": dict(scale_minus=100.0),
            "repair_50pct_of_damage_every_5y_from_60": dict(pulse=(60.0, 5.0, 0.5))}
    for name, kw in scen.items():
        death, tg, fi = network_model(A, 2000, G0, gp, gm, R, 0.1, 150.0, np.random.default_rng(3), **kw)
        mids, h, dd = lt_hazard(death, 150.0)
        ok = dd >= 15
        dfin = np.minimum(death, 150.0)
        lo, hi = np.percentile(dfin, [10, 90])
        out[name] = dict(restricted_mean_150=float(dfin.mean()), median=float(np.median(dfin)),
                         frac_alive_150=float(np.mean(~np.isfinite(death))),
                         doubling_10th_to_90th_pct_death_ages=gompertz_slope(mids[ok], h[ok], lo, min(hi, 149))[1],
                         FI_at={str(a): float(np.interp(a, tg, fi)) for a in (40, 60, 80, 100)})
        curves[name] = (mids, h, ok)
        fis[name] = (tg, fi)
        print(f"[{time.time()-T0:.0f}s] network {name}")
    for name in scen:
        out[name]["delta_mean"] = out[name]["restricted_mean_150"] - out["baseline"]["restricted_mean_150"]
    return out, curves, fis


# ======================================================= (c) reliability with repair
def reliability(m=1, n=10, k=0.01, q=0.5, repair_rates=(0.0, 0.001, 0.005, 0.02)):
    """Each of m blocks has n parallel elements, failure rate k each, initially working with prob q
    (binomial). Failed elements repaired at rate rho each. Block fails when all elements failed (absorbing).
    System hazard = m * block hazard (blocks independent). ILLUSTRATIVE parameters."""
    from scipy.stats import binom
    ages = np.linspace(0.0, 150.0, 601)
    out, curves = {}, {"age": ages}
    for rho in repair_rates:
        Q = np.zeros((n + 1, n + 1))          # state = number working; 0 absorbing
        for j in range(1, n + 1):
            Q[j, j - 1] = j * k
            if j < n:
                Q[j, j + 1] = (n - j) * rho
            Q[j, j] = -Q[j].sum()
        p0 = binom.pmf(np.arange(n + 1), n, q)
        p0[0] = 0
        p0 /= p0.sum()                        # condition on block alive at birth
        P = np.array([p0 @ expm(Q * t) for t in ages])
        S_block = 1 - P[:, 0]
        h_block = (P[:, 1] * k) / S_block
        h = m * h_block
        curves[f"rho={rho}"] = h
        out[f"rho={rho}"] = dict(hazard_at={str(a): float(np.interp(a, ages, h)) for a in (10, 40, 80, 120, 150)},
                                 doubling_20_60=gompertz_slope(ages, h, 20, 60)[1],
                                 doubling_60_100=gompertz_slope(ages, h, 60, 100)[1],
                                 late_limit_no_repair=m * k if rho == 0 else None)
    # Gavrilov closed form check (no repair, early ages): R e^{alpha x}, alpha = kq(n-1)/(1-q)
    out["gavrilov_alpha"] = k * q * (n - 1) / (1 - q)
    out["gavrilov_doubling"] = np.log(2) / out["gavrilov_alpha"]
    out["params"] = dict(m=m, n=n, k=k, q=q, note="ILLUSTRATIVE")
    return out, curves


def main():
    res = {}
    res["resilience"], (mids, rate_hat, p), ou_curves = resilience()
    print(f"[{time.time()-T0:.0f}s] resilience done")
    res["frailty_network"], net_curves, fis = frailty()
    res["reliability"], rel_curves = reliability()
    res["runtime_s"] = time.time() - T0
    save_results("resilience", res)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    ax[0, 0].plot(mids, rate_hat, "o", label="estimated from synthetic weekly series")
    xx = np.linspace(40, 160, 50)
    ax[0, 0].plot(xx, np.polyval(p, xx), "--", label=f"linear extrapolation, zero at {-p[1]/p[0]:.0f}")
    ax[0, 0].axhline(0, c="k", lw=0.5)
    ax[0, 0].set(xlabel="Age", ylabel="Recovery rate (per week)", title="(a) Critical slowing down (SYNTHETIC)")
    ax[0, 0].legend(fontsize=7)
    for name in ["baseline", "recovery_x1.2", "noise_x0.8", "decline_slower_tc+20"]:
        ax[0, 1].semilogy(ou_curves["age"], ou_curves[name], label=name)
    ax[0, 1].set(xlim=(20, 170), ylim=(1e-5, 50), xlabel="Age", ylabel="Hazard (per yr)",
                 title="(a) Hazard from loss of resilience (ILLUSTRATIVE)")
    ax[0, 1].legend(fontsize=7)
    for name, (m_, h_, ok) in net_curves.items():
        ax[1, 0].semilogy(m_[ok], h_[ok], "o-", ms=3, label=name)
    ax[1, 0].set(xlabel="Age", ylabel="Hazard (per yr)", title="(b) Frailty network model (ILLUSTRATIVE)")
    ax[1, 0].legend(fontsize=6.5)
    for name in ["rho=0.0", "rho=0.001", "rho=0.005", "rho=0.02"]:
        ax[1, 1].semilogy(rel_curves["age"], rel_curves[name], label=f"repair {name}")
    ax[1, 1].set(xlabel="Age", ylabel="Hazard (per yr)", title="(c) Reliability model with repair (ILLUSTRATIVE)")
    ax[1, 1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/resilience.png", dpi=140); plt.close(fig)
    save_csv("resilience_ou_hazard", pd.DataFrame(ou_curves), "ILLUSTRATIVE OU-barrier hazard; recovery-time line from "
             "Pyrkov et al. 2021 Discussion values (3 wk at 40, 6 wk at 90).")
    rows = []
    for name, (m_, h_, ok) in net_curves.items():
        for a, hv, o in zip(m_, h_, ok):
            rows.append(dict(scenario=name, age=a, hazard=hv, enough_deaths=o))
    save_csv("resilience_network_hazard", pd.DataFrame(rows), "ILLUSTRATIVE simplified frailty network model.")
    fi_df = pd.DataFrame({"age": fis["baseline"][0], **{k: v[1] for k, v in fis.items()}}).iloc[::10]
    save_csv("resilience_network_FI", fi_df, "Mean frailty index (30 most-connected non-mortality nodes), ILLUSTRATIVE.")
    save_csv("resilience_reliability", pd.DataFrame(rel_curves), "ILLUSTRATIVE reliability model with repair.")
    print(f"resilience done in {time.time()-T0:.0f}s")


if __name__ == "__main__":
    main()
