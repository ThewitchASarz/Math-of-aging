"""
Why aging exists: Hamilton's forces of selection, mutation accumulation (Medawar), antagonistic pleiotropy
(Williams) and the disposable soma (Kirkwood) as an optimisation problem.

All life histories here are ILLUSTRATIVE (synthetic), chosen to look roughly human-like or wolf/dog-like.
They are not fitted to data.  The mathematics (Euler-Lotka, selection gradients) is standard; the gradient
formula follows Moorad, Promislow & Silvertown (2019, TREE, Box 1), opened at
https://www.zoology.ubc.ca/~otto/veg/Readings/Moorad2019.pdf
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq, minimize_scalar
import sympy as sp

from common import FIG, save_results, save_csv

DX = 0.05
X = np.arange(0, 120 + DX, DX)


def trap(y, x=X):
    return np.trapezoid(y, x)


def survival(mu):
    return np.exp(-np.concatenate([[0], np.cumsum(0.5 * (mu[1:] + mu[:-1]) * DX)]))


def euler_lotka_r(l, m):
    """Solve int e^{-r x} l(x) m(x) dx = 1 for r."""
    f = lambda r: trap(np.exp(-r * X) * l * m) - 1.0
    return brentq(f, -1.5, 3.0)


def forces(l, m, r):
    """Hamilton's forces.  Mortality: s_mu(x) = int_x^inf e^{-ry} l m dy / T  (minus sign dropped)
    Fertility: s_m(x) = e^{-rx} l(x) / T.  T = generation time = int y e^{-ry} l m dy."""
    w = np.exp(-r * X) * l * m
    T = trap(X * w)
    tail = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * DX)])
    s_mu = (tail[-1] - tail) / T
    s_m = np.exp(-r * X) * l / T
    return s_mu, s_m, T


def life_history(kind):
    """ILLUSTRATIVE life histories. Returns (mu, m, params)."""
    if kind == "human_like":
        lam, a, b, alpha, omega = 0.01, 1e-4, 0.08, 16.0, 45.0
    else:  # wolf/dog-like
        lam, a, b, alpha, omega = 0.15, 2e-3, 0.23, 2.0, 10.0
    mu = lam + a * np.exp(b * X)
    fert = ((X >= alpha) & (X <= omega)).astype(float)
    l = survival(mu)
    m0 = 1.0 / trap(l * fert)                     # stationary population (r = 0) at baseline
    return mu, m0 * fert, dict(lam=lam, a=a, b=b, alpha=alpha, omega=omega, m0=m0)


def mutation_accumulation(kind, U_T=1e-3, delta=0.05, n_iter=200):
    """Mutation-selection balance for age-specific deleterious alleles.
    An allele adding hazard delta at age x has selection coefficient s = delta * S(x)/T, S(x)=int_x e^{-ry} l m.
    Equilibrium frequency q = u/s (capped at 1).  Summed over loci acting at age x (total rate U per age class):
       extra hazard  dmu(x) = n u delta / (delta S(x)/T) = U T / S(x)   (independent of delta), capped at n*delta.
    Here U_T := U*T lumps the constants (ILLUSTRATIVE).  Iterate because S depends on l, which depends on dmu."""
    mu0, m, p = life_history(kind)
    cap = 2.0                                       # cap on extra hazard (all loci fixed) - ILLUSTRATIVE
    dmu = np.zeros_like(X)
    for _ in range(n_iter):
        l = survival(mu0 + dmu)
        r = euler_lotka_r(l, m)
        w = np.exp(-r * X) * l * m
        tail = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * DX)])
        S = tail[-1] - tail
        new = np.minimum(U_T / np.maximum(S, 1e-300), cap)
        if np.max(np.abs(new - dmu)) < 1e-10:
            break
        dmu = 0.5 * dmu + 0.5 * new
    return mu0, dmu, r


def ap_threshold(kind, phi=0.05, A_grid=None):
    """Antagonistic pleiotropy: allele raises fertility by fraction phi at all fertile ages and adds an extra
    hazard delta from age A onwards.  Find the largest delta that still lets it invade (r_mutant > r_wild)."""
    mu, m, p = life_history(kind)
    l = survival(mu)
    r0 = euler_lotka_r(l, m)
    out = {}
    for A in A_grid:
        def gain(delta):
            l2 = survival(mu + delta * (X >= A))
            return euler_lotka_r(l2, m * (1 + phi)) - r0
        if gain(2.0) > 0:
            out[A] = np.inf                       # invades even with +2/yr hazard (nothing left to lose)
        else:
            out[A] = brentq(gain, 0.0, 2.0)
    return r0, out


# ---------------------------------------------------------------- disposable soma
def soma_fitness(rho, lam, crit, bmax=0.25, a=1e-3, alpha=2.0, m0=1.0, juvenile_lam=None):
    """Fraction rho of resources to somatic maintenance.  Intrinsic hazard a*exp(b x), b = bmax(1-rho).
    Fertility m0(1-rho) from maturity alpha onward (no menopause, to favour long life as much as possible).
    Extrinsic hazard lam at all ages.  crit = 'R0' (lifetime reproductive success) or 'r' (Euler-Lotka)."""
    b = bmax * (1 - rho) ** 3
    mu = lam + (a * np.exp(b * X) if b > 0 else a)
    l = survival(mu)
    m = m0 * (1 - rho) * (X >= alpha)
    if crit == "R0":
        return trap(l * m)
    return euler_lotka_r(l, m)


def soma_optimum(lam, crit, **kw):
    r = minimize_scalar(lambda rho: -soma_fitness(rho, lam, crit, **kw), bounds=(0.0, 0.999),
                        method="bounded", options=dict(xatol=1e-5))
    return r.x


def symbolic_uniform_lambda():
    """PROVEN: with uniform extra hazard d, l -> l e^{-d x}, and r -> r - d solves Euler-Lotka, so
    e^{-r x} l(x) (and hence every Hamilton force) is unchanged."""
    x, r, d = sp.symbols("x r d", real=True)
    l = sp.Function("l")(x)
    lhs = sp.exp(-(r - d) * x) * l * sp.exp(-d * x)
    return sp.simplify(lhs - sp.exp(-r * x) * l) == 0


def main():
    res = {"note": "All life histories ILLUSTRATIVE (synthetic)."}

    # A. Hamilton forces
    fr = {"age": X}
    res["hamilton"] = {}
    for kind in ["human_like", "dog_like"]:
        mu, m, p = life_history(kind)
        l = survival(mu)
        r = euler_lotka_r(l, m)
        s_mu, s_m, T = forces(l, m, r)
        fr[f"{kind}_force_mortality"] = s_mu
        fr[f"{kind}_force_fertility"] = s_m
        res["hamilton"][kind] = dict(params=p, r=r, generation_time=T,
                                     force_mu_at={str(a): float(np.interp(a, X, s_mu))
                                                  for a in [0, p["alpha"], (p["alpha"] + p["omega"]) / 2,
                                                            p["omega"], p["omega"] + 5]},
                                     force_mu_zero_after=p["omega"])
    res["hamilton"]["uniform_extrinsic_invariance_symbolic"] = bool(symbolic_uniform_lambda())

    # B. Mutation accumulation
    mu0, dmu, r = mutation_accumulation("human_like")
    res["mutation_accumulation"] = dict(
        U_T=1e-3, cap=2.0, r=r,
        extra_hazard_at={str(a): float(np.interp(a, X, dmu)) for a in [20, 30, 40, 45, 50, 60]},
        age_extra_hazard_exceeds_0p1=float(X[np.argmax(dmu > 0.1)]),
        formula="dmu(x) = U T / S(x),  S(x) = int_x^inf e^{-ry} l(y) m(y) dy")

    # C. Antagonistic pleiotropy
    A_grid = [10, 20, 30, 40, 45, 50]
    r0, thr = ap_threshold("human_like", phi=0.05, A_grid=A_grid)
    res["antagonistic_pleiotropy"] = dict(phi_fertility_gain=0.05,
                                          max_tolerated_late_hazard_by_onset_age={str(k): v for k, v in thr.items()})

    # D. Disposable soma: optimum maintenance vs extrinsic mortality
    lam_grid = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.4])
    ds = {}
    rows = []
    # For the r criterion, scale fertility so r = 0 at lam = 0.05 (otherwise r >> 0 and growth discounting
    # dominates). The uniform-lam invariance holds for any m0.
    rho05 = soma_optimum(0.05, "R0")
    m0_r = 1.0 / soma_fitness(rho05, 0.05, "R0")
    for crit in ["R0", "r"]:
        for lam in lam_grid:
            rho = soma_optimum(lam, crit, m0=(1.0 if crit == "R0" else m0_r))
            b = 0.25 * (1 - rho) ** 3
            ds[f"{crit}_lam{lam}"] = dict(rho_opt=rho, gompertz_b=b, doubling=np.log(2) / b if b > 0 else np.inf)
            rows.append(dict(criterion=crit, extrinsic_hazard=lam, rho_opt=rho, gompertz_b=b,
                             doubling_time=np.log(2) / b))
    res["disposable_soma"] = dict(
        model="hazard = lam + a e^{b x}, b = bmax(1-rho)^3, fertility m0(1-rho) from age 2; bmax=0.25, a=1e-3",
        results=ds,
        never_full_repair="rho_opt < 1 always: at rho=1 fertility is 0, so fitness is 0 (PROVEN within model).")

    # E. Dog breeds: size as a perturbation of an evolved optimum (SPECULATIVE toy)
    # Nunney 2024 (RSOS) fit across breeds: lifespan T = 13.3 - 0.08 W (W kg).
    # Toy: growth to size W diverts maintenance: b(W) = b_wolf * (1 + k (W - W0)/W0); calibrate k to that slope.
    a_d = 0.0116; b_ref = np.log(2) / 3.0            # dog-like Gompertz (ILLUSTRATIVE; from saturating_removal target)
    W0 = 30.0

    def e0_dog(W, k):
        b = b_ref * (1 + k * (W - W0) / W0)
        return trap(np.exp(-(a_d / b) * (np.exp(b * X) - 1)))
    k_fit = brentq(lambda k: (e0_dog(40, k) - e0_dog(20, k)) / 20 + 0.08, 0.0, 5.0)
    Wg = np.arange(2, 81, 2.0)
    e0W = np.array([e0_dog(W, k_fit) for W in Wg])
    res["dog_size_toy"] = dict(k=k_fit, slope_target_per_kg=-0.08, note="SPECULATIVE toy: size raises Gompertz "
                               "slope b (consistent with Kraus et al. 2013: large dogs age faster). Calibrated, not predicted.",
                               b_at_10kg=b_ref * (1 + k_fit * (10 - W0) / W0),
                               b_at_60kg=b_ref * (1 + k_fit * (60 - W0) / W0),
                               e0_at={str(int(W)): float(e) for W, e in zip(Wg, e0W) if W in (10, 30, 60)})
    save_results("evolution", res)

    # figures
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for kind, c in [("human_like", "C0"), ("dog_like", "C1")]:
        ax[0].plot(X, fr[f"{kind}_force_mortality"] / fr[f"{kind}_force_mortality"][0], c, label=kind)
    ax[0].set(xlim=(0, 70), xlabel="Age", ylabel="Force of selection on mortality (rel. to birth)",
              title="Hamilton's force (ILLUSTRATIVE life histories)")
    ax[0].legend()
    ax[1].semilogy(X, mu0, label="baseline hazard")
    ax[1].semilogy(X, mu0 + dmu, label="+ mutation-selection balance load")
    ax[1].set(xlim=(0, 70), ylim=(1e-3, 5), xlabel="Age", title="Mutation accumulation (ILLUSTRATIVE)")
    ax[1].legend(fontsize=8)
    df = pd.DataFrame(rows)
    for crit, sty in [("R0", "o-"), ("r", "s--")]:
        d = df[df.criterion == crit]
        ax[2].plot(d.extrinsic_hazard, d.doubling_time, sty, label=f"fitness = {crit}")
    ax[2].set(xscale="log", xlabel="Extrinsic hazard (per yr)", ylabel="Evolved mortality doubling time (yr)",
              title="Disposable soma optimum (ILLUSTRATIVE)")
    ax[2].legend()
    fig.tight_layout(); fig.savefig(f"{FIG}/evolution.png", dpi=150); plt.close(fig)
    save_csv("evolution_forces", pd.DataFrame(fr).iloc[::20], "Hamilton forces, ILLUSTRATIVE life histories.")
    save_csv("evolution_mutation_accumulation", pd.DataFrame(dict(age=X, mu0=mu0, extra=dmu)).iloc[::20],
             "Mutation-selection balance hazard, ILLUSTRATIVE parameters (U*T=1e-3, cap 2/yr).")
    save_csv("evolution_disposable_soma", df, "Disposable soma optimum vs extrinsic hazard, ILLUSTRATIVE model.")
    save_csv("evolution_dog_size_toy", pd.DataFrame(dict(weight_kg=Wg, e0=e0W)),
             "SPECULATIVE toy calibrated to Nunney 2024 slope -0.08 yr/kg.")
    print("evolution done")


if __name__ == "__main__":
    main()
