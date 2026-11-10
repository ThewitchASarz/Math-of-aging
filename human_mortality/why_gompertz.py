"""Part 3: why Gompertz? Derive and simulate
 (1) Gavrilov-Gavrilova reliability model (series of redundant blocks, Poisson initial functional elements),
 (2) Vaupel gamma-frailty model (individual Gompertz -> population logistic with plateau b/s2),
 (3) Strehler-Mildvan vitality model (linear vitality loss + exponentially distributed challenges -> Gompertz;
     implies ln a = ln K - b/B, the compensation-type correlation).
All parameter values in the simulations are ILLUSTRATIVE (synthetic), chosen to make the shape visible,
except the frailty example which uses the female gamma-Gompertz fit from results.json['fits'].
"""
import json
import numpy as np
import sympy as sp
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import update_results, RESULTS, FIG, DATA, mu_ggm

rng = np.random.default_rng(42)
out = {}


def empirical_hazard(T, edges):
    """Life-table hazard estimate from simulated lifetimes: deaths / person-years in each age bin."""
    T = np.sort(T)
    mu = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        alive = T[T > lo]
        if len(alive) < 50:
            mu.append(np.nan); continue
        deaths = np.sum(alive <= hi)
        py = np.sum(np.minimum(alive, hi) - lo)
        mu.append(deaths / py)
    return np.array(mu)


# ---------------- (1) reliability model ----------------
x, k, lam, m = sp.symbols("x k lambda m", positive=True)
i = sp.symbols("i", integer=True, nonnegative=True)
# block with i functional elements (each fails at rate k, independent): S_i = 1 - (1-e^{-kx})^i
# Poisson(lam) mixture: S_block = sum_i e^-lam lam^i/i! [1 - (1-e^{-kx})^i] = 1 - exp(-lam e^{-kx})
S_block_series = sp.summation(sp.exp(-lam) * lam**i / sp.factorial(i) * (1 - (1 - sp.exp(-k * x))**i), (i, 0, sp.oo))
S_block_closed = 1 - sp.exp(-lam * sp.exp(-k * x))
chk = sp.simplify(S_block_series - S_block_closed)
mu_org = sp.simplify(-m * sp.diff(sp.log(S_block_closed), x))
out["reliability_symbolic"] = {"S_block_series_minus_closed_simplified": str(chk),
                               "organism_hazard": str(mu_org),
                               "small_x_limit": "m*lam*k*exp(-lam)*exp(lam*k*x)  (Gompertz: R=m lam k e^-lam, slope lam k)",
                               "large_x_limit": str(sp.limit(mu_org, x, sp.oo))}

K_, LAM, M = 0.02, 5.0, 20  # ILLUSTRATIVE
N = 200_000
n_el = rng.poisson(LAM, size=(N, M))
# block lifetime = max of n_el exponentials (rate k); 0 elements -> dead at birth (time 0)
# max of n iid Exp(k): inverse CDF  -ln(1 - U^{1/n})/k
U = rng.random((N, M))
with np.errstate(divide="ignore", invalid="ignore"):
    tb = np.where(n_el > 0, -np.log1p(-U ** (1.0 / np.maximum(n_el, 1))) / K_, 0.0)
T = tb.min(axis=1)
T = T[T > 0]  # condition on being alive at birth (all blocks have >=1 element)
edges = np.arange(0, 301, 5.0)
mid = 0.5 * (edges[1:] + edges[:-1])
emp = empirical_hazard(T, edges)
xx = mid
y = LAM * np.exp(-K_ * xx)
theory = M * LAM * K_ * np.exp(-K_ * xx) * np.exp(-y) / (-np.expm1(-y))
gomp_approx = M * LAM * K_ * np.exp(-LAM) * np.exp(LAM * K_ * xx)
ok = ~np.isnan(emp)
rel_err = np.abs(emp[ok] / theory[ok] - 1)
out["reliability_simulation"] = {
    "params_ILLUSTRATIVE": {"k_element_failure_rate": K_, "lambda_mean_functional_elements": LAM,
                            "m_blocks": M, "N_individuals": N, "seed": 42},
    "n_alive_at_birth": int(len(T)),
    "gompertz_slope_theory_lam*k": LAM * K_, "plateau_theory_m*k": M * K_,
    "median_abs_rel_error_sim_vs_exact_hazard": float(np.median(rel_err)),
    "max_abs_rel_error_sim_vs_exact_hazard_bins_with_>=50_alive": float(np.max(rel_err)),
    "fitted_loglinear_slope_ages_0_40_sim": float(np.polyfit(xx[(xx < 40) & ok], np.log(emp[(xx < 40) & ok]), 1)[0]),
    "exact_mean_loglinear_slope_ages_0_40": float((np.log(M*LAM*K_*np.exp(-K_*40)/np.expm1(LAM*np.exp(-K_*40)))
                                                    - np.log(M*LAM*K_/np.expm1(LAM)))/40),
    "local_log_slope_formula": "d ln mu/dx = k*( y/(1-exp(-y)) - 1 ),  y = lam*exp(-k x); ~ k(lam e^{-kx} - 1) for large y: "
                               "the slope is largest at birth (~ lam k - k) and declines steadily with age",
    "sim_hazard_last_bins": [float(v) for v in emp[ok][-3:]],
}
pd.DataFrame({"age_mid": xx, "sim_hazard": emp, "exact_hazard": theory, "gompertz_approx": gomp_approx}) \
    .to_csv(f"{DATA}/reliability_model.csv", index=False)

# ---------------- (2) gamma frailty ----------------
a_, b_, z_, s2_ = sp.symbols("a b z sigma2", positive=True)
# individual: mu(x|z) = z a e^{bx}; H0 = (a/b)(e^{bx}-1); z ~ Gamma(shape 1/s2, scale s2)
H0 = a_ / b_ * (sp.exp(b_ * x) - 1)
shape = 1 / s2_
dens = z_**(shape - 1) * sp.exp(-z_ / s2_) / (sp.gamma(shape) * s2_**shape)
S_pop = sp.integrate(dens * sp.exp(-z_ * H0), (z_, 0, sp.oo), conds="none")
S_pop_closed = (1 + s2_ * H0)**(-1 / s2_)
mu_pop = sp.simplify(-sp.diff(sp.log(S_pop_closed), x))
out["frailty_symbolic"] = {
    "S_pop_from_integral_equals_closed": bool(sp.simplify(sp.log(S_pop) - sp.log(S_pop_closed)) == 0)
        if not isinstance(S_pop, sp.Piecewise) else "piecewise (check numerically)",
    "population_hazard": str(mu_pop),
    "limit_x_to_inf": str(sp.limit(mu_pop, x, sp.oo)),
    "mean_frailty_of_survivors": "zbar(x) = 1/(1 + s2 H0(x))",
}
fits = json.load(open(RESULTS))["fits"]["models"]["female"]["gamma_gompertz"]["params"]
A, B, S2 = fits["a"], fits["b"], fits["s2"]
Nf = 1_000_000
z = rng.gamma(1 / S2, S2, size=Nf)
E = rng.exponential(size=Nf)
Tf = np.log1p(E * B / (A * z)) / B  # solve z (a/b)(e^{bT}-1) = E
edges_f = np.arange(40, 131, 2.0)
mid_f = 0.5 * (edges_f[1:] + edges_f[:-1])
emp_f = empirical_hazard(Tf, edges_f)
th_f = A * np.exp(B * mid_f) / (1 + S2 * A / B * np.expm1(B * mid_f))
okf = ~np.isnan(emp_f)
# mean frailty among survivors at 100 vs theory
H100 = A / B * np.expm1(B * 100)
out["frailty_simulation"] = {
    "params": {"a": A, "b": B, "s2": S2, "source": "female gamma-Gompertz fit to SSA 2023 ages 70-110", "N": Nf, "seed": 42},
    "plateau_theory_b/s2": B / S2,
    "median_abs_rel_error_sim_vs_theory": float(np.median(np.abs(emp_f[okf] / th_f[okf] - 1))),
    "mean_frailty_survivors_at_100_sim": float(z[Tf > 100].mean()),
    "mean_frailty_survivors_at_100_theory": float(1 / (1 + S2 * H100)),
    "individual_hazard_doubling_time_ln2/b": float(np.log(2) / B),
}
pd.DataFrame({"age_mid": mid_f, "sim_hazard": emp_f, "theory_pop_hazard": th_f,
              "individual_gompertz_z1": A * np.exp(B * mid_f)}).to_csv(f"{DATA}/frailty_model.csv", index=False)

# ---------------- (3) Strehler-Mildvan ----------------
# vitality V(x) = V0 (1 - B x); challenges arrive at rate K, magnitude ~ Exp(mean D); death if magnitude > V(x)
# => mu(x) = K exp(-V(x)/D) = K e^{-V0/D} exp((V0 B / D) x): Gompertz a = K e^{-V0/D}, b = V0 B/D
# => ln a = ln K - b / B  (negative linear relation between ln a and b; all lines meet at x = 1/B, mu = K)
Ksm, V0, Bsm, D = 2.0, 10.0, 1 / 110.0, 1.0  # ILLUSTRATIVE
Nsm = 200_000
# simulate: thinning — each challenge kills with prob exp(-V(x)/D); do via time-stepping
dt = 0.05
alive = np.ones(Nsm, bool); Tsm = np.full(Nsm, np.nan)
for step in range(int(110 / dt)):
    t0 = step * dt
    V = V0 * (1 - Bsm * t0)
    p_kill = 1 - np.exp(-Ksm * np.exp(-V / D) * dt)
    dead = alive & (rng.random(Nsm) < p_kill)
    Tsm[dead] = t0 + dt / 2
    alive &= ~dead
Tsm = Tsm[~np.isnan(Tsm)]
edges_s = np.arange(0, 111, 5.0); mid_s = 0.5 * (edges_s[1:] + edges_s[:-1])
emp_s = empirical_hazard(np.concatenate([Tsm, np.full(alive.sum(), 1e9)]), edges_s)
th_s = Ksm * np.exp(-V0 * (1 - Bsm * mid_s) / D)
oks = ~np.isnan(emp_s)
slope_sim = np.polyfit(mid_s[oks], np.log(emp_s[oks]), 1)[0]
out["strehler_mildvan"] = {
    "params_ILLUSTRATIVE": {"K": Ksm, "V0": V0, "B": Bsm, "D": D},
    "gompertz_b_theory_V0B/D": V0 * Bsm / D, "gompertz_b_sim": float(slope_sim),
    "gompertz_a_theory_K e^-V0/D": Ksm * np.exp(-V0 / D),
    "median_abs_rel_error_sim_vs_theory": float(np.median(np.abs(emp_s[oks] / th_s[oks] - 1))),
    "compensation_relation": "ln a = ln K - b/B ; populations differing only in V0/D share the convergence point x*=1/B where mu=K",
}
# real-data illustration: male vs female Gompertz lines cross at:
fm = json.load(open(RESULTS))["fits"]
out["strehler_mildvan"]["US2023_male_female_gompertz_crossing_age"] = fm["sex_convergence_age_gompertz"]
update_results("why_gompertz", out)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
ax = axes[0]
ax.semilogy(xx, emp, "o", ms=4, label="simulation (200k organisms)")
ax.semilogy(xx, theory, "-", label="exact hazard (derived)")
ax.semilogy(xx[xx < 45], gomp_approx[xx < 45], "--", label="Gompertz approx. R e^{λk x} (valid kx<<1)")
ax.set_xlim(0, 200)
ax.axhline(M * K_, color="k", ls=":", label="plateau m k")
ax.set_title("Reliability model (ILLUSTRATIVE params:\nm=20 blocks, Poisson λ=5 elements, k=0.02/yr)")
ax.set_xlabel("age"); ax.set_ylabel("hazard (per year)"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
ax = axes[1]
ax.semilogy(mid_f, emp_f, "o", ms=4, label="simulation (1M individuals)")
ax.semilogy(mid_f, th_f, "-", label="population hazard a e^{bx}/(1+σ²H₀)")
ax.semilogy(mid_f, A * np.exp(B * mid_f), "--", label="individual Gompertz (z=1)")
ax.axhline(B / S2, color="k", ls=":", label="plateau b/σ²")
ax.set_title("Gamma-frailty model\n(params = female gamma-Gompertz fit to SSA 2023, ages 70-110)")
ax.set_xlabel("age"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
ax = axes[2]
ax.semilogy(mid_s, emp_s, "o", ms=4, label="simulation (200k)")
ax.semilogy(mid_s, th_s, "-", label="K exp(-V(x)/D)")
ax.set_title("Strehler-Mildvan vitality model\n(ILLUSTRATIVE params)")
ax.set_xlabel("age"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
plt.tight_layout(); plt.savefig(f"{FIG}/why_gompertz_models.png", dpi=150)
pd.DataFrame({"age_mid": mid_s, "sim_hazard": emp_s, "theory": th_s}).to_csv(f"{DATA}/strehler_mildvan_model.csv", index=False)

if __name__ == "__main__":
    print(json.dumps(out, indent=1, default=float))
