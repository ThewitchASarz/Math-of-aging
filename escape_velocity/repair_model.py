"""Part 3: periodic damage-repair model.

K damage types D_k(t) (measured in units of log-hazard). Between therapies
    dD_k/dt = r_k + g_k * D_k        (g_k = 0: linear accumulation; g_k > 0: self-amplifying)
Every T years a therapy multiplies each D_k by (1 - f_k).
Hazard: mu(t) = c + a * exp(sum_k D_k(t)).
With g_k = 0 and no therapy, sum_k D_k(age) = b * age, which reproduces the Gompertz law.
All damage-type splits, f, T and g values are ILLUSTRATIVE inputs.
"""
import numpy as np
import pandas as pd
import sympy as sp

from common import C, DATA, INK2, MUTED, US2024_INJURY_AADR, finish, human_params, plt, save_json

np.random.seed(0)
out = {}
P = human_params(0.085, US2024_INJURY_AADR)
A, B, CM = P["a"], P["b"], P["c"]
K = 7
START_AGE = 50.0

# ---------------------------------------------------------- closed forms (sympy)
D, r, g, f, T, n = sp.symbols("D r g f T n", positive=True)
# linear: fixed point of D -> (1-f)(D + rT)
fp_lin = sp.solve(sp.Eq(D, (1 - f) * (D + r * T)), D)[0]
# exponential: D -> (1-f)(D e^{gT} + r/g (e^{gT}-1)); contraction factor (1-f)e^{gT}
fp_exp = sp.solve(sp.Eq(D, (1 - f) * (D * sp.exp(g * T) + r / g * (sp.exp(g * T) - 1))), D)[0]
out["fixed_point_linear_post_repair"] = str(sp.simplify(fp_lin))
out["fixed_point_linear_pre_repair"] = str(sp.simplify(fp_lin + r * T))
out["fixed_point_exponential_post_repair"] = str(sp.simplify(fp_exp))
out["bounded_condition_exponential"] = "(1-f)*exp(g*T) < 1  <=>  f > 1 - exp(-g*T)"


def iterate_map(f_, T_, r_, g_, D0=0.0, n_=4000):
    Dn = D0
    for _ in range(n_):
        grow = np.exp(g_ * T_) if g_ > 0 else 1.0
        add = (r_ / g_ * (grow - 1)) if g_ > 0 else r_ * T_
        Dn = (1 - f_) * (Dn * grow + add)
        if Dn > 1e6:
            return np.inf
    return Dn

# numeric check of both closed forms
chk = []
for f_, T_, r_, g_ in [(0.3, 2, 0.0121, 0), (0.5, 5, 0.0121, 0.05), (0.2, 5, 0.0121, 0.05), (0.2, 4, 0.0121, 0.05)]:
    num = iterate_map(f_, T_, r_, g_)
    if g_ == 0:
        cf = (1 - f_) * r_ * T_ / f_
    else:
        cf = float(fp_exp.subs({f: f_, T: T_, r: r_, g: g_})) if (1 - f_) * np.exp(g_ * T_) < 1 else np.inf
    chk.append(dict(f=f_, T=T_, r=r_, g=g_, iterated=num, closed_form=cf,
                    contraction=(1 - f_) * (np.exp(g_ * T_) if g_ else 1)))
out["closed_form_checks"] = chk


# ---------------------------------------------------------- simulation
def simulate(f_vec, T_, r_vec, g_vec, start_age=START_AGE, max_years=30000, pts=40):
    """Remaining life expectancy from start_age; therapy every T_ years from start_age.

    Returns remaining LE, and arrays (age, hazard) for the first 150 years.
    """
    f_vec, r_vec, g_vec = map(np.asarray, (f_vec, r_vec, g_vec))
    # damage at start age with no prior therapy
    Dk = np.where(g_vec > 0, r_vec / np.where(g_vec > 0, g_vec, 1) * np.expm1(g_vec * start_age), r_vec * start_age)
    Hcum, LE, t0 = 0.0, 0.0, 0.0
    ages, haz = [], []
    tau = np.linspace(0, T_, pts)
    while t0 < max_years:
        Dk = (1 - f_vec) * Dk                      # therapy at start of cycle
        grow = np.where(g_vec > 0, np.exp(np.outer(tau, g_vec)), 1.0)
        add = np.where(g_vec > 0, (r_vec / np.where(g_vec > 0, g_vec, 1)) * (grow - 1), np.outer(tau, r_vec))
        Dt = Dk * grow + add                       # (pts, K)
        logz = np.log(A) + Dt.sum(1)
        mu = CM + np.exp(np.minimum(logz, 50))
        Hc = Hcum + np.concatenate([[0], np.cumsum(0.5 * (mu[1:] + mu[:-1]) * np.diff(tau))])
        S = np.exp(-Hc)
        LE += np.trapezoid(S, tau)
        if t0 < 150:
            ages.extend(start_age + t0 + tau); haz.extend(mu)
        Hcum, Dk, t0 = Hc[-1], Dt[-1], t0 + T_
        if Hcum > 40:
            break
    return LE, np.array(ages), np.array(haz)


r_eq = np.full(K, B / K)       # ILLUSTRATIVE: equal shares of the Gompertz slope
g0 = np.zeros(K)
scen = {
    "no therapy": (np.zeros(K), 1.0, r_eq, g0),
    "all 7 types, f=1.0 every year": (np.ones(K), 1.0, r_eq, g0),
    "all 7 types, f=0.5 every 5 y": (np.full(K, 0.5), 5.0, r_eq, g0),
    "all 7 types, f=0.2 every 10 y": (np.full(K, 0.2), 10.0, r_eq, g0),
    "6 of 7 f=1.0 yearly, 1 type f=0": (np.r_[np.ones(6), 0.0], 1.0, r_eq, g0),
    "6 of 7 f=1.0 yearly, 1 type f=0.5": (np.r_[np.ones(6), 0.5], 1.0, r_eq, g0),
}
res, traj = [], {}
for lab, (fv, T_, rv, gv) in scen.items():
    le, ag, hz = simulate(fv, T_, rv, gv)
    res.append(dict(scenario=lab, remaining_LE_from_50=le, expected_age_at_death=START_AGE + le))
    traj[lab] = (ag, hz)
res = pd.DataFrame(res)
out["scenarios_age50"] = res.to_dict(orient="records")

# analytic check for 'no therapy' (pure Gompertz-Makeham from 50)
from common import period_e
out["no_therapy_LE_analytic"] = period_e(50, A, B, CM)

# unrepaired share sweep: one type carries share s of the slope and gets f=0; rest f=1 yearly
shares = [0.0, 0.02, 0.05, 1 / 7, 0.25, 1 / 3, 0.5, 1.0]
sw = []
for s_ in shares:
    rv = np.r_[np.full(6, B * (1 - s_) / 6), B * s_]
    fv = np.r_[np.ones(6), 0.0]
    le, _, _ = simulate(fv, 1.0, rv, g0)
    sw.append(dict(unrepaired_share=s_, residual_slope=B * s_,
                   residual_doubling_time_yr=(np.log(2) / (B * s_)) if s_ > 0 else np.inf,
                   remaining_LE_from_50=le))
sw = pd.DataFrame(sw)
sw.to_csv(f"{DATA}/fig_unrepaired_share.csv", index=False)
out["unrepaired_share_sweep"] = sw.to_dict(orient="records")

# self-amplifying case: g=0.05 for all types (ILLUSTRATIVE), r scaled so that
# untreated log-hazard rise from 30 to 80 still equals b*50 on average
g_amp = 0.03
from scipy.optimize import brentq
# calibrate r so that the untreated remaining LE at 50 matches the Gompertz baseline (33.3 y)
target = period_e(50, A, B, CM)
r_amp = brentq(lambda rr: simulate(np.zeros(K), 1.0, np.full(K, rr), np.full(K, g_amp))[0] - target,
               1e-4, B / K, xtol=1e-6)
out["self_amplifying_r_per_type"] = r_amp
grid = []
for T_ in (1, 2, 5, 10):
    for f_ in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8):
        thr = 1 - np.exp(-g_amp * T_)
        rv = np.full(K, r_amp)  # calibrated, g ILLUSTRATIVE
        le, _, _ = simulate(np.full(K, f_), T_, rv, np.full(K, g_amp), max_years=20000)
        grid.append(dict(T=T_, f=f_, threshold_f=thr, bounded=f_ > thr, remaining_LE_from_50=le))
le_unt, _, _ = simulate(np.zeros(K), 1.0, np.full(K, r_amp), np.full(K, g_amp))
out["self_amplifying_untreated_LE_from_50"] = le_unt
grid = pd.DataFrame(grid)
grid.to_csv(f"{DATA}/repair_threshold_grid.csv", index=False)
out["self_amplifying_grid_g0.03"] = grid.to_dict(orient="records")

# ---------------------------------------------------------- figures
fig, ax = plt.subplots(figsize=(8, 4.2))
show = ["no therapy", "all 7 types, f=0.5 every 5 y", "6 of 7 f=1.0 yearly, 1 type f=0", "all 7 types, f=1.0 every year"]
rows = []
for i, lab in enumerate(show):
    ag, hz = traj[lab]
    m = ag <= 150
    ax.plot(ag[m], hz[m], color=C[i], lw=1.5, label=lab)
    rows.append(pd.DataFrame({"scenario": lab, "age": ag[m], "hazard": hz[m]}))
pd.concat(rows).to_csv(f"{DATA}/fig_repair_hazard.csv", index=False)
ax.set_yscale("log"); ax.set_ylim(1e-4, 5)
ax.set_xlabel("age (years); therapy starts at 50"); ax.set_ylabel("annual death rate (log scale)")
ax.set_title("Periodic repair: hazard stays bounded only if every damage type is repaired")
ax.legend(fontsize=8, loc="upper left")
finish(fig, "fig_repair_hazard", "ILLUSTRATIVE: 7 equal damage types summing to Gompertz slope 0.085; level calibrated to US 2024 rate at 70.")

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(sw.unrepaired_share, sw.remaining_LE_from_50 + 50, color=C[0], marker="o", ms=5)
for _, rw in sw.iterrows():
    ax.annotate(f"{rw.remaining_LE_from_50 + 50:.0f}", (rw.unrepaired_share, rw.remaining_LE_from_50 + 50),
                textcoords="offset points", xytext=(5, 4), fontsize=8, color=INK2)
ax.set_yscale("log")
ax.set_xlabel("share of the aging slope carried by the one unrepairable damage type")
ax.set_ylabel("expected age at death (log scale)")
ax.set_title("The least-repaired link sets the lifespan")
finish(fig, "fig_unrepaired_share", "ILLUSTRATIVE model; other damage fully removed every year from age 50; background hazard 53.3/100k/yr.")

fig, ax = plt.subplots(figsize=(6.5, 4))
Ts = np.linspace(0.1, 20, 200)
for i, gg in enumerate((0.01, 0.03, 0.1)):
    ax.plot(Ts, 1 - np.exp(-gg * Ts), color=C[i], label=f"self-amplification rate g = {gg}/yr")
pd.DataFrame({"T": Ts, **{f"fmin_g{gg}": 1 - np.exp(-gg * Ts) for gg in (0.01, 0.03, 0.1)}}).to_csv(
    f"{DATA}/fig_repair_threshold.csv", index=False)
ax.set_xlabel("years between treatments, T"); ax.set_ylabel("minimum fraction repaired, f")
ax.set_title("Damage stays bounded above each curve: f > 1 - exp(-gT)")
ax.legend(fontsize=8)
finish(fig, "fig_repair_threshold", "Proven threshold for self-amplifying damage; g values ILLUSTRATIVE.")

save_json("results_repair.json", out)
print(pd.DataFrame(chk).to_string())
print(res.round(1).to_string()); print("analytic no-therapy", out["no_therapy_LE_analytic"])
print(sw.round(3).to_string())
print(grid.round(3).to_string()); print("untreated amplifying", le_unt)
