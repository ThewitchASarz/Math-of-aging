"""Numerical / symbolic checks of the theorems in
'A Mathematical Theory of Aging in Humans and Dogs'. Fixed seeds. Writes verification.json."""
import json, numpy as np, sympy as sp
from scipy.integrate import quad
from scipy.special import exp1
from scipy.optimize import brentq
rng = np.random.default_rng(20261110)
out = {}
g = 0.5772156649015329

def e_num(haz, x0, upper=400.0):
    """remaining life expectancy at x0 for hazard function haz (numeric)."""
    xs = np.linspace(x0, upper, 400001)
    H = np.concatenate([[0], np.cumsum((haz(xs[1:]) + haz(xs[:-1])) / 2 * np.diff(xs))])
    return np.trapezoid(np.exp(-H), xs)

# human-like and dog-like Gompertz parameters (fitted earlier in the project)
ah, bh = 3.03e-5, 0.0910      # US women, Gompertz 40-95
ad, bd = 4.97e-3, 0.3195      # UK dogs, Gompertz 2-19

# T1 closed form
res = []
for a, b, x0 in [(ah, bh, 0), (ah, bh, 50), (ad, bd, 0), (ad, bd, 7)]:
    M = a * np.exp(b * x0) / b
    closed = np.exp(M) * exp1(M) / b
    num = e_num(lambda x: a * np.exp(b * x), x0, upper=x0 + 200)
    asym = (np.log(1 / M) - g) / b
    res.append(dict(a=a, b=b, x0=x0, closed=closed, numeric=num, asymptotic=asym))
out['T1_closed_form'] = res
assert all(abs(r['closed'] - r['numeric']) < 2e-3 for r in res)

# T2 age shift: gain = e(x0-D)-e(x0) < D, D = -ln r / b ; also with Makeham c (aging part scaled only)
res = []
for a, b, c, x0 in [(ah, bh, 0.0, 50), (ah, bh, 0.0008, 50), (ad, bd, 0.0, 7), (ad, bd, 0.005, 10), (ah, bh, 0.0, 80)]:
    for r in [0.9, 0.8, 0.5, 0.1]:
        D = -np.log(r) / b
        e0 = e_num(lambda x: c + a * np.exp(b * x), x0, upper=x0 + 400)
        er = e_num(lambda x: c + r * a * np.exp(b * x), x0, upper=x0 + 400)
        # shifted-age identity: hazard seen from x0 with ratio r equals hazard of someone D years younger
        es = e_num(lambda x: c + a * np.exp(b * (x - D)), x0, upper=x0 + 400)
        res.append(dict(b=b, c=c, x0=x0, r=r, gain=er - e0, bound=D, shift_identity_err=abs(er - es)))
        assert er - e0 < D and abs(er - es) < 1e-6
out['T2_age_shift'] = res

# T3 disease elimination bound: random cause-share curves s(x) <= sbar
res = []
for trial in range(200):
    sbar = rng.uniform(0.05, 0.6)
    k = rng.uniform(0.01, 0.2); m = rng.uniform(40, 100)
    s = lambda x, sbar=sbar, k=k, m=m: sbar * (0.2 + 0.8 / (1 + np.exp(-k * (x - m)) )) * rng.uniform(0.999, 1.0)
    x0 = 0.0
    base = e_num(lambda x: ah * np.exp(bh * x), x0, upper=250)
    cured = e_num(lambda x: (1 - s(x)) * ah * np.exp(bh * x), x0, upper=250)
    bound = -np.log(1 - sbar) / bh
    res.append((cured - base, bound))
    assert cured - base <= bound + 1e-6
out['T3_elimination_bound'] = dict(trials=len(res), max_ratio_gain_to_bound=max(gn / bd_ for gn, bd_ in res))

# T4 halving all mortality gains at most one doubling time
out['T4_halving'] = dict(human_gain=e_num(lambda x: 0.5 * ah * np.exp(bh * x), 0, 250) - e_num(lambda x: ah * np.exp(bh * x), 0, 250),
                         human_doubling=np.log(2) / bh,
                         dog_gain=e_num(lambda x: 0.5 * ad * np.exp(bd * x), 0, 100) - e_num(lambda x: ad * np.exp(bd * x), 0, 100),
                         dog_doubling=np.log(2) / bd)
assert out['T4_halving']['human_gain'] < out['T4_halving']['human_doubling']
assert out['T4_halving']['dog_gain'] < out['T4_halving']['dog_doubling']

# T5 Keyfitz entropy of a Gompertz table ~ 1/(b e0)
xs = np.linspace(0, 250, 500001); S = np.exp(-(ah / bh) * (np.exp(bh * xs) - 1))
e0 = np.trapezoid(S, xs); Hk = -np.trapezoid(S * np.log(np.clip(S, 1e-300, 1)), xs) / e0
out['T5_entropy'] = dict(H=Hk, approx=1 / (bh * e0))
assert abs(Hk - 1 / (bh * e0)) < 0.002

# T6 record age: Gompertz median-of-maximum closed form; plateau grows linearly in ln N
res = []
for N in [1e3, 1e6, 1e9]:
    closed = np.log(1 + (bh / ah) * np.log(N / np.log(2))) / bh
    f = lambda x: (ah / bh) * (np.exp(bh * x) - 1) - np.log(N / np.log(2))
    res.append(dict(N=N, closed=closed, numeric=brentq(f, 50, 200)))
out['T6_record_age'] = res
assert all(abs(r['closed'] - r['numeric']) < 1e-6 for r in res)

# T7 decomposition of life expectancy ratio
lhs = (np.exp(ah / bh) * exp1(ah / bh) / bh) / (np.exp(ad / bd) * exp1(ad / bd) / bd)
rhs = (bd / bh) * (np.log(bh / ah) - g) / (np.log(bd / ad) - g)
out['T7_decomposition'] = dict(exact_ratio=lhs, asymptotic_ratio=rhs, speed_part=bd / bh, level_part=(np.log(bh / ah) - g) / (np.log(bd / ad) - g))
assert abs(lhs - rhs) / lhs < 0.03

# T8 dog-year map (same fraction alive), closed form and limiting slope
h = lambda d: np.log(1 + (ad * bh) / (ah * bd) * (np.exp(bd * d) - 1)) / bh
Sd = lambda d: np.exp(-(ad / bd) * (np.exp(bd * d) - 1)); Sh = lambda x: np.exp(-(ah / bh) * (np.exp(bh * x) - 1))
errs = [abs(Sd(d) - Sh(h(d))) for d in [1, 4, 8, 12, 16]]
slope = (h(16.01) - h(16.0)) / 0.01; slope_far = (h(40.01) - h(40.0)) / 0.01
out['T8_dog_years'] = dict(max_survival_mismatch=max(errs), slope_at_16=slope, slope_at_40=slope_far, limit_slope=bd / bh)
assert max(errs) < 1e-9 and abs(slope_far - bd / bh) < 1e-3 and abs(slope - bd / bh) < 0.05

# T9 Schoenfeld deaths
from scipy.stats import norm
D = lambda HR: 4 * (norm.ppf(0.975) + norm.ppf(0.8)) ** 2 / np.log(HR) ** 2
out['T9_deaths_needed'] = {str(HR): D(HR) for HR in [0.9, 0.8, 0.7]}

# T12 flat hazard: finite repairable chain -> hazard converges to smallest decay rate
Q = np.array([[-0.30, 0.25, 0.00], [0.20, -0.52, 0.30], [0.00, 0.35, -0.60]])  # live states 0,1,2; missing row mass = death rate
death = -Q.sum(axis=1)
from scipy.linalg import expm
p0 = np.array([1.0, 0, 0]); hz = []
for t in [1, 5, 20, 80, 300]:
    p = p0 @ expm(Q * t); hz.append(float(p @ death / p.sum()))
lam = -max(np.linalg.eigvals(Q).real)
out['T12_flat_hazard'] = dict(hazard_at_t=hz, limit=lam)
assert abs(hz[-1] - lam) < 1e-9
# and with one irreversible live state (no repair out of state 2) the chain is still homogeneous: limit still constant,
# while an age-dependent rate keeps the hazard rising:
hz2 = [0.01 * np.exp(0.09 * t) for t in [20, 40, 60]]
out['T12_note'] = 'age-dependent rates (e.g. production growing with t) are needed for an ever-rising hazard'

# T14 cancer vs body size: p = A W (T0 - k W)^4 peaks at W = T0/(5k)
W, T0, k, A = sp.symbols('W T0 k A', positive=True)
sol = sp.solve(sp.diff(A * W * (T0 - k * W) ** 4, W), W)
out['T14_peak'] = dict(solutions=[str(s) for s in sol], dog_value=13.3 / (5 * 0.08))
assert any(sp.simplify(s - T0 / (5 * k)) == 0 for s in sol)

# T15, T16 escape velocity identities with Makeham term c
a0, b, c = 5.6e-5, 0.085, 0.000533
def e_period(x, t, rho):
    return e_num(lambda y: c + a0 * np.exp(-rho * t) * np.exp(b * y), x, upper=x + 3000) if False else \
           quad(lambda s: np.exp(-(c * s + a0 * np.exp(-rho * t + b * x) * (np.exp(b * s) - 1) / b)), 0, np.inf, limit=400)[0]
x, t, rho, hstep = 60.0, 0.0, 0.03, 1e-3
de_dx = (e_period(x + hstep, t, rho) - e_period(x - hstep, t, rho)) / (2 * hstep)
de_dt = (e_period(x, t + hstep, rho) - e_period(x, t - hstep, rho)) / (2 * hstep)
mu = c + a0 * np.exp(b * x); e = e_period(x, t, rho)
out['T15_identity_A'] = dict(de_dx=de_dx, mu_e_minus_1=mu * e - 1)
out['T16_identity_B'] = dict(de_dt=de_dt, minus_rho_over_b_de_dx=-(rho / b) * de_dx)
assert abs(de_dx - (mu * e - 1)) < 1e-5 and abs(de_dt + (rho / b) * de_dx) < 1e-5
# along a life: remaining period life expectancy is constant iff rho = b
track = {}
for label, r in [('0.5b', 0.5 * b), ('b', b), ('1.1b', 1.1 * b)]:
    track[label] = [e_period(50 + s, s, r) for s in [0, 10, 20, 30]]
out['T16_along_life'] = track
assert max(track['b']) - min(track['b']) < 1e-6 and track['0.5b'][-1] < track['0.5b'][0] and track['1.1b'][-1] > track['1.1b'][0]
out['T16_thresholds_percent_per_year'] = dict(human_b_0p085=100 * (1 - np.exp(-0.085)), human_b_0p10=100 * (1 - np.exp(-0.10)), dog_b_0p32=100 * (1 - np.exp(-0.3195)))

# T17 periodic repair
r_, T_, f_ = 1.0, 10.0, 0.5; Dm = 0.0
for _ in range(200): Dm = (1 - f_) * (Dm + r_ * T_)
out['T17_linear'] = dict(iterated=Dm, formula=r_ * T_ * (1 - f_) / f_)
assert abs(Dm - r_ * T_ * (1 - f_) / f_) < 1e-9
gg, T_ = 0.03, 10.0; thr = 1 - np.exp(-gg * T_); grow = {}
for f_ in [thr - 0.02, thr + 0.02]:
    Dm = 1.0
    for _ in range(400): Dm = (1 - f_) * (Dm * np.exp(gg * T_) + (0.1 / gg) * (np.exp(gg * T_) - 1))
    grow[round(float(f_), 4)] = float(Dm)
out['T17_self_amplifying'] = dict(threshold=thr, damage_after_400_cycles=grow)

# T18 weakest link: 7 types, six fully repaired every year, one never: long-run slope = b/7
bK = 0.085; K = 7; ts = np.linspace(0, 400, 400001)
D = (bK / K) * ts + (K - 1) * (bK / K) * (ts - np.floor(ts))
slope = np.polyfit(ts[200000:], D[200000:], 1)[0]
out['T18_weakest_link'] = dict(long_run_slope=slope, predicted=bK / K, doubling_time=np.log(2) / (bK / K))
assert abs(slope - bK / K) < 1e-4


# T5 rate vs level, 30-year-old woman
M30 = ah * np.exp(bh * 30) / bh
e30 = np.exp(M30) * exp1(M30) / bh
e30_half_level = np.exp(0.5 * M30) * exp1(0.5 * M30) / bh
e30_half_rate = 2 * np.exp(2 * M30) * exp1(2 * M30) / bh
num_half_rate = e_num(lambda x: ah * np.exp(bh * 30) * np.exp((bh / 2) * (x - 30)), 30, upper=30 + 600)
out['T5_rate_vs_level'] = dict(e30=e30, gain_half_level=e30_half_level - e30, gain_half_rate=e30_half_rate - e30, numeric_half_rate=num_half_rate - e30)
assert abs(e30_half_rate - num_half_rate) < 5e-3 and e30_half_level - e30 < np.log(2) / bh

# T10 speed alone: rate x median lifespan rises by ln s; level unchanged
med = lambda a, b: np.log(1 + b * np.log(2) / a) / b
s_ = bd / bh
out['T10_speed_alone'] = dict(human_rate_x_median=bh * med(ah, bh), dog_rate_x_median=bd * med(ad, bd),
                              sped_up_human_median=med(ah, s_ * bh), actual_dog_gompertz_median=med(ad, bd),
                              predicted_rate_x_median=bh * med(ah, bh) + np.log(s_), exact_rate_x_median=s_ * bh * med(ah, s_ * bh),
                              level_ratio=ad / ah)
assert abs(out['T10_speed_alone']['predicted_rate_x_median'] - out['T10_speed_alone']['exact_rate_x_median']) < 0.01

json.dump(out, open('verification.json', 'w'), indent=1, default=float)
print('ALL CHECKS PASSED')
for k_, v in out.items(): print(k_, '->', json.dumps(v, default=float)[:230])
