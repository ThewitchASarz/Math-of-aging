"""
Molecular clocks and counters.
 (a) Telomere shortening / Hayflick limit as a counting process        (simulated)
 (b) Somatic mutation accumulation, Cagan et al. 2022 inverse scaling   (simulated)
 (c) Epigenetic clocks from accumulating stochastic variation, entropy  (simulated, own elastic net)
 (d) mtDNA heteroplasmy as clonal drift (Wright-Fisher)                 (simulated)
 (e) Protein aggregation with clearance (nucleation-elongation-fragmentation)  (solved)
Inputs marked ILLUSTRATIVE are synthetic.  Sources for empirical numbers are in findings.md.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sympy as sp

from common import FIG, save_results, save_csv

RNG = np.random.default_rng(20261110)


# ============================================================ (a) telomeres
def telomeres(n_cells=20000):
    """92 telomere ends per cell. Initial length per end ~ N(L0, sd) (ILLUSTRATIVE L0=10 kb, sd=1.5 kb).
    Each division removes U(50,200) bp per end (range from Harley et al. 1990 via BioNumbers 104276).
    Senescence (Hayflick) when the SHORTEST end falls below Lc (ILLUSTRATIVE 3 kb).
    Counting process: divisions D = min_j (L_j - Lc)/mean loss approx."""
    L0, sd, Lc = 10000.0, 1500.0, 3000.0
    L = RNG.normal(L0, sd, size=(n_cells, 92))
    D = np.zeros(n_cells, int)
    alive = np.ones(n_cells, bool)
    while alive.any():
        L[alive] -= RNG.uniform(50, 200, size=(alive.sum(), 92))
        D[alive] += 1
        alive &= L.min(axis=1) >= Lc
    det_mean_end = (L0 - Lc) / 125.0
    # extreme-value expectation of the shortest of 92 ends
    z_min = np.mean(np.sort(RNG.normal(size=(20000, 92)), axis=1)[:, 0])
    det_shortest = (L0 + z_min * sd - Lc) / 125.0
    # leukocyte arithmetic: years for an adult telomere to lose (L_adult - Lc) at 20-40 bp/yr (secondary source)
    years = {f"{rate}bp_per_yr": (7000 - Lc) / rate for rate in (20, 40)}
    return dict(params=dict(L0=L0, sd=sd, Lc=Lc, loss_bp_per_div="U(50,200)", ends=92,
                            note="L0, sd, Lc ILLUSTRATIVE"),
                hayflick_mean=float(D.mean()), hayflick_sd=float(D.std()), hayflick_cv=float(D.std() / D.mean()),
                deterministic_single_end=det_mean_end, deterministic_shortest_end=det_shortest,
                years_from_7kb_to_3kb_at_population_rate=years), D


# ============================================================ (b) somatic mutations
def somatic_mutations():
    """Cagan et al. 2022: rate per genome per year ~ ELB / L with ELB = 3206.4 (95% CI 2683.9-3728.9);
    human 47/yr, mouse 796/yr (opened).  Dog rate not readable in opened text -> predicted here."""
    ELB, ELB_lo, ELB_hi = 3206.4, 2683.9, 3728.9
    L_h = 83.67                      # human 80th-percentile lifespan used by Cagan (83.7, opened)
    pred_h = ELB / L_h
    dog_L = np.array([12.0, 14.0, 16.0])   # ILLUSTRATIVE range for a dog 80th-percentile lifespan
    pred_dog = {f"L={L:.0f}": dict(rate=ELB / L, lo=ELB_lo / L, hi=ELB_hi / L) for L in dog_L}
    # simulation: gamma-Poisson accumulation in 2000 crypts per species (overdispersion ILLUSTRATIVE, CV 0.2)
    ages_h = np.linspace(0, 100, 101)
    ages_d = np.linspace(0, 20, 101)
    def sim(rate, ages):
        r_i = rate * RNG.gamma(25, 1 / 25, size=2000)            # crypt-to-crypt rate variation, CV 0.2
        return RNG.poisson(np.outer(r_i, np.diff(ages, prepend=0))).cumsum(axis=1)
    mh = sim(47.0, ages_h)
    md = sim(ELB / 14.0, ages_d)
    # multistage consequence: P(crypt has k specific drivers) ~ prod (g_i * n(t)) / k!  -> same at end of life
    return dict(ELB=ELB, human_rate_observed=47.0, human_rate_predicted=pred_h, mouse_rate_observed=796.0,
                dog_rate_predicted=pred_dog, dog_rate_note="UNVERIFIED: measured dog rate not read in this session",
                human_burden_at_80=float(np.interp(80, ages_h, mh.mean(0))),
                dog_burden_at_14_if_L14=float(np.interp(14, ages_d, md.mean(0))),
                ratio_dog_over_human_rate_if_L14=(ELB / 14.0) / 47.0), (ages_h, mh, ages_d, md)


# ============================================================ (c) epigenetic clocks
def elastic_net(Xtr, ytr, lam, alpha=0.5, n_iter=300):
    """Coordinate-descent elastic net on standardised X (glmnet objective)."""
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    Z = (Xtr - mu) / sd
    y0 = ytr.mean()
    r = ytr - y0
    n, p = Z.shape
    w = np.zeros(p)
    zz = (Z ** 2).sum(0) / n
    for _ in range(n_iter):
        w_old = w.copy()
        for j in range(p):
            rj = r + Z[:, j] * w[j]
            rho = Z[:, j] @ rj / n
            w[j] = np.sign(rho) * max(abs(rho) - lam * alpha, 0) / (zz[j] + lam * (1 - alpha))
            r = rj - Z[:, j] * w[j]
        if np.max(np.abs(w - w_old)) < 1e-6:
            break
    return lambda Xn: y0 + ((Xn - mu) / sd) @ w, int((w != 0).sum())


def stochastic_clock(n_feat=500, sigma=0.01, bounded=True, ground=None, seed=1):
    """Meyer & Schumacher 2024 construction: shared random ground state in [0,1]; each year add N(0, sigma^2)
    to every feature independently; clip to [0,1] if bounded.  Samples aged 1-100. ILLUSTRATIVE sizes."""
    rng = np.random.default_rng(seed)
    g = rng.uniform(0, 1, n_feat) if ground is None else ground

    def sample(n):
        ages = rng.integers(1, 101, n)
        Xs = np.empty((n, n_feat))
        for i, a in enumerate(ages):
            x = g + rng.normal(0, 0.01, n_feat) * 0       # ground state (no extra initial noise here)
            steps = rng.normal(0, sigma, (a, n_feat))
            if bounded:
                for s in steps:
                    x = np.clip(x + s, 0, 1)
            else:
                x = x + steps.sum(0)
            Xs[i] = x
        return Xs, ages.astype(float)
    return g, sample


def epigenetics():
    out = {}
    # 1) pure-noise clocks
    for bounded in (True, False):
        g, sample = stochastic_clock(bounded=bounded)
        Xtr, ytr = sample(300)
        Xva, yva = sample(150)
        Xte, yte = sample(300)
        best = None
        for lam in [0.01, 0.03, 0.1, 0.3, 1.0]:
            f, nz = elastic_net(Xtr, ytr, lam)
            err = np.mean((f(Xva) - yva) ** 2)
            if best is None or err < best[0]:
                best = (err, lam, f, nz)
        pred = best[2](Xte)
        out["bounded" if bounded else "unbounded"] = dict(
            pearson_r_test=float(np.corrcoef(pred, yte)[0, 1]),
            median_abs_err_yr=float(np.median(np.abs(pred - yte))), lambda_=best[1], n_nonzero=best[3])
        if bounded:
            clock_pred, clock_age = pred, yte
            # "intervention": individuals whose noise per year was 30% lower (same ground state)
            rng = np.random.default_rng(9)
            Xs = []
            for a in yte:
                x = g.copy()
                for s in rng.normal(0, 0.7 * 0.01, (int(a), g.size)):
                    x = np.clip(x + s, 0, 1)
                Xs.append(x)
            pred_slow = best[2](np.array(Xs))
            fit = np.polyfit(yte, pred_slow, 1)
            out["bounded"]["noise_x0.7_predicted_age_slope"] = float(fit[0])
            out["bounded"]["coef_sign_vs_ground_state"] = "features near 0 drift up, near 1 drift down: regression to the mean"
    # 2) per-CpG two-state model (Tarkhov et al.) and entropy
    #    each CpG: methylation fraction beta(t) = beq + (q - beq) e^{-k t}, k = pm + pd, beq = pm/k.
    #    finite cells (n=200 per sample) -> binomial sampling.  Rates per year ILLUSTRATIVE U(0, 0.01).
    rng = np.random.default_rng(3)
    n_cpg = 2000
    q = rng.choice([0.02, 0.98], n_cpg) * 1.0                   # sites start near fully un/methylated
    pm, pd = rng.uniform(0, 0.01, n_cpg), rng.uniform(0, 0.01, n_cpg)
    k, beq = pm + pd, pm / (pm + pd)
    ages = np.arange(0, 101, 5.0)
    beta = beq + (q - beq) * np.exp(-np.outer(ages, k))
    obs = rng.binomial(200, beta) / 200.0
    H = lambda b: -(np.clip(b, 1e-9, 1) * np.log2(np.clip(b, 1e-9, 1)) +
                    np.clip(1 - b, 1e-9, 1) * np.log2(np.clip(1 - b, 1e-9, 1)))
    ent = H(obs).mean(1)
    # symbolic check of the two-state solution
    t, a_, b_, q_ = sp.symbols("t p_m p_d q", positive=True)
    sol = a_ / (a_ + b_) + (q_ - a_ / (a_ + b_)) * sp.exp(-(a_ + b_) * t)
    ode_ok = sp.simplify(sp.diff(sol, t) - (a_ * (1 - sol) - b_ * sol)) == 0
    out["two_state"] = dict(ode_solution_verified=bool(ode_ok),
                            mean_entropy_bits_at={str(int(a)): float(e) for a, e in zip(ages, ent) if a in (0, 25, 50, 75, 100)},
                            rates_note="ILLUSTRATIVE per-year flip rates U(0,0.01)")
    return out, (clock_age, clock_pred, ages, ent)


# ============================================================ (d) mtDNA drift
def mtdna(N=1000, tau_days=10.0, n_traj=400000, max_years=120, thresh=0.6):
    """Neutral Wright-Fisher drift of a single new mutant among N mtDNA copies in a post-mitotic cell,
    one 'generation' = one turnover time tau.  N, tau ILLUSTRATIVE.  Returns the distribution of first-passage
    time to heteroplasmy >= thresh (60% used as an ILLUSTRATIVE biochemical threshold)."""
    gens_per_yr = 365.0 / tau_days
    G = int(max_years * gens_per_yr)
    k = np.ones(n_traj, int)
    idx = np.arange(n_traj)
    hit = np.full(n_traj, np.inf)
    for g in range(1, G + 1):
        k = RNG.binomial(N, k / N)
        reached = k >= thresh * N
        hit[idx[reached]] = g / gens_per_yr
        keep = (k > 0) & ~reached
        k, idx = k[keep], idx[keep]
        if k.size == 0:
            break
    T = hit[np.isfinite(hit)]
    p_reach = T.size / n_traj
    # cells: new mutant copies arise at rate nu per cell per year; F(t) = 1 - exp(-nu * int_0^t P(T<=t-s) ds)
    tg = np.linspace(0, 100, 1001)
    cdf = np.array([(T <= t).sum() / n_traj for t in tg])
    cum = np.concatenate([[0], np.cumsum(0.5 * (cdf[1:] + cdf[:-1]) * np.diff(tg))])
    nu = 1.0           # ILLUSTRATIVE: one new mutant copy (of a deleterious class) per cell per year
    F = 1 - np.exp(-nu * cum)
    return dict(N=N, tau_days=tau_days, threshold=thresh, nu_per_cell_per_year=nu,
                note="N, tau, nu ILLUSTRATIVE",
                p_single_mutant_reaches_threshold=p_reach, theory_p_reach_approx=1 / (thresh * N),
                median_years_to_threshold_given_reached=float(np.median(T)),
                frac_cells_over_threshold={str(a): float(np.interp(a, tg, F)) for a in (12, 20, 40, 60, 80)},
                frac_cells_dog_like_nu_x5_at_12=float(1 - np.exp(-5 * nu * np.interp(12, tg, cum))),
                ), (tg, F)


# ============================================================ (e) aggregation with clearance
def aggregation():
    """Moment equations (Knowles et al. 2009 structure) with first-order clearance c of aggregates and
    monomer m held fixed by synthesis/degradation balance:
        dP/dt = kn m^nc + km M - c P          (P: number of aggregates; km: fragmentation)
        dM/dt = 2 kp m P - c M                (M: mass in aggregates; kp: elongation)
    Linear system: eigenvalues -c +- kappa, kappa = sqrt(2 kp km m).  Runaway iff kappa > c.
    Clearance declines with age: c(t) = c0 (1 - t/150)  (ILLUSTRATIVE).  Units dimensionless/years."""
    kn, nc, kp, km, m = 1e-6, 2, 0.5, 0.25, 1.0
    kappa = np.sqrt(2 * kp * km * m)
    lam_ = sp.symbols("lambda")
    c_, k1, k2 = sp.symbols("c k_m g", positive=True)     # g = 2 kp m
    Msym = sp.Matrix([[-c_, k1], [k2, -c_]])
    eig = [sp.simplify(e) for e in Msym.eigenvals()]
    out = dict(kappa=kappa, eigenvalues_symbolic=[str(e) for e in eig],
               params=dict(kn=kn, nc=nc, kp=kp, km=km, m=m, note="ILLUSTRATIVE"))
    t_eval = np.linspace(0, 120, 1201)
    curves = {"age": t_eval}
    for name, c0, mscale in [("baseline", 1.0, 1.0), ("clearance_x1.5", 1.5, 1.0), ("monomer_x0.5", 1.0, 0.5)]:
        mm = m * mscale
        cfun = lambda t: c0 * max(1 - t / 150.0, 0.0)
        rhs = lambda t, y: [kn * mm ** nc + km * y[1] - cfun(t) * y[0], 2 * kp * mm * y[0] - cfun(t) * y[1]]
        sol = solve_ivp(rhs, (0, 120), [0.0, 0.0], t_eval=t_eval, rtol=1e-8, atol=1e-14)
        Mt = sol.y[1]
        kap = np.sqrt(2 * kp * km * mm)
        t_cross = 150 * (1 - kap / c0) if kap < c0 else 0.0
        curves[name] = Mt
        out[name] = dict(c0=c0, kappa=kap, age_clearance_falls_below_kappa=t_cross,
                         M_at={str(a): float(np.interp(a, t_eval, Mt)) for a in (40, 60, 80, 100, 120)})
    return out, curves


def main():
    res = {}
    tel, D = telomeres()
    res["telomeres"] = tel
    som, (ah, mh, ad, md) = somatic_mutations()
    res["somatic_mutations"] = som
    epi, (cage, cpred, eages, ent) = epigenetics()
    res["epigenetic"] = epi
    mt, (tg, F) = mtdna()
    res["mtdna"] = mt
    agg, curves = aggregation()
    res["aggregation"] = agg
    save_results("clocks", res)

    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    ax[0, 0].hist(D, bins=40, color="C0")
    ax[0, 0].set(xlabel="Divisions to senescence", ylabel="Cells", title="(a) Telomere counter (ILLUSTRATIVE L0, Lc)")
    ax[0, 1].plot(ah, mh.mean(0), label="human, 47/yr (Cagan 2022)")
    ax[0, 1].plot(ad, md.mean(0), label="dog, PREDICTED ELB/14 yr")
    ax[0, 1].axhline(3206.4, ls=":", c="k", label="end-of-lifespan burden 3206")
    ax[0, 1].set(xlabel="Age (years)", ylabel="Substitutions per crypt genome", title="(b) Somatic mutations")
    ax[0, 1].legend(fontsize=7)
    ax[0, 2].plot(cage, cpred, ".", ms=3)
    ax[0, 2].plot([0, 100], [0, 100], "k:")
    ax[0, 2].set(xlabel="True age (synthetic)", ylabel="Clock prediction",
                 title=f"(c) Clock trained on PURE NOISE, r={epi['bounded']['pearson_r_test']:.2f}")
    ax[1, 0].plot(eages, ent, "o-")
    ax[1, 0].set(xlabel="Age", ylabel="Mean methylation entropy (bits)", title="(c) Epigenetic entropy (ILLUSTRATIVE)")
    ax[1, 1].plot(tg, F)
    ax[1, 1].set(xlabel="Age (years)", ylabel="Fraction of cells above 60% mutant",
                 title="(d) mtDNA clonal drift (ILLUSTRATIVE N, tau)")
    for k in ["baseline", "clearance_x1.5", "monomer_x0.5"]:
        ax[1, 2].semilogy(curves["age"], np.maximum(curves[k], 1e-12), label=k)
    ax[1, 2].set(xlabel="Age (years)", ylabel="Aggregate mass (arb.)", title="(e) Aggregation vs clearance (ILLUSTRATIVE)")
    ax[1, 2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/clocks.png", dpi=140); plt.close(fig)
    save_csv("clocks_telomere_hayflick", pd.DataFrame(dict(divisions=D)), "Simulated Hayflick counts, ILLUSTRATIVE L0/Lc.")
    save_csv("clocks_somatic", pd.DataFrame(dict(age_h=ah, human=mh.mean(0), age_d=ad, dog_predicted=md.mean(0))),
             "Cagan 2022 rates (human 47/yr) and PREDICTED dog rate ELB/14; https://pmc.ncbi.nlm.nih.gov/articles/PMC9021023")
    save_csv("clocks_noise_clock", pd.DataFrame(dict(age=cage, predicted=cpred)), "Clock trained on synthetic pure noise.")
    save_csv("clocks_entropy", pd.DataFrame(dict(age=eages, entropy_bits=ent)), "ILLUSTRATIVE two-state model.")
    save_csv("clocks_mtdna", pd.DataFrame(dict(age=tg, frac_cells=F)), "ILLUSTRATIVE WF drift N=1000, tau=10 d.")
    save_csv("clocks_aggregation", pd.DataFrame(curves), "ILLUSTRATIVE aggregation model with declining clearance.")
    print("clocks done")


if __name__ == "__main__":
    main()
