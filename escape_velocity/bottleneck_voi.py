"""Part 5: bottleneck model and value-of-information ranking.

EVERY numeric input in this file is a JUDGEMENT input (not data), chosen by the stream
orchestrator and stated in findings.md. The point of the model is transparency:
anyone can change an input and see whether the ranking changes. Robustness to the
inputs is tested by randomly perturbing all of them.

Two outcome measures, both functions of six binary unknowns U_i (1 = favourable):
  M1 = F + exp(.) years from 2026 until a therapy is shown in humans to slow aging and is approved
       for use (log-additive model):  log M1 = log Y0 + sum_i log(k_i^{U_i} * j_i^{1-U_i}) + eps
  M2 = long-run ceiling on the rate of mortality improvement, as a fraction of the
       escape-velocity rate b (rho/b; 1.0 = escape velocity):  log M2 = log Z0 + sum_i ... + eps
"""
import numpy as np
import pandas as pd

from common import C, DATA, INK2, MUTED, finish, plt, save_json

rng = np.random.default_rng(42)
out = {}

U = pd.DataFrame([
    # key, question, prior p(favourable), M1 factor if yes, if no, M2 factor if yes, if no
    ("U1", "Few upstream causes (vs many independent ones)", 0.30, 0.70, 1.20, 2.0, 0.7),
    ("U2", "Mouse lifespan results translate to humans",      0.40, 0.60, 1.40, 1.3, 0.9),
    ("U3", "Biomarkers can stand in for lifespan in trials",  0.35, 0.50, 1.30, 1.1, 1.0),
    ("U4", "Partial reprogramming is safe in whole bodies",   0.25, 0.85, 1.00, 2.0, 0.9),
    ("U5", "Regulators approve an aging indication",          0.50, 0.70, 1.30, 1.1, 0.95),
    ("U6", "No damage type is beyond repair (weakest link)",  0.30, 1.00, 1.00, 2.5, 0.6),
], columns=["key", "question", "p", "m1_yes", "m1_no", "m2_yes", "m2_no"])
F, Y0, Z0, SIG = 5.0, 17.0, 0.15, 0.4    # JUDGEMENT: floor years (one trial + review); median extra years; median ceiling rho/b; residual log-sd


def analyse(tab, n=200_000, rng=rng):
    Ub = rng.random((n, len(tab))) < tab.p.values
    eps1, eps2 = rng.normal(0, SIG, n), rng.normal(0, SIG, n)
    l1 = np.log(Y0) + np.where(Ub, np.log(tab.m1_yes.values), np.log(tab.m1_no.values)).sum(1) + eps1
    l2 = np.log(Z0) + np.where(Ub, np.log(tab.m2_yes.values), np.log(tab.m2_no.values)).sum(1) + eps2
    res = []
    for i, k in enumerate(tab.key):
        sw1 = l1[Ub[:, i]].mean() - l1[~Ub[:, i]].mean()
        sw2 = l2[Ub[:, i]].mean() - l2[~Ub[:, i]].mean()
        p = tab.p.values[i]
        # first-order variance share (exact for this additive-in-log model; MC estimate here)
        s1 = p * (1 - p) * sw1 ** 2 / l1.var()
        s2 = p * (1 - p) * sw2 ** 2 / l2.var()
        res.append(dict(key=k, swing_M1_ratio=np.exp(sw1), swing_M2_ratio=np.exp(sw2),
                        var_share_M1=s1, var_share_M2=s2, voi_score=0.5 * (s1 + s2)))
    res = pd.DataFrame(res)
    return res, l1, l2, Ub


res, l1, l2, Ub = analyse(U)
res = U[["key", "question", "p"]].merge(res, on="key").sort_values("voi_score", ascending=False)
res["rank"] = np.arange(1, len(res) + 1)
out["voi_ranking"] = res.to_dict(orient="records")
q = lambda a, pcts, add=0.0: {f"p{pp}": float(add + np.exp(np.percentile(a, pp))) for pp in pcts}
out["M1_years_to_first_validated_therapy_quantiles"] = q(l1, (10, 25, 50, 75, 90), F)
out["M2_ceiling_rho_over_b_quantiles"] = q(l2, (10, 25, 50, 75, 90))
out["P_M2_ge_1_escape_velocity_ceiling"] = float(np.mean(l2 >= 0))

# scenario reading (conditional medians) -- JUDGEMENT, not forecasts
scen = {
    "pessimistic: U2,U3,U5 unfavourable": (~Ub[:, 1]) & (~Ub[:, 2]) & (~Ub[:, 4]),
    "central: all draws": np.ones(len(l1), bool),
    "optimistic: U2,U3,U5 favourable": Ub[:, 1] & Ub[:, 2] & Ub[:, 4],
    "repair-ceiling: U1,U4,U6 favourable": Ub[:, 0] & Ub[:, 3] & Ub[:, 5],
}
out["scenarios"] = {k: dict(share_of_draws=float(m.mean()),
                            M1_median_years=float(F + np.exp(np.median(l1[m]))),
                            M1_p10_p90=[float(F + np.exp(np.percentile(l1[m], 10))), float(F + np.exp(np.percentile(l1[m], 90)))],
                            M2_median=float(np.exp(np.median(l2[m]))),
                            P_M2_ge_1=float(np.mean(l2[m] >= 0)))
                    for k, m in scen.items()}

# robustness: perturb every judgement input by up to +/-50% (factors toward 1 in log space)
ranks = {k: [] for k in U.key}
for _ in range(500):
    T = U.copy()
    T["p"] = np.clip(T.p * rng.uniform(0.5, 1.5, len(T)), 0.05, 0.95)
    for col in ("m1_yes", "m1_no", "m2_yes", "m2_no"):
        T[col] = np.exp(np.log(T[col]) * rng.uniform(0.5, 1.5, len(T)))
    r_, *_ = analyse(T, n=20_000, rng=rng)
    order = r_.sort_values("voi_score", ascending=False).key.tolist()
    for k in U.key:
        ranks[k].append(order.index(k) + 1)
rob = pd.DataFrame({k: dict(mean_rank=np.mean(v), share_top2=np.mean(np.array(v) <= 2),
                            share_top3=np.mean(np.array(v) <= 3)) for k, v in ranks.items()}).T
rob.index.name = "key"
rob = rob.reset_index()
out["robustness"] = rob.to_dict(orient="records")
full = res.merge(rob, on="key")
full.to_csv(f"{DATA}/fig_voi_ranking.csv", index=False)

fig, ax = plt.subplots(figsize=(8, 3.8))
y = np.arange(len(full))[::-1]
ax.barh(y + 0.18, full.var_share_M1, height=0.34, color=C[0], label="share of uncertainty in time to first therapy (M1)")
ax.barh(y - 0.18, full.var_share_M2, height=0.34, color=C[1], label="share of uncertainty in long-run ceiling (M2)")
ax.set_yticks(y); ax.set_yticklabels([f"{k}: {qq}" for k, qq in zip(full.key, full.question)], fontsize=8)
ax.set_xlabel("first-order variance share (higher = resolving it changes the outlook more)")
ax.set_title("Which unknown matters most? (judgement inputs)")
ax.legend(fontsize=7, loc="lower right")
finish(fig, "fig_voi_ranking", "JUDGEMENT model: every prior and effect size is an orchestrator judgement, not data; see findings.md.")

save_json("results_voi.json", out)
print(full.round(3).to_string())
print(out["M1_years_to_first_validated_therapy_quantiles"], out["M2_ceiling_rho_over_b_quantiles"], out["P_M2_ge_1_escape_velocity_ceiling"])
for k, v in out["scenarios"].items():
    print(k, {kk: (np.round(vv, 3) if not isinstance(vv, list) else np.round(vv, 1)) for kk, vv in v.items()})
