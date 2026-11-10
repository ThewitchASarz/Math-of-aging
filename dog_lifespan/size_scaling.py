"""Part 2: why big dogs die young, and how that squares with across-species scaling.

Inputs (see CSV headers for URLs): Montoya 2023 e0 by size group, O'Neill 2013 bodyweight coefficients,
McMillan 2024 size-class medians, de Magalhaes 2007 and Fan 2016 allometric exponents,
plus the size-class Gompertz calibration from fit_mortality.py.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp

from common import read, update_results, save_csv, load_results, FIG, life_expectancy, LN2

size = read("montoya2023_size_e0.csv")
oneill = read("oneill2013_bodyweight_coef.csv")
lit = read("scaling_literature.csv").set_index("group")
cal = pd.DataFrame(load_results()["part1_mortality_fits"]["size_class_calibration"]["rows"])
out = {}

# ---- 1. within-dog allometry from Montoya size groups
# ASSUMPTION: representative adult mass of a group = geometric mean of its weight bounds
# (toy lower bound 2 kg and giant upper bound 70 kg are this stream's assumptions, not in the paper).
size["M_rep_kg"] = np.sqrt(size.kg_lo * size.kg_hi)
X = np.log(size.M_rep_kg.values); Y = np.log(size.e0.values)
slope, intercept = np.polyfit(X, Y, 1)
lin_slope, lin_int = np.polyfit(size.M_rep_kg.values, size.e0.values, 1)
# sensitivity to the two assumed bounds
sens = []
for lo_toy in (1.0, 2.0, 3.0):
    for hi_giant in (60.0, 70.0, 90.0):
        m = np.sqrt(np.r_[lo_toy, size.kg_lo[1:]] * np.r_[size.kg_hi[:-1], hi_giant])
        sens.append(np.polyfit(np.log(m), Y, 1)[0])
# exclude toy (non-monotone: toy < small)
slope_no_toy = np.polyfit(X[1:], Y[1:], 1)[0]
out["within_dog_allometry_montoya"] = dict(
    exponent_lnE0_on_lnM=slope, exponent_range_over_bound_assumptions=[min(sens), max(sens)],
    exponent_excluding_toy=slope_no_toy,
    linear_years_per_kg=lin_slope,
    Fan2016_breed_exponent=float(lit.loc["dog_breeds_within_species", "b"]),
    note="5 size-group points only; representative masses are assumptions (geometric means of bounds).")

# ---- 2. O'Neill 2013 individual-level adjusted coefficients -> years per kg
mid = np.r_[5, 15, 25, 35, 45, 60.0]   # ASSUMPTION: band midpoints, top band >=50 kg set to 60 kg
w = 1 / ((oneill.hi - oneill.lo).replace(0, np.nan) / 3.92) ** 2
w = w.fillna(w.max())
yrs_per_kg = np.polyfit(mid, oneill.coef_years.values, 1, w=np.sqrt(w.values))[0]
out["oneill2013_years_per_kg"] = dict(slope_years_per_kg=yrs_per_kg,
                                      difference_50plus_vs_under10_years=-4.0,
                                      note="Weighted line through published band coefficients; band midpoints assumed. Dogs >=3 y only.")

# ---- 3. elasticity of the Gompertz slope with mass (scenario A of the calibration)
elas_b = np.polyfit(np.log(size.M_rep_kg.values[1:]), np.log(cal.A_b.values[1:]), 1)[0]
out["gompertz_slope_vs_mass"] = dict(
    elasticity_dlnb_dlnM_small_to_giant=elas_b,
    b_small=cal.A_b[1], b_giant=cal.A_b[4], MRDT_small=cal.A_MRDT[1], MRDT_giant=cal.A_MRDT[4],
    note="If size acts only through the rate of aging (Kraus 2013), the doubling time falls from ~2.2 y (small) to ~1.4 y (giant).")

# ---- 4. years attributable to size
wts = size.mean_annual_n.values / size.mean_annual_n.sum()
ref = size.e0.max()                      # best group (small, 5.5-11 kg)
gap = ref - size.e0.values
attrib = float(np.sum(wts * gap))
attrib_vs_medium = float(np.sum(wts * np.clip(size.e0[size["size"] == "medium"].values[0] - size.e0.values, 0, None)))
out["years_attributable_to_size"] = dict(
    reference_group="small (5.5-11 kg), e0=%.2f" % ref,
    per_group_gap_years=dict(zip(size["size"], gap.round(2))),
    population_weighted_mean_years_lost=attrib,
    population_weighted_mean_years_lost_above_medium=attrib_vs_medium,
    giant_minus_small_years=float(size.e0.iloc[4] - size.e0.iloc[1]),
    fraction_of_small_dog_life=attrib / ref,
    weights="Montoya mean annual dogs per size group (purebreds only)",
    caveat="Size groups also differ in breed-specific disease and face shape; this is an upper-bound association, not a causal effect.")

# ---- 5. across-species vs within-species: maths
M = sp.symbols("M", positive=True)
a_s, b_s = float(lit.loc["mammals_no_bats", "a"]), float(lit.loc["mammals_no_bats", "b"])
tmax = a_s * (M * 1000) ** b_s          # mass in kg -> g
ratio_across = float((tmax.subs(M, 60) / tmax.subs(M, 3)))
ratio_within = float(60 ** slope / 3 ** slope)
ratio_within_fan = float(20 ** lit.loc["dog_breeds_within_species", "b"])
out["across_vs_within"] = dict(
    mammal_tmax_60kg_over_3kg=ratio_across, dog_e0_60kg_over_3kg_montoya=ratio_within,
    dog_ratio_from_Fan_exponent=ratio_within_fan,
    mammal_tmax_predicted_for_30kg_years=float(tmax.subs(M, 30)),
    kleiber_metabolic_exponent=0.75,
    explanation=("Across species, lifespan exponent ~+0.15 to +0.19 (de Magalhaes 2007). Within dogs, ~-0.1. "
                 "Both hold because they are slopes of different relationships: ln T_ij = alpha + beta_between*ln M_j "
                 "+ beta_within*(ln M_ij - ln M_j). beta_between reflects evolved maintenance (large species evolved "
                 "slower aging, b_species ~ M^-beta). beta_within reflects recent artificial selection for size via "
                 "growth pathways (IGF-1) without time to evolve matching maintenance; here b rises with mass."))
print(out)

# ---- figures
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].errorbar(size.M_rep_kg, size.e0, yerr=[size.e0 - size.lo, size.hi - size.e0], fmt="o", label="Montoya 2023 e0 by size group (US)")
mm = np.geomspace(2, 80, 100)
ax[0].plot(mm, np.exp(intercept) * mm ** slope, "-", label=f"fit: e0 ∝ M^{slope:.3f}")
ax[0].plot(mid[:], 13.0 + oneill.coef_years.values, "s", mfc="none",
           label="O'Neill 2013 (UK) adjusted band effect,\nplotted as 13.0 + coef (offset illustrative)")
ax[0].set_xscale("log"); ax[0].set_xlabel("adult body mass (kg, log scale; group masses assumed)"); ax[0].set_ylabel("life expectancy / longevity (years)")
ax[0].set_title("Within dogs: bigger = shorter-lived"); ax[0].legend(fontsize=7)
# across species
mg = np.geomspace(0.01, 5000, 200)
ax[1].loglog(mg, a_s * (mg * 1000) ** b_s, "-", label=f"mammals (no bats): tmax = {a_s}·M_g^{b_s} (de Magalhães 2007)")
ax[1].loglog(size.M_rep_kg, size.e0, "o", label="dogs: mean life expectancy by size")
ax[1].loglog(mm, np.exp(intercept) * mm ** slope, "-", color="C1")
ax[1].set_xlabel("adult body mass (kg)"); ax[1].set_ylabel("years")
ax[1].set_title("Across species bigger = longer-lived; dogs run the other way"); ax[1].legend(fontsize=7)
ax[1].text(0.02, 0.02, "Note: species line is MAXIMUM lifespan; dog points are MEAN life expectancy.\nCompare slopes, not heights.", transform=ax[1].transAxes, fontsize=7)
plt.tight_layout(); plt.savefig(f"{FIG}/size_vs_lifespan.png", dpi=150); plt.close()
save_csv(size[["size", "M_rep_kg", "e0", "lo", "hi"]], "size_vs_lifespan_dogs.csv",
         ["Series for figures/size_vs_lifespan.png; e0 from Montoya 2023; M_rep_kg = geometric mean of group bounds (ASSUMPTION)"])
save_csv(pd.DataFrame(dict(M_kg=mg, mammal_tmax=a_s * (mg * 1000) ** b_s)), "size_vs_lifespan_mammal_line.csv",
         ["de Magalhaes 2007 mammals excluding bats tmax = 3.34*M_g^0.193"])

fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(cal["size"], cal.A_MRDT)
ax.set_ylabel("mortality-rate doubling time (years)"); ax.set_title("If size acts only on the rate of aging (calibrated)")
for i, v in enumerate(cal.A_MRDT): ax.text(i, v + 0.03, f"{v:.2f}", ha="center")
plt.tight_layout(); plt.savefig(f"{FIG}/mrdt_by_size.png", dpi=150); plt.close()
save_csv(cal[["size", "A_b", "A_MRDT", "e0_observed"]], "mrdt_by_size.csv", ["Series for figures/mrdt_by_size.png (calibration scenario A)"])

update_results("part2_size", out)
print("done")
