"""Part 4: causes of death and healthspan.

Quantitative inputs available from opened sources:
  O'Neill 2013 (UK, 5,095 deaths): cause of death by age group (<3 y vs >=3 y).
  McMillan 2024: hazard ratios by face shape and size x face shape.
  Banzato 2019: frailty index (33 deficits) and mortality in 401 dogs.
  Bannasch 2021: inbreeding and body weight vs insurance morbidity across breeds.
Age-specific cause-of-death *rates* are derived by combining O'Neill's age-group cause shares with the
Teng 2022 age-specific hazard (ASSUMPTION: cause mix within each age band is the same in both UK datasets).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import read, update_results, save_csv, FIG, load_results

c = read("oneill2013_causes_by_age.csv").fillna(0)
tot_u3, tot_o3 = 489, 4606   # totals in O'Neill Table 3
c["share_under3"] = c.under3_n / tot_u3
c["share_3plus"] = c.age3plus_n / tot_o3
c["ratio_3plus_over_under3"] = c.share_3plus / c.share_under3.replace(0, np.nan)
print(c[["cause", "share_under3", "share_3plus"]].round(3).to_string(index=False))

# cause-specific hazards at two ages: share x all-cause hazard (Teng 2022)
teng = read("teng2022_uk_dog_life_table.csv")
mu = -np.log(1 - teng.qx.values[:-1])   # drop open 20+ interval (q=1)
mu_young = float(np.mean(mu[0:3]))         # ages 0-2
mu_old = float(np.mean(mu[3:17]))          # ages 3-16 (crude average)
c["hazard_under3_per_1000"] = 1000 * c.share_under3 * mu_young
c["hazard_3plus_per_1000"] = 1000 * c.share_3plus * mu_old

mc = read("mcmillan2024_size_face.csv")
out = dict(
    oneill2013_shares={r.cause: dict(under3=round(r.share_under3, 4), age3plus=round(r.share_3plus, 4)) for r in c.itertuples()},
    neoplastic_share_under3_vs_3plus=[float(c.loc[c.cause == "Neoplastic", "share_under3"].iloc[0]),
                                      float(c.loc[c.cause == "Neoplastic", "share_3plus"].iloc[0])],
    external_causes_share_under3=float(c[c.cause.isin(["Road traffic accident", "Trauma", "Behavioural abnormality", "Dangerous Dogs Act"])].share_under3.sum()),
    teng_mean_hazard_age0to2=mu_young, teng_mean_hazard_age3to16=mu_old,
    face_shape_HR_vs_mesocephalic={r.group: r.HR for r in mc[mc.group_type == "face"].itertuples()},
    frailty_Banzato2019={"n_dogs": 401, "deficits": 33, "HR_FI_0.2to0.4_vs_lt0.2": 9.21, "HR_FI_gt0.4_vs_lt0.2": 18.06,
                         "FI_slope_per_year": 0.016, "mean_FI_survivors": 0.12, "mean_FI_died": 0.29,
                         "url": "https://air.unimi.it/retrieve/handle/2434/717203/1431369/s41598-019-52585-9.pdf"},
    inbreeding_Bannasch2021=dict(n_breeds_inbreeding=227, mean_Fadj=0.249, spearman_inbreeding_morbidity=0.29,
                                 spearman_weight_morbidity=0.41, predicted_morbidity_5kg_vs_60kg_increase_pct=37,
                                 brachycephalic_morbidity_vs_other_per_10000_DYAR=[1829, 1508],
                                 note="Published regression coefficients could not be reconciled with the paper's worked predictions; only correlations and the stated 37% prediction are used.",
                                 url="https://escholarship.org/content/qt5b61d7p0/qt5b61d7p0.pdf?t=rxtyal"),
    multimorbidity_DAP_2026=dict(claim="dogs with >=5 chronic conditions >2x as likely to die in a given year as dogs with 1",
                                 n_dogs=">50,000", source="https://www.dvm360.com/view/multiple-chronic-diseases-compound-to-shorten-dogs-lifespans-study-finds",
                                 note="secondary report of Foreman et al. GeroScience 2026; hazard ratios not seen"),
    fleming2011_size_pattern="larger breeds: more musculoskeletal, GI and especially cancer deaths; smaller: more metabolic (qualitative, via ScienceDaily)")
update_results("part4_causes", out)

keep = c[(c.share_under3 + c.share_3plus) > 0.04].sort_values("share_3plus")
fig, ax = plt.subplots(figsize=(7.5, 5))
y = np.arange(len(keep))
ax.barh(y - 0.2, 100 * keep.share_under3, 0.4, label="died < 3 years (n=489)")
ax.barh(y + 0.2, 100 * keep.share_3plus, 0.4, label="died ≥ 3 years (n=4,606)")
ax.set_yticks(y); ax.set_yticklabels(keep.cause, fontsize=8); ax.set_xlabel("% of deaths")
ax.set_title("UK dogs: young dogs die of accidents and behaviour,\nolder dogs of cancer, joints and nerves (O'Neill 2013)")
ax.legend(fontsize=8); plt.tight_layout(); plt.savefig(f"{FIG}/causes_by_age.png", dpi=150); plt.close()
save_csv(c, "causes_by_age.csv", ["Series for figures/causes_by_age.png; counts from O'Neill 2013 Table 3 https://rvc-repository.worktribe.com/OutputFile/1409027",
                                  "hazard columns = share x Teng 2022 mean all-cause hazard (ASSUMPTION: same cause mix across datasets)"])
print("done")
