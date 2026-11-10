"""Writes results.json['headline']: the human parameters other streams should reuse, with units, data year
and source. Run after fit_mortality_laws.py, le_sensitivity.py, limits_extremes.py, no_aging.py."""
import json
from common import RESULTS, SSA_URL, update_results

r = json.load(open(RESULTS))
f = r["fits"]
h = {"data": "US SSA period life table, data year 2023 (2026 Trustees Report)", "source_url": SSA_URL,
     "cross_check": "NCHS Data Brief 521: US 2023 e0 male 75.8, female 81.1, total 78.4 "
                    "(https://cdc.gov/nchs/products/databriefs/db521.htm)",
     "hazard_definition": "mu(x) = c + a*exp(b*x), x = age in years, mu per year",
     "recommended_model": "gompertz_makeham fitted to ages 30-95",
     "units": {"a": "per year", "b": "per year of age", "c": "per year", "doubling_time": "years", "e0": "years"}}
for sex in ["male", "female"]:
    gm = f["models"][sex]["gompertz_makeham"]
    g = f["models"][sex]["gompertz"]
    o = f["observed"][sex]
    h[sex] = {
        "e0_period_2023": o["e0_published"], "e65_period_2023": o["e65_published"],
        "median_age_at_death": round(o["median_age_at_death"], 2), "modal_adult_age_at_death": o["modal_age_at_death_adult"],
        "GM_a": gm["params"]["a"], "GM_a_ci95": gm["ci95"]["a"],
        "GM_b": gm["params"]["b"], "GM_b_ci95": gm["ci95"]["b"],
        "GM_c_baseline_hazard": gm["params"]["c"], "GM_c_ci95": gm["ci95"]["c"],
        "GM_doubling_time_years": gm["doubling_time_years"], "GM_doubling_time_ci95": gm["doubling_time_ci95"],
        "Gompertz_40_95_a": g["params"]["a"], "Gompertz_40_95_b": g["params"]["b"],
        "Gompertz_40_95_doubling_time_years": g["doubling_time_years"],
        "empirical_loglinear_slope_40_90": o["loglinear_slope_40_90"],
        "empirical_doubling_time_40_90_years": o["doubling_time_40_90_years"],
        "observed_hazard_age20": o["hazard_age20_per_year"], "observed_hazard_age30": o["hazard_age30_per_year"],
        "keyfitz_entropy_H": r["le_sensitivity"][sex]["keyfitz_entropy_H"],
        "no_aging_life_expectancy_hazard_frozen_at_30": r["no_aging"][sex]["age30"]["A_constant_from_birth"]["life_expectancy"],
        "no_aging_life_expectancy_hazard_frozen_at_20": r["no_aging"][sex]["age20"]["A_constant_from_birth"]["life_expectancy"],
    }
h["doubling_time_note"] = ("The mortality doubling time depends on method: 6.8-7.1 y (GM, slope of the senescent part "
                           "after removing c), 7.6-8.3 y (pure Gompertz 40-95), 7.9-8.7 y (raw log-linear 40-90). "
                           "Quote ~7-9 years, not a single number.")
update_results("headline", h)
print(json.dumps(h, indent=1))
