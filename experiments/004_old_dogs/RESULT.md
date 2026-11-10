# Experiment 004: result

**The doubling law holds for dogs up to about age 15 and then overstates risk: by 19 percent at age 16 and by 50 to 95 percent at ages 17 to 19.** Hidden frailty of the size that explains this does not break the fingerprint test, so the advice from experiment 001 to follow dogs until death stands. Two predictions were confirmed and two refuted.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | A model with levelling-off fits the UK counts significantly better | likelihood-ratio statistic 100.6, p about 10⁻²³ | Confirmed |
| P2 | The plain law fitted on ages 2 to 13 forecasts ages 14 to 16 within 15% | within 5% at 14 and 15; 19% too high at 16 | **Refuted** |
| P3 | The same fit overstates the death rate at age 18 by more than 20% | 80% too high | Confirmed |
| P4 | With hidden frailty, a 30% risk cut is misread as a slope change in over 10% of 1,300-dog, 10-year trials | 6% | **Refuted** |

P2 failed by a narrow margin: the departure starts one year earlier than predicted. P4 failed in the useful direction: the slope test is sturdier than expected.

## Findings

**The plain law, fitted to ages 2 to 13, against what was observed in the UK table**

| Age | 14 | 15 | 16 | 17 | 18 | 19 |
| --- | --- | --- | --- | --- | --- | --- |
| Forecast divided by observed | 1.01 | 1.05 | 1.19 | 1.50 | 1.80 | 1.95 |

The US table, which stops at age 16, shows the same thing: 0.98, 1.06 and 1.22 at ages 14, 15 and 16.

**Risk keeps rising at the oldest ages, only more slowly.** Observed yearly death rates at ages 16 to 19 are 0.90, 0.99, 1.15 and 1.48. "Levelling off" overstates it; "bending" is accurate.

**A frailty model fits.** If each dog carries a hidden personal risk factor with variance 0.10, every dog's own risk doubles every 2.00 years (against 2.17 for the population), and the population curve bends because the frailest dogs die first.

**Effect on the fingerprint test.** Simulated trials with that hidden frailty, analysed without it:

| Design | No effect: flagged | 30% risk cut: flagged | 50% risk cut: flagged | 20% slowing: recognised |
| --- | --- | --- | --- | --- |
| 580 dogs, 3 years | 3% | 4% | 6% | 8% |
| 1,300 dogs, 4 years | 5% | 6% | 7% | 30% |
| 1,300 dogs, 10 years | 5% | 6% | 10% | 87% |
| 5,000 dogs, 3 years | 4% | 4% | 5% | 39% |

False alarms stay near the designed 5 percent, and a long trial still recognises a true slowing of aging 87 percent of the time.

## What this means

1. **The doubling law is safe to use for dogs up to about 15.** Results in this repository that concern dogs older than that, including the oldest ages in the dog-years table, lean on a law that overstates their risk.
2. **"Follow dogs until death" survives hidden frailty of the fitted size.**
3. **It would not survive much larger hidden differences between dogs.** The check below found that tripling the frailty variance pushes false alarms to about 20 percent for a 50 percent risk cut. Breed differences could plausibly be that large. A long trial should model differences between dogs, at least by breed and size.

## Independent check

A separate agent refitted both models with its own code and reproduced every parameter, the likelihood-ratio statistic, the forecast ratios and the trial numbers. It verified the formula for frailty among dogs alive at entry by brute-force simulation. It found no errors, and corrected three points of emphasis, all adopted above:

- Risk at 17 to 19 is still rising, so "levels off" was too strong. Ages 18 and 19 rest on 218 and 78 deaths.
- The 50 to 95 percent overstatement is for the law fitted to ages 2 to 13. Fitted to all ages, the plain law overstates ages 17 to 19 by about 35 to 71 percent.
- Our 10 percent false-alarm figure for a 50 percent cut in the 10-year design was probably a high draw; the checker's pooled estimate is about 7 percent. With frailty variance 0.3 it found 20 percent.

## A caveat on the data

The UK table is built from deaths recorded in 2016 to 2020, not by following dogs through life. The bend at old ages mixes true aging with breed mix (small, long-lived breeds dominate the oldest deaths), changes in breed popularity over time, and errors in recorded ages. The data show that this population's death rate grows more slowly at the oldest ages. They do not show that an individual dog's risk bends.

## Files

`old_dogs.py` (about 40 seconds), `results.json`, `PREDICTION.md`.
