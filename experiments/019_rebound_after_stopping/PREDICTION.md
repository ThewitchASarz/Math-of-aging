# Experiment 019: prediction

Written before any code for this experiment was run.

## Question

In a trial that doses for one year and watches for three, what do the standard analyses report when death risk rebounds above the untreated level after dosing stops, or when only part of the benefit lasts, and does the mix of dogs that survive the dosing year bias the after-dosing comparison?

## Why it matters

Experiments 003, 006 and 016 simulated two things that can happen after a treatment stops: the benefit vanishes, or it lasts in full. Published switch experiments in animals show two more (`papers/literature-check.md`). Part of the benefit can last (Hahn et al. 2019, mice). Or risk can rise above the untreated level for a while (McCracken et al. 2020, flies). A trial that stops dosing needs to know what its analyses would say in those worlds.

There is also a published warning about comparing groups after the treatment period (Bartlett et al. 2020). The dogs alive after the dosing year are survivors. If the treatment kept frail dogs alive, the treated survivors are frailer on average than the placebo survivors, and the after-dosing comparison can look worse for the treatment than it is. This experiment measures how large that is in this design.

Nothing here says any treatment rebounds. No rebound has been reported in dogs in the sources this project has read. The sizes used are assumptions, chosen to see what a trial could detect.

## Method, fixed in advance

**The trial.** The rapamycin trial's published design, as in experiment 016, Part A: 580 dogs, half treated; death risk follows the doubling law with level 0.0066 and slope 0.29 per year; half the dogs enter at age 7 to 8, a quarter at 8 to 9, a quarter at 9 to 10; 12 months of dosing; 36 months of follow-up. These inputs are taken from `interventions/triad_published_inputs.py` (Coleman et al. 2025, as recorded in this repository). No dropout.

**The model.** A treated dog's risk is the untreated risk at its age times a multiplier. The multiplier is one number during the dosing year and other numbers afterwards. Placebo dogs have multiplier 1 throughout.

**True effects simulated (all sizes are ASSUMPTIONS):**

| Name | Dosing year | Months 12 to 18 | Months 18 to 36 |
| --- | --- | --- | --- |
| No effect | 1 | 1 | 1 |
| Fades (as in experiment 016) | 0.5 | 1 | 1 |
| Lasting rebound | 0.5 | 1.5 | 1.5 |
| Short rebound | 0.5 | 2 | 1 |
| Small lasting rebound | 0.8 | 1.25 | 1.25 |
| Rebound with no benefit first | 1 | 1.5 | 1.5 |
| Half the benefit lasts | 0.5 | 0.75 | 0.75 |
| Half of a smaller benefit lasts | 0.7 | 0.85 | 0.85 |

**Hidden differences between dogs.** For the survivor question, each dog also gets a hidden personal risk factor (average 1, gamma-shaped) that multiplies its risk for life. Two spreads (variances): 0.10, the tentative value fitted to the UK table in experiment 004, and 0.30, an ASSUMPTION used as the wide case in experiments 004, 008 and 012. The "fades" and "no effect" worlds are run with each. In these worlds the treatment truly does nothing to any dog after the dosing year.

**Analyses of each simulated trial (all reused from earlier experiments):**
1. The any-difference test: a two-sided log-rank test at the 5% level over all three years. We also record which group died more when it fires.
2. The two-period test of experiment 006: is the risk ratio after dosing different from 1? We also record whether the estimate points above 1 (worse after stopping) or below 1 (benefit lasts).
3. The slope test of experiment 001. Reported, not scored: experiment 016 already showed it must not be used on this design.

**Side calculations, reported and not scored.** The after-dosing risk ratio between the surviving treated and surviving placebo dogs, by direct calculation, for the "fades" world with hidden differences. The same for an older group (entry ages 10 to 14, UK death rates), where more dogs die in the dosing year.

1,000 simulated trials per case, fixed seeds.

## Predictions

1. **P1. A rebound makes the main test read a helpful-while-given treatment as harmful.** In the "lasting rebound" world (risk halved in the dosing year, then 1.5 times the untreated level), the any-difference test fires in 20 to 45% of trials. In at least 90% of the trials where it fires, the treated group is the one that died more.
2. **P2. The two-period test sees a lasting rebound but not a short one.** It reports an after-dosing difference pointing above 1 in 60 to 85% of trials in the "lasting rebound" world, and in fewer than 35% in the "short rebound" world.
3. **P3. A benefit that lasts in part is seen by the main test more often than it is recognised as lasting.** In the "half the benefit lasts" world, the any-difference test fires in 50 to 75% of trials, and the two-period test reports a lasting benefit in 20 to 45%.
4. **P4. Survivor bias is too small to matter in this design.** In the "fades" world with hidden differences of spread 0.30, the two-period test reports an after-dosing difference in fewer than 8% of trials (5% is chance), and the directly calculated after-dosing risk ratio between survivors, averaged over months 12 to 36, is below 1.03.

## What would change our mind

If P1 fails because the test rarely fires, a rebound of this size hides inside the trial and would go unnoticed, which is a different warning. If P2's first part fails downward, the two-period test of experiment 006 is not enough to catch a rebound and the trial would need more dogs or a different analysis. If P4 fails, the two-period test needs a correction for survivor mix before it is used on any stop-and-watch trial, and experiments 006 and 016 need a caveat.

These are computational experiments. They show what the mathematics and the published design imply. They do not replace a real trial, they say nothing about whether any treatment works or rebounds, and they do not question the trial's own power statement, which answers a different question.
