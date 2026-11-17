# Experiment 016: prediction

Written before any code for this experiment was run.

## Question

When the two dog trials now running are simulated with their own published designs, what could each one detect, and what share of placebo dogs dies in each?

## Why it matters

Experiments 001, 003 and 010 used this project's own guesses for the trials: entry ages 7 to 12, the UK death rates, and in experiment 001 treatment for the whole trial. The literature check (`papers/literature-check.md`) found the real designs differ. The rapamycin trial (TRIAD) publishes its own death-rate inputs, enrols dogs aged 7 to 10, and doses for 12 months of its 36. Loyal's STAY study enrols dogs aged 10 or older and appears to dose throughout. This experiment repeats the earlier ones with those designs, so the numbers describe the trials as designed.

These simulations answer a different question from the trials' own power statements. The rapamycin trial's design paper is sized for a large, lasting cut in risk. Nothing here says that calculation is wrong.

## Method, fixed in advance

**Part A: the rapamycin trial as published.** 580 dogs, half treated. Death risk follows the doubling law with level 0.0066 and slope 0.29 per year. Half the dogs enter at age 7 to 8, a quarter at 8 to 9, a quarter at 9 to 10. Treatment lasts 12 months; every dog is followed for 36 months. These inputs are published in the design paper (Coleman et al. 2025) and are taken here from `interventions/triad_published_inputs.py` and `papers/literature-check.md`. **This run could not open the paper itself**, so they rest on the earlier reading recorded in this repository.

True effects simulated: none; a 20% and a 50% risk cut during the dosing year that then fades; a 20% and a 38% risk cut that lasts all three years; aging slowed by 20% during the dosing year with the gain kept; aging stopped during the dosing year with the gain kept.

Three analyses of each simulated trial:
1. An any-difference test (two-sided log-rank test at the 5% level). ASSUMPTION: the paper names a Cox model; the exact test is ours.
2. The slope test of experiment 001.
3. The two-period test of experiment 006 (is there a benefit after dosing stops?).

One side calculation: early looks at the data after 200 and 400 dogs, with a standard cautious stopping rule (O'Brien-Fleming type). ASSUMPTION: the rule is ours; the published rule could not be read. Reported, not scored.

**Part B: a trial like STAY.** 1,300 dogs, half treated throughout, followed for 4 years, entry at age 10 or older. Published: the number of dogs, the minimum age, about 4 years, nearly all sizes. ASSUMPTIONS: death rates from this project's fit to the UK table; no dropout; and the entry ages, for which two guesses are run: evenly spread from 10 to 14, and shaped like the dog population aged 10 to 16 (fewer dogs at each older age, in step with survival). True effects: none; risk cut 20% and 30%; aging slowed 20% and 30%. Analyses 1 and 2. Also reported, not scored: the size of risk cut that adds one year of life.

**Part C: experiment 010 with the rapamycin trial's dogs.** The same dogs as Part A, but treated throughout and followed for 10 years, with 580 and with 1,300 dogs. True effects: none; a 50% risk cut; aging slowed 20%. Analysis 2.

At least 400 simulated trials per case, fixed seeds.

## Predictions

1. **P1.** In the rapamycin trial as published, a 20% slowing of aging during the dosing year (gain kept) is close to invisible: the any-difference test detects it in fewer than 10% of trials, and the slope test reports a slope change in fewer than 10%.
2. **P2.** In the rapamycin trial as published, even a 50% cut in risk that lasts only for the dosing year is detected by the any-difference test in 25 to 50% of trials. The slope test misreads it as a change in the rate of aging in more than 20% of trials.
3. **P3.** In the STAY-like trial, 60 to 85% of placebo dogs die within 4 years under both guesses about entry ages. A lasting 20% risk cut is detected in more than 85% of trials. A true 20% slowing of aging is recognised by the slope test in fewer than 35% of trials.
4. **P4.** With the rapamycin trial's dogs treated throughout and followed for 10 years, the slope test recognises a 20% slowing in 80 to 95% of 1,300-dog trials and in fewer than 65% of 580-dog trials.

## What would change our mind

If P1 or P3's last part is refuted upward, the running trials can see more of a slope effect than experiments 001 and 003 said, and those statements need correcting. If P4 is refuted, the design recommended in experiment 010 depends on which death-rate inputs are used, and it should be stated as a range across inputs.

These are computational experiments. They show what the mathematics and published inputs imply. They do not replace a real trial, and they say nothing about whether any treatment works.
