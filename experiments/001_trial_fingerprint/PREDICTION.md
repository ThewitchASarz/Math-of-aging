# Experiment 001: prediction

Written before any code for this experiment was run.

## Question

A treatment can lower the risk of death in two ways. It can cut risk by a fixed fraction, which shifts the mortality curve sideways and delays death (a *level* effect). Or it can slow the yearly growth of risk, which changes the slope of the curve (a *rate* effect). Only the second is slowed aging (Theorems 2 and 5 of the theory paper).

Can the dog trials now running tell these two apart?

## Method, fixed in advance

- Simulate two-arm trials with dog mortality from this repository's fit to the UK life table.
- Designs: a 580-dog trial with 3 years of follow-up (dogs enrolled at ages 7 to 12) and a 1,300-dog trial with 4 years of follow-up (ages 10 to 14). These mirror the published sizes of the Dog Aging Project rapamycin trial and Loyal's STAY trial. Entry-age ranges are our assumption. Two hypothetical designs for comparison: 5,000 dogs for 3 years, and 1,300 dogs followed for 10 years.
- True effects: level cuts of 20, 30 and 50 percent; rate slowdowns of 10, 20 and 30 percent; and no effect.
- Analysis: fit by maximum likelihood, with each dog's own entry age. The fingerprint test is a likelihood-ratio test of "the slope did not change" against "the slope changed", at the 5 percent level.
- At least 200 simulated trials per case. Fixed seed.

## Predictions

1. **P1.** In the 580-dog, 3-year design, a true 20 percent slowing of aging is recognised as a slope change in fewer than 30 percent of trials.
2. **P2.** In the 1,300-dog, 4-year design, the same effect is recognised in fewer than 50 percent of trials.
3. **P3.** Follow-up matters more than headcount: 1,300 dogs followed for 10 years recognise it more often than 5,000 dogs followed for 3 years.
4. **P4.** The test does not cry wolf. When the truth is a pure level cut, it reports a slope change in 3 to 8 percent of trials.

## What would change our mind

If P1 and P2 are refuted, the running trials can already separate the two kinds of effect, and Prediction 2 of the theory paper is testable with their data as planned.
