# Experiment 004: prediction

Written before any code for this experiment was run.

## Question

Every theorem and experiment so far assumes that a dog's risk of death keeps doubling on a fixed schedule. Does that hold for old dogs, or does risk level off? And if it levels off because some dogs are frailer than others, does that break the fingerprint test of experiment 001?

## Method, fixed in advance

- **Data.** The UK dog life table (Teng and others 2022), which publishes deaths and dogs alive at each age from 2 to 19. The US table publishes rates only and is used for a simpler curvature check.
- **Models.** The plain doubling law, and the same law with hidden frailty (each dog's risk is multiplied by a personal factor drawn from a gamma distribution), which makes population risk level off at old ages. Both fitted by maximum likelihood to the UK counts.
- **Out-of-range test.** Fit the plain law to ages 2 to 13 only, then compare its forecast with observed death rates at ages 14 to 19.
- **Effect on experiment 001.** Simulate trials in which dogs have hidden frailty of the fitted size, analyse them with the plain model as before, and measure how often a pure risk cut is misread as a slope change.

## Predictions

1. **P1.** On the UK counts, the model with levelling-off fits significantly better than the plain doubling law (likelihood-ratio test, p below 0.01).
2. **P2.** Below age 14 the plain law is adequate: fitted on ages 2 to 13, it forecasts death rates at ages 14 to 16 within 15 percent.
3. **P3.** The same fit overstates the death rate at age 18 by more than 20 percent.
4. **P4.** With hidden frailty of the fitted size, a pure 30 percent risk cut is misread as a slope change in more than 10 percent of 1,300-dog, 10-year trials, against the 5 percent the test is designed for.

## What would change our mind

If P4 is refuted, the "follow dogs to death" recommendation from experiment 001 survives hidden frailty. If P4 is confirmed, long trials need an analysis that models frailty, and the recommendation must say so.

## A caveat known in advance

The UK table is built from recorded deaths in a single period, not by following dogs through life. Levelling-off at the oldest ages could partly reflect that construction, or errors in recorded ages, and not biology.
