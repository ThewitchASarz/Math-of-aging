# Experiment 007: prediction

Written before any code for this experiment was run.

## Question

Experiment 004 found that the slope test survives hidden differences between dogs of the size fitted to the UK table, and its check found that differences three times larger would break it. Dogs differ a great deal by breed and size. Are real differences big enough to break the test, and does modelling size fix it?

## Method, fixed in advance

- **Part 1: how different are dogs?** Convert published life expectancies by breed (UK, 18 breeds plus crossbreds) and by size class (US, five classes) into risk multipliers, using the doubling-law relation between life expectancy and level of risk. Report the spread as a variance comparable with the frailty variance of experiment 004 (0.10 fitted; 0.3 broke the test).
- **Part 2: trials with a mix of sizes.** Simulate trials in which each dog belongs to a size class with its own rate of aging, using this repository's size-class calibration (larger dogs age faster) and the US population shares of the classes. Analyse them two ways: ignoring size, and with a separate baseline for each size class.
- Designs: 1,300 dogs for 4 years (entry ages 10 to 14) and 1,300 dogs for 10 years (entry ages 7 to 12). Effects: none, risk cuts of 30 and 50 percent, and aging slowed by 20 percent. At least 300 simulated trials per case.

The size-class rates are a calibration to published life expectancies, not a fit to full life tables, and the class shares among trial dogs are assumed equal to population shares.

## Predictions

1. **P1.** The spread of risk between breeds corresponds to a variance above 0.2, more than twice the fitted frailty variance.
2. **P2.** Ignoring size, a pure 50 percent risk cut is misread as a slope change in more than 10 percent of 1,300-dog, 10-year trials.
3. **P3.** With a separate baseline for each size class, false alarms return to 3 to 8 percent in every case.
4. **P4.** The size-aware analysis keeps its power: a true 20 percent slowing of aging is recognised in more than 80 percent of 1,300-dog, 10-year trials.

## What would change our mind

If P2 is refuted, a mix of sizes does not threaten the slope test and no adjustment is needed. If P3 is refuted, adjusting for size class is not enough.
