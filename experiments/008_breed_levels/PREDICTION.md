# Experiment 008: prediction

Written before any code for this experiment was run.

## Question

Are differences in risk between dog breeds, across the plausible range (spread 0.15 to 0.71), large enough to make the slope test misread a plain cut in risk as slowed aging, and does giving each breed its own baseline fix it?

## Why it matters

Experiment 004's check found that hidden differences between dogs with a spread of 0.3 break the slope test. Experiment 007 found that published breed life expectancies give a spread between 0.15 and 0.71, depending on which breeds are counted. A long trial that reads a risk cut as "slowed aging" would report the wrong kind of result.

## Method, fixed in advance

- **Breeds.** The 19 groups (18 breeds plus crossbreds) of Teng et al. 2022, already in `dog_lifespan/data/teng2022_breed_e0.csv`. Each breed's level of risk is set from its published life expectancy by the same conversion as experiment 007. Three versions of the breed list give the three spreads measured there: all 19 groups (0.71), without the French Bulldog (0.29), and without the four shortest-lived breeds (0.15).
- **ASSUMPTIONS.** Breeds differ only in their level of risk and share one rate of aging (the UK fit). The mix of breeds at birth follows the table's counts. Dogs enter a trial only if alive at their entry age, so short-lived breeds are rarer among trial dogs than at birth. There are no other hidden differences between dogs.
- **Designs.** 1,300 dogs for 4 years (entry ages 10 to 14) and 1,300 dogs for 10 years (entry ages 7 to 12), as in experiments 001 and 007.
- **True effects.** None; a 50 percent cut in risk; aging slowed by 20 percent.
- **Two analyses of every simulated trial.** Ignoring breed (the slope test of experiment 001), and breed-aware (one baseline level of risk per breed, one shared rate of aging).
- At least 300 simulated trials per case, fixed seeds.

## Predictions

1. **P1.** At the lowest spread (0.15), ignoring breed is still safe: a pure 50 percent risk cut is misread as a slope change in fewer than 10 percent of 1,300-dog, 10-year trials.
2. **P2.** At the higher spreads, ignoring breed breaks the test in the 10-year design: a pure 50 percent risk cut is misread as a slope change in more than 10 percent of trials at spread 0.29, and in more than 15 percent at spread 0.71.
3. **P3.** The breed-aware analysis repairs it: in every case with no true change of slope (no effect, or a 50 percent risk cut; both designs; all three spreads) a slope change is reported in 3 to 8 percent of trials.
4. **P4.** The breed-aware analysis keeps its power: a true 20 percent slowing of aging is recognised in more than 80 percent of 1,300-dog, 10-year trials at every spread.

## What would change our mind

If P2 is refuted, breed differences of the published size do not threaten the slope test and experiment 007's open question closes with "no". If P3 is refuted, a baseline per breed is not enough and a long trial needs a different analysis. If P4 is refuted, the fix costs so much power that the 1,300-dog design of experiment 010 has to grow.

These are computational experiments. They show what the mathematics and the published breed figures imply. They do not replace a real trial.
