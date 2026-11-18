# Experiment 013: prediction

Written before any code for this experiment was run.

## Question

If dog breeds differ in how fast their death risk rises with age, and not only in their level of risk, does the slope test misread a plain cut in risk as a change in the rate of aging, and does giving each common breed its own level and its own rate repair it?

## Why it matters

Experiment 008 found that breed differences in the level of risk do not disturb the slope test. It assumed every breed's risk rises at one shared rate. That is not established. Published work disagrees on whether large dogs die younger because their risk starts higher or because it rises faster (Kraus et al. 2013, Yordy et al. 2020, da Silva & Cross 2023, as recorded in `papers/literature-check.md`). Experiment 008 itself warned that a baseline per breed would not repair differences in the rate. A long trial that plans a slope test needs to know whether it must also plan for this.

"Rate of aging" here means the rate at which death risk rises with age. Nothing here says any treatment works.

## Method, fixed in advance

**Breeds.** The 19 groups (18 breeds plus crossbreds) of Teng et al. 2022, already in `dog_lifespan/data/teng2022_breed_e0.csv`, with the breed lists of experiment 008: without the French Bulldog (the main list, because that breed's figure is unreliable), and all 19 groups.

**How each breed's risk is set (ASSUMPTIONS).** Each breed's risk follows the doubling law with its own level and its own rate, chosen so that the breed's life expectancy equals the published figure. A share f of the difference between breeds is put on the rate:

- A reference level is found so that the UK rate (0.3195 per year, `dog_lifespan/results.json`) gives the average life expectancy of the breed list.
- For each breed, the "all-rate" rate is the one that gives the breed's published life expectancy at that reference level.
- The breed's rate is the UK rate times (all-rate rate / UK rate) to the power f. Its level is then solved so that its life expectancy matches the published figure.
- f = 1: breeds differ only in rate ("all on the rate"). f = 0.5: "half and half". f = 0 is experiment 008.

No breed-by-breed rates were found in a source we could open, so f is an assumption, and f = 1 is meant as the extreme case. The breed mix at birth follows the table's counts, and a dog enters a trial only if alive at its entry age, as in experiment 008. No other hidden differences, no dropout.

**Designs.** 1,300 dogs for 10 years (entry ages 7 to 12), the design of experiments 008 and 010. For the main case only, also 1,300 dogs for 4 years (entry ages 10 to 14). Entry ages are assumptions.

**True effects.** None; a 50 percent cut in risk; the rise of risk slowed by 20 percent in every breed.

**Three analyses of every simulated trial.**

1. Ignoring breed: the slope test of experiment 001.
2. A level per breed, one shared rate: the breed-aware analysis of experiment 008.
3. A level and a rate per breed: every breed expected to supply at least 3 percent of the dogs entering the trial gets its own level and its own rate; all rarer breeds are pooled into one group with its own level and rate. The treatment multiplies the level by one number and the rate by one number, shared across breeds. The slope test asks whether the second number differs from 1.

**Cases.** Main list with f = 1, main list with f = 0.5, all 19 groups with f = 1, each with the three effects in the 10-year design; and the main list with f = 1 in the 4-year design. At least 500 simulated trials per case, fixed seeds.

## Predictions

1. **P1.** With no true effect, breed differences in rate do not create false alarms: every one of the three analyses reports a slope change in 3 to 8 percent of trials, in all four combinations of breed list, f and design.
2. **P2.** With a true 50 percent risk cut in the main case (without the French Bulldog, f = 1, 10 years), ignoring breed misreads it as a slope change in 8 to 20 percent of trials, with the middle estimate of the rate multiplier between 1.02 and 1.10 (pointing to faster aging). A level per breed with a shared rate does not repair it: more than 8 percent. In the half-and-half case, ignoring breed gives under 10 percent.
3. **P3.** A level and a rate per breed repairs it: in every case with a 50 percent risk cut (all four combinations) a slope change is reported in 3 to 8 percent of trials.
4. **P4.** The repair keeps most of the power: with a level and a rate per breed, a true 20 percent slowing is recognised in 75 to 90 percent of 1,300-dog, 10-year trials in the main case, and in 15 to 35 percent of 4-year trials.

## What would change our mind

If P2 is refuted on the low side (under 8 percent), breed differences in rate of this size are not a practical threat and no special analysis is needed. If it is refuted on the high side, the threat is larger than hidden differences of the kind tested in experiment 012. If P3 is refuted, a long trial needs a different repair. If P4 is refuted on the low side, the repair costs so much power that the trial sizes of experiments 010 and 011 would have to grow.

These are computational experiments. They show what the mathematics and the published breed figures imply under stated assumptions. They do not replace a real trial.
