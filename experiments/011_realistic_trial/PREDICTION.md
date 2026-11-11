# Experiment 011: prediction

Written before any code for this experiment was run.

## Question

How many dogs does the 10-year trial of experiment 010 need to recognise a 20 percent slowing of aging once dogs drop out, risk bends at old ages, and the dogs are a mix of sizes?

## Why it matters

Experiment 010 gave a design: about 1,300 dogs, treated throughout, followed for 10 years. It assumed that no dog leaves the trial, that risk follows the doubling law at every age, and that all dogs are alike. A real trial has none of these. If the number moves a lot, the design has to be bigger.

## Method, fixed in advance

- **Size mix.** Dogs belong to the five size classes of experiment 007, each with its own rate of aging (this repository's calibration to published life expectancies). ASSUMPTION: the mix among trial dogs equals the mix in the US population table.
- **Bending at old ages.** Each dog carries a hidden personal risk factor with the spread fitted to the UK table in experiment 004 (0.102). The frailest dogs die first, so the risk of the survivors grows more slowly at old ages. ASSUMPTION: this spread applies inside every size class. Some of the bend in the UK table is itself due to the size mix, so this probably counts part of it twice.
- **Dropout.** Dogs leave the trial at a constant yearly rate, the same in both groups and unrelated to health. Their deaths are not seen. ASSUMPTION: 0, 5 and 10 percent a year. We tried to open the design paper of the rapamycin dog trial for a planned figure and could not, so no published dropout rate is used.
- **Two named scenarios.** "Realistic": size mix, bending and 5 percent dropout a year. "Harsh": the same with 10 percent a year. The names are ours; neither is a measurement.
- **Trial.** Dogs enter at ages 7 to 12, half are treated from enrolment to death, follow-up is 10 years. Sizes: 1,300, 1,600, 2,000 and 2,600 dogs.
- **True effects.** Aging slowed by 20 percent; no effect; a 50 percent cut in risk.
- **Analysis.** The slope test with a separate baseline for each of four size groups, as recommended by experiment 007. The analysis knows nothing about the hidden risk factor or the reason for dropout.
- At least 300 simulated trials per case, fixed seeds. "Enough" means the slowing is recognised in 80 percent of trials.

## Predictions

1. **P1.** Size mix and bending alone cost little. With no dropout, 1,300 dogs recognise the slowing in 80 to 95 percent of trials.
2. **P2.** Dropout costs more. With the size mix and 10 percent dropout a year (no bending), 1,300 dogs recognise the slowing in 60 to 78 percent of trials.
3. **P3.** In the realistic scenario, 1,300 dogs reach 70 to 82 percent and 1,600 dogs are enough. In the harsh scenario, the smallest trial on the grid that is enough has 2,000 dogs.
4. **P4.** The test stays honest. In both scenarios, with 1,300 dogs, a slope change is reported in 3 to 9 percent of trials when there is no effect and when the true effect is a plain 50 percent cut in risk.

## What would change our mind

If P3 is refuted upward (2,600 dogs are not enough in the harsh scenario), the design of experiment 010 is too small by a factor of two or more and its headline needs a warning. If P4 is refuted, dropout or bending makes the test report slowed aging that is not there, and a long trial needs a different analysis.

These are computational experiments. They show what the mathematics and the published data imply under stated assumptions. They do not replace a real trial.
