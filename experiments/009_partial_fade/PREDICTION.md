# Experiment 009: prediction

Written before any code for this experiment was run.

## Question

In a trial that doses for one year and then watches, what does the two-period test report when the benefit fades slowly or only part of the way, and does a longer or a larger trial read it better?

## Why it matters

Experiment 006 built the two-period test and tried it on two cases only: a benefit that vanishes the day dosing stops, and one that lasts in full. Its independent check warned that a partial fade was not tested. Experiment 019 tested one partial case (half the benefit lasts) in one design (580 dogs, 3 years). Published switch experiments in animals show that part of a benefit can last (Hahn et al. 2019, mice, as recorded in `papers/literature-check.md`). Nobody knows how quickly a benefit would fade in a dog.

The two-period test gives one number for the whole time after dosing. A benefit that is slowly on its way out and a benefit that lasts for good could give the same number. A trial that plans this test needs to know how often each is reported as "lasting", and whether more dogs or more years of watching changes the answer.

Nothing here says any treatment works, lasts or fades. Every effect size and every fading speed is an assumption, chosen to see what a trial could detect.

## Method, fixed in advance

**The model.** Death risk follows the doubling law with the rapamycin trial's published level 0.0066 and slope 0.29 per year. Half the dogs enter at age 7 to 8, a quarter at 8 to 9, a quarter at 9 to 10. These come from `interventions/triad_published_inputs.py` (Coleman et al. 2025, as recorded in this repository). Half the dogs are treated for 12 months. A treated dog's risk is the untreated risk at its age times a multiplier. No dropout, no hidden differences between dogs.

**Three designs.**

| Name | Dogs | Years followed | Status |
| --- | --- | --- | --- |
| D1 | 580 | 3 | The rapamycin trial's published design |
| D2 | 580 | 10 | ASSUMPTION: the same trial watched for longer |
| D3 | 5,000 | 3 | ASSUMPTION: a much larger trial |

**True effects simulated (all sizes and speeds are ASSUMPTIONS).** "u" is years since dosing stopped.

| Name | Multiplier in the dosing year | Multiplier after dosing |
| --- | --- | --- |
| No effect | 1 | 1 |
| Fades at once | 0.5 | 1 |
| Half lasts | 0.5 | 0.75 for good |
| A fifth lasts | 0.5 | 0.9 for good |
| Half of a smaller benefit lasts | 0.7 | 0.85 for good |
| Quick fade (halves every 6 months) | 0.5 | 1 - 0.5 x 2^(-u / 0.5) |
| Slow fade (halves every 2 years) | 0.5 | 1 - 0.5 x 2^(-u / 2) |

In the two fading cases the benefit shrinks smoothly and is never exactly gone. With the quick fade the risk cut is about 3 percent two years after dosing. With the slow fade it is about 25 percent two years after dosing and about 2 percent nine years after.

**Analyses of each simulated trial (reused from earlier experiments).**
1. The two-period test of experiment 006: is the risk ratio after dosing different from 1? "Reported as lasting" means the test fires and the estimate is below 1. We also record the typical (median) after-dosing ratio it estimates.
2. The any-difference test: a two-sided log-rank test at the 5 percent level over the whole trial. Reported, not scored.

1,000 simulated trials per case for the 580-dog designs and 400 for the 5,000-dog design. Fixed seeds.

## Predictions

1. **P1. More years or more dogs both rescue a benefit that half lasts.** In the "half lasts" world the two-period test reports a lasting benefit in 25 to 45 percent of D1 trials, in more than 85 percent of D2 trials and in more than 95 percent of D3 trials.
2. **P2. A small lasting part is missed in every design.** In the "a fifth lasts" world it reports a lasting benefit in fewer than 15 percent of D1 trials and in fewer than 50 percent of both D2 and D3 trials.
3. **P3. A slow fade is read as a lasting benefit, and longer watching does not clear that up.** In the "slow fade" world the test reports a lasting benefit in 45 to 70 percent of D1 trials, with a typical after-dosing ratio of 0.60 to 0.72. In D2 the typical ratio moves up to 0.78 to 0.92 and the test still reports a lasting benefit in 35 to 70 percent of trials.
4. **P4. A large trial reports a benefit as lasting when it is almost gone.** In the "quick fade" world the test reports a lasting benefit in fewer than 20 percent of D1 trials but in 50 to 85 percent of D3 trials.

Not scored, but checked: with no effect, and with a benefit that fades at once, the test should report an after-dosing difference in about 5 percent of trials in every design.

## What would change our mind

If P1 fails for D2, watching longer is not a cure for the weak result of experiment 019, and only more dogs would do. If P3 or P4 fails downward, the two-period test is safer against slow fades than we think, and "lasting" can be taken more at face value. If they hold, a trial should say in advance that "lasting" means "still there on average during the time we watched", and should look at whether the after-dosing ratio drifts over time (queued as experiment 021).

These are computational experiments. They show what the mathematics and the published design imply. They do not replace a real trial, they say nothing about whether any treatment works, and they do not question the trial's own power statement, which answers a different question.
