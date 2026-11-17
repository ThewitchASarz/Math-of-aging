# Experiment 014: prediction

Written before any code for this experiment was run.

## Question

In a 1,300-dog, 10-year trial, when dogs leave the trial for reasons linked to their health or to the treatment, do the standard analyses report slowed aging or a cut in risk that is not there?

## Why it matters

Experiment 011 let dogs drop out of the trial, but only for reasons that had nothing to do with their health, and equally in both groups. Dropout of that kind costs information and nothing else. Real dropout is often not like that. An owner may withdraw a dog when it falls ill, so its death is never recorded. An owner may stop a treatment that seems to cause side effects. If the deaths that go missing are not a fair sample, an analysis can be wrong, not only weaker. This experiment measures how wrong, for several kinds of unfair dropout.

Nothing here says what dropout in a real dog trial looks like. No published dropout rate or pattern for a dog longevity trial has been found by this project. Every dropout pattern below is an assumption, chosen to see which kinds matter.

## Method, fixed in advance

**The trial.** The "realistic" trial of experiment 011: 1,300 dogs, half treated for the whole trial, entry ages 7 to 12, 10 years of follow-up, a mix of five size classes (this repository's size-class calibration and the US size shares), and a hidden personal risk factor of spread 0.102 (experiment 004's tentative fit). Every dog also has a background chance of leaving of 5 percent a year, unrelated to health. All of this is reused from `experiments/011_realistic_trial/realistic_trial.py`.

**Dropout patterns added on top (all ASSUMPTIONS):**

| Name | What happens |
| --- | --- |
| Reference | Nothing added. This is experiment 011's realistic scenario. |
| Treated dogs leave more | Treated dogs leave at 10 percent a year, untreated at 5, unrelated to health. |
| Sick dogs leave, equally | In both groups, 30 percent of the dogs that are about to die are withdrawn first, at a random time in the 6 months before death. Their deaths are not recorded. |
| Sick treated dogs leave more | The same, but 30 percent in the treated group and 15 percent in the untreated group. |
| Sick untreated dogs leave more | The mirror image: 15 percent treated, 30 percent untreated. |
| Frail treated dogs leave | A treated dog leaves at 5 percent a year plus 5 percent times its hidden risk factor. Untreated dogs leave at 5 percent. |
| Missed deaths drift, to 15 percent | In the treated group only, the share of sick dogs withdrawn before death rises in a straight line from 0 at enrolment to 15 percent at year 10. |
| Missed deaths drift, to 30 percent | The same, rising to 30 percent at year 10. |

**True effects.** No effect (all eight patterns). A plain 50 percent cut in risk, and a 20 percent slowing of aging (the first, second, third, fourth and sixth patterns for slowing; the second, third, fourth and sixth for the risk cut).

**Analyses of each simulated trial**, both with a separate baseline for each of four size groups, as in experiments 007 and 011, and both treating every dog that leaves as simply "no longer watched":
1. The slope test: does the rate at which risk rises differ between the groups? (likelihood-ratio test, 5 percent level). We also record which way the estimate points.
2. The risk test: does the level of risk differ between the groups, with the slope held equal? (same kind of test). We record the estimated risk ratio.

800 simulated trials per case, fixed seeds.

## Predictions

1. **P1. Dropout that is linked to treatment only, or to health equally in both groups, does not fool the slope test.** With "treated dogs leave more" and with "sick dogs leave, equally", the slope test fires in 3 to 8 percent of trials when the treatment does nothing, and in 3 to 8 percent when it plainly halves risk.
2. **P2. Unequal withdrawal of sick dogs creates a benefit that is not there, in the level of risk and not in the slope.** With "sick treated dogs leave more" and a treatment that does nothing, the risk test fires in more than 50 percent of trials, the middle estimate of the risk ratio is between 0.78 and 0.88, and the slope test still fires in only 3 to 8 percent.
3. **P3. The slope test is fooled only when the share of missed deaths drifts over time in one group.** With a treatment that does nothing: "frail treated dogs leave" gives under 9 percent; a drift to 15 percent gives 6 to 15 percent; a drift to 30 percent gives 15 to 40 percent, and in at least 90 percent of the trials where it fires the estimate points to slowed aging.
4. **P4. Unfair dropout of these kinds costs some power but not most of it.** A true 20 percent slowing of aging is recognised in 65 to 82 percent of trials under each of "treated dogs leave more", "sick dogs leave, equally" and "sick treated dogs leave more" (the reference is about 79 to 82 percent in experiment 011).

## What would change our mind

If P1 fails, experiment 011's trial size is not safe even with fair-looking dropout, and its false-alarm figures need a warning. If P2's first part fails, unequal withdrawal is less of a danger to the ordinary "did it cut deaths" result than we think. If P3 fails upward for "frail treated dogs leave", hidden differences between dogs matter more than experiments 004 and 012 found. If the drift cases stay near 5 percent, the slope test is safer than we expect.

These are computational experiments. They show what the mathematics and the stated assumptions imply. They do not replace a real trial, they say nothing about whether any treatment works, and they do not describe the dropout of any real trial.
