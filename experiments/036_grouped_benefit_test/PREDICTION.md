# Experiment 036: prediction

Written before any code for this experiment was written or run.

## Question

When a trial asks "does the benefit differ between sturdy and frail dogs?" by splitting a frailty score into five groups in both arms, is the test sound, and does it notice a difference as often as the straight-line test?

## Why it matters

Experiments 032 to 035 found that a frailty score measured at entry does not repair the slope test. Its main use is as a warning: a test of "does the benefit differ by score?" tells the trial when a finding of slowed aging should not be relied on. So the warning test itself has to be sound.

Two versions exist so far.

- **Straight-line.** The benefit is allowed to change in a straight line with the score. It noticed sturdy dogs helped more in 88 percent of simulated trials at score quality 0.7 (experiment 034).
- **Five groups, treated dogs only.** A separate risk cut for each fifth of the score. Experiment 034 found it faulty: with a treatment that does nothing and a perfect rank score it fired in 15 percent of trials, where 5 is expected. The check traced this to the grouping being applied to treated dogs only, and got about 6 percent when both arms were grouped. That figure was the check's own and was never saved or tested properly.

A straight-line test can only see a benefit that rises or falls steadily with the score. A grouped test can see any pattern, but spends more settings to do so. This experiment tests the both-arms grouped version properly, and asks what each test sees when the pattern is not a steady one.

This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works, whether any treatment helps some dogs more, or whether any real frailty score is good enough.

## Method, fixed in advance

**The model and the trial** are those of experiments 032, 034 and 035, unchanged. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r

times exp(b * t), where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread 0.3 among dogs alive at age 2: ASSUMPTION) and r is the treatment's cut in risk for that dog (1 in untreated dogs). No dog's rate of aging changes in any case. a and b are the UK fit in `dog_lifespan/results.json`. The trial has 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

**Frailty scores (ASSUMPTION throughout)**, as in experiment 034: the rank score (the dog's rank in hidden risk among dogs of its own entry age, read with error) at quality 0.5, 0.7, 0.9 and 1.0, and the risk score (the dog's hidden risk itself, read with error) at quality 0.7 and 1.0.

**True effects (sizes are ASSUMPTIONS).**

- No effect.
- Risk cut by 50 percent in every dog.
- Sturdy dogs helped more: risk cut by 70 percent in the sturdier half and 30 percent in the frailer half.
- Sturdy dogs helped a little more: cuts of 60 and 40 percent.
- Middle dogs helped more: risk cut by 70 percent in the middle half of the hidden-risk ranking and 30 percent in the sturdiest quarter and the frailest quarter. This is a pattern with no steady trend.

**Tests of "does the benefit differ by score?"** Each compares two fitted models with a 5 percent likelihood-ratio test. Neither model has a change in the rate of aging. Entry age has its own setting in all of them, as experiment 034 recommends.

- **Straight-line** (experiment 034's): risk depends on the score in a straight line; the benefit is one cut for all treated dogs, against a cut that changes in a straight line with the score. One extra setting.
- **Five groups, treated only** (experiment 034's faulty version, kept for comparison): risk depends on the score in a straight line in both arms; one cut, against a cut per fifth of the score. Four extra settings.
- **Five groups, both arms** (the version under test): risk depends on the score in a straight line plus a level for each fifth of the score, in both arms; one cut, against a cut per fifth. Four extra settings.
- **Three groups, both arms**: the same with thirds. Two extra settings. Reported, not scored.

Fifths and thirds are cut on the scores of all dogs in the trial. At least 600 simulated trials per case, fixed seeds.

## Predictions

1. **P1. Grouping both arms repairs the false alarms.** With no effect, the five-groups-both-arms test fires in 3 to 8 percent of trials for every one of the six scores. The treated-only version fires in 11 to 20 percent with the rank score at quality 1.0.
2. **P2. An equal benefit is not read as an unequal one.** With risk cut by 50 percent in every dog, both the straight-line test and the five-groups-both-arms test fire in 3 to 9 percent of trials with the rank score at quality 0.7 and with the risk score at quality 0.7.
3. **P3. For a steady pattern the grouped test is the weaker warning.** With sturdy dogs helped more (70 and 30), the five-groups-both-arms test notices the difference in 22 to 40 percent of trials with the rank score at quality 0.5, 62 to 80 percent at quality 0.7 and at least 97 percent at quality 1.0. At quality 0.7 that is 6 to 24 points below the straight-line test on the same trials.
4. **P4. For a pattern with no steady trend only the grouped test sees it.** With middle dogs helped more, the straight-line test fires in fewer than 12 percent of trials with the rank score at quality 0.7 and at quality 1.0. The five-groups-both-arms test notices the difference in at least 95 percent at quality 1.0 and in 20 to 45 percent at quality 0.7.

Reported, not scored: the 60 and 40 case; the three-group test; every other score and quality.

## What would change our mind

If P1 fails for the both-arms test, grouping both arms is not the repair the check of experiment 034 suggested, and a grouped warning test should not be used until the cause is found. If P2 fails, hidden differences the score does not see make an equal benefit look unequal, and the warning test cries wolf. If P3 fails high (the grouped test is within a few points of the straight-line test), the grouped test costs little and could be the default. If P4 fails because the straight-line test also sees the middle pattern, or because the grouped test misses it even with a perfect score, the case for running both tests falls away.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
