# Experiment 044: prediction

Written before any code for this experiment was written or run.

## Question

When a trial runs two tests of "does the benefit differ by frailty score?" side by side, how often does the pair raise a false alarm, does a setting for hidden differences bring each test back to 5 percent, is a bent-line test a useful middle way, and how many dogs does each test need to notice a 60 and 40 split?

## Why it matters

Experiments 032 to 038 found that a slope test can report slowed aging that is not there when a treatment helps sturdy dogs more than frail ones. A frailty score measured at entry does not repair this. Its use is as a warning: a test of "does the benefit differ by score?" tells the trial when a finding of slowed aging should not be relied on.

Experiment 036 found two sound versions of that warning test. The straight-line test is sharper when the benefit rises or falls steadily with the score. The five-group test sees patterns that are not steady. It suggested running both. That left four things unmeasured.

- Two tests are two chances of a false alarm. By how much does the pair exceed 5 percent?
- The check of experiment 036 found each test a little above 5 percent (about 6 to 7) when nothing differs, and guessed the reason: the analyses have no setting for hidden differences between dogs that the score does not see. Is that the reason?
- A bent-line test (the score and its square) has two extra settings, between the straight-line test's one and the five-group test's four. Does it see both kinds of pattern?
- With cuts of 60 and 40 percent and a score of quality 0.7, both tests miss the difference in most 1,300-dog trials. How many dogs would each need?

This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works, whether any treatment helps some dogs more, or whether any real frailty score is good enough.

## Method, fixed in advance

**The model and the trial** are those of experiment 036, unchanged. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r

times exp(b * t), where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread 0.3 among dogs alive at age 2: ASSUMPTION) and r is the treatment's cut in risk for that dog (1 in untreated dogs). No dog's rate of aging changes in any case. a and b are the UK fit in `dog_lifespan/results.json`. Half the dogs are treated from enrolment, entry ages are 7 to 12 (ASSUMPTION), the trial runs 10 years, and there is no dropout. The trial has 1,300 dogs unless stated.

**Frailty scores (ASSUMPTION throughout)**, as in experiments 034 and 036: the rank score (the dog's rank in hidden risk among dogs of its own entry age, read with error) and the risk score (the dog's hidden risk itself, read with error). Quality runs from 0 (useless) to 1 (perfect).

**True effects (sizes are ASSUMPTIONS)**, as in experiment 036: no effect; risk cut by 50 percent in every dog; sturdy dogs helped more (cuts of 70 and 30 percent, and of 60 and 40 percent, in the sturdier and frailer halves); middle dogs helped more (70 percent in the middle half of the hidden-risk ranking, 30 percent in the sturdiest and frailest quarters).

**Tests of "does the benefit differ by score?"** Each compares two fitted models with a likelihood-ratio test at the 5 percent level. No model has a change in the rate of aging. Entry age has its own setting in all of them.

- **Straight-line** (experiment 036's): the cut may change in a straight line with the score. One extra setting.
- **Five groups, both arms** (experiment 036's): a risk level for each fifth of the score in both arms; one cut against a cut per fifth. Four extra settings.
- **Bent-line** (new here): risk depends on the score and its square in both arms; one cut against a cut that changes with the score and its square. Two extra settings.
- **Either of two**: the trial raises the warning if the straight-line test or the five-group test fires. Also reported with each test run at the 2.5 percent level, the textbook correction for running two tests.
- **With a spread setting**: each of the three tests again, with one more setting in both of its models for hidden differences the score leaves over (a gamma-shaped spread among dogs at entry, as in experiment 035).

**Runs.**

- False alarms: no effect and the equal 50 percent cut, with the rank score and the risk score at quality 0.7, 1,300 dogs, at least 800 simulated trials per case.
- Sharpness at 1,300 dogs: 70 and 30, 60 and 40, and middle dogs helped more with the rank score at quality 0.7; middle dogs helped more also at quality 1.0. At least 400 trials per case.
- Number of dogs: 60 and 40 with the rank score at quality 0.7, at 1,300, 2,600, 4,000, 5,500 and 7,000 dogs, tests without the spread setting, at least 400 trials per size. The number of dogs a test needs for 80 percent is read off by drawing a straight line between the two sizes on either side of 80 percent. If a test has not reached 80 percent at 7,000 dogs the answer is "more than 7,000".

Fixed seeds.

## Predictions

1. **P1. Running both tests roughly doubles the false alarms.** In each of the four false-alarm cases the either-of-two rule fires in 8 to 12 percent of trials. With each test run at the 2.5 percent level it fires in 3.5 to 6.5 percent.
2. **P2. A spread setting brings each test back to about 5 percent.** Pooled over the four false-alarm cases, the straight-line test without the spread setting fires in 5.5 to 7.5 percent of trials. With the spread setting the straight-line test and the five-group test each fire in 4 to 6 percent, pooled, and the straight-line figure is lower than without the setting.
3. **P3. The bent-line test is a true middle way.** With the rank score at quality 0.7 and 1,300 dogs: for sturdy dogs helped more (70 and 30) it notices in 74 to 86 percent of trials, between the five-group test and the straight-line test. For middle dogs helped more it notices in 28 to 48 percent, more than the five-group test on the same trials, and in at least 95 percent at quality 1.0.
4. **P4. Noticing a 60 and 40 split takes several thousand dogs.** With the rank score at quality 0.7 the straight-line test needs 3,300 to 5,200 dogs to notice in 80 percent of trials, the bent-line test needs 4,300 to 6,500, and the five-group test needs more than 5,500.

Reported, not scored: every figure with the spread setting other than those in P2; the either-of-two rule's sharpness; each false-alarm case on its own.

## What would change our mind

If P1 fails low (the pair fires in under 8 percent), the two tests mostly fire on the same trials and running both costs little. If it fails high, the pair needs the correction more than expected. If P2 fails because the tests are already at 5 percent without the setting, the check's larger figure in experiment 036 was chance. If it fails because the spread setting does not help, the small excess has another cause and stays an open issue. If P3 fails, the bent-line test is not worth adding. If P4 fails low, a warning test for a 60 and 40 split is within reach of a trial of the size now running; if it fails high, it is further out of reach than thought.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
