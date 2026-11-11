# Experiment 012: prediction

Written before any code for this experiment was run.

## Question

If the slope test is given an extra setting for hidden differences between dogs, does it stop misreading a plain cut in risk as slowed aging, and how much of its power to recognise real slowing does it keep?

## Why it matters

Experiments 004 and 008 found one thing that breaks the slope test. When dogs differ widely in hidden risk (spread 0.3), a plain 50 percent cut in risk is reported as a change in the rate of aging in about 22 percent of 10-year trials, against the 5 percent the test is built for. Nobody knows how large hidden differences inside breeds really are. A long trial needs an analysis that is safe either way. The idea of a hidden risk factor that bends the mortality curve is long published (see `papers/literature-check.md`); this experiment only measures what adding it to this test does.

## Method, fixed in advance

- **Trial.** 1,300 dogs, entry ages 7 to 12 (ASSUMPTION), half treated from enrolment, followed for 10 years. Mortality from the UK fit in `dog_lifespan/results.json`. No dropout.
- **Four worlds.** (A) No hidden differences. (B) Hidden differences with spread 0.1, about the size fitted to the UK table in experiment 004. (C) Spread 0.3, the control of experiment 008. In B and C the hidden risk factor follows a gamma distribution among dogs alive at age 2, simulated with experiment 004's function. (D) Spread 0.3 with a different shape, a log-normal distribution, which the analysis does not assume. The spreads in C and D are ASSUMPTIONS: no measurement exists.
- **True effects.** No effect; a plain 50 percent cut in risk; aging slowed by 20 percent.
- **Two analyses of every simulated trial.** The plain slope test of experiment 001. And the new test: the same comparison, with one extra setting for the spread of a gamma-shaped hidden risk factor, estimated from the trial itself, and allowing for the fact that frail dogs are already rarer among older dogs at entry.
- At least 300 simulated trials per case, fixed seeds.

## Predictions

1. **P1.** The problem is reproduced and the new test repairs it. In world C with a 50 percent risk cut, the plain test reports a slope change in more than 15 percent of trials and the new test in 3 to 8 percent.
2. **P2.** The new test is honest everywhere its assumed shape is right. In worlds A, B and C, with no effect or a 50 percent risk cut, it reports a slope change in 3 to 8 percent of trials in all six cases.
3. **P3.** The repair costs power. In world A, where the extra setting is not needed, the plain test recognises a 20 percent slowing in more than 85 percent of trials and the new test in 60 to 85 percent, at least 5 points lower. In worlds B and C the new test still recognises it in at least 60 percent of trials.
4. **P4.** The repair does not depend on guessing the shape exactly. In world D with a 50 percent risk cut, the plain test reports a slope change in more than 15 percent of trials and the new test in fewer than 10 percent.

## What would change our mind

If P1 or P2 is refuted, adding a hidden-difference setting is not a fix and the open threat from experiment 008 stands. If the new test's power falls below 60 percent, the fix is too expensive for a 1,300-dog trial and the design of experiments 010 and 011 has to grow or choose between safety and power. If P4 is refuted, the fix only works when the analyst guesses the shape of the hidden differences, which no one can check.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial.
