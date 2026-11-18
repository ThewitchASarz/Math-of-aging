# Experiment 017: prediction

Written before any code for this experiment was run.

## Question

Does the slope test with a setting for hidden differences between dogs (experiment 012) still work in a 4-year trial, in a 10-year trial with dropout and a mix of dog sizes, and when the treatment helps frail dogs more than sturdy ones?

## Why it matters

Experiment 012 found a repair for the one thing that breaks the slope test. When dogs differ widely in hidden risk, a plain cut in risk can be misread as a change in the rate of aging. Adding one setting for the spread of hidden risk brought false alarms back to about 5 percent and cost no power. But that was tested in easy conditions only: a 10-year trial in which nearly every dog dies, no dropout, one size of dog, and a treatment that helps every dog by the same fraction. Experiment 012's own limits (2, 3 and 4) say the result should not be carried to short trials, to trials with dropout and sizes, or to treatments that help frail dogs more, without testing. A trial that writes this analysis into its plan needs to know where it holds.

The idea of a hidden risk factor is long published (see `papers/literature-check.md`). This experiment only measures how one analysis behaves on simulated dogs. It says nothing about whether any treatment works.

## Method, fixed in advance

**The model.** A dog of size class c that enters at age x0 has, t years later, a death risk of

    z * a * exp(b_c * x0) * r * exp(b_c * s * t)

where z is the dog's hidden personal risk factor, r is the treatment's cut in risk (1 = none) and s is the treatment's change in the rate of aging (1 = none). Among dogs alive at age 2, z has average 1 and a gamma shape with spread (variance) 0.3. This spread is an ASSUMPTION, the wide case of experiments 004, 008 and 012; nobody has measured it. A dog enters a trial only if alive at its entry age, so frail dogs are rarer among older entrants.

**Three parts. Every trial has 1,300 dogs, half treated from enrolment.**

- **Part 1: a 4-year trial.** Entry ages 10 to 14 (ASSUMPTION, as in experiments 001 and 016), 4 years of follow-up, one size, no dropout, death rates from the UK fit in `dog_lifespan/results.json`. True effects: none; a 50 percent risk cut; aging slowed 20 percent. One more case with no hidden differences and aging slowed 20 percent, to measure what the extra setting costs when it is not needed.
- **Part 2: a 10-year trial with dropout and a mix of sizes.** Entry ages 7 to 12 (ASSUMPTION), five size classes with their own rates of aging (this repository's size-class calibration, as in experiments 007 and 011), size mix as in the US dog population (ASSUMPTION), 5 percent of dogs leaving each year for reasons unrelated to health (ASSUMPTION). The same spread of 0.3 inside every size class (ASSUMPTION). Same three true effects. Both analyses give each of four size groups its own level and rate.
- **Part 3: a treatment that helps frail dogs more.** The design of experiment 012 (10 years, entry ages 7 to 12, one size, no dropout, UK fit). No dog's rate of aging changes. Among treated dogs, the frailer half (hidden risk above the middle value for dogs of the same entry age) has its risk cut by 70 percent and the sturdier half by 30 percent (ASSUMPTION). A mirror case with the two cuts swapped. And the frail-helped-more case again in the 4-year trial of Part 1.

**Two analyses of every simulated trial.** The plain slope test of experiment 001 (with size groups in Part 2, as in experiment 007). And the test of experiment 012: the same, with one extra setting for the spread of a gamma-shaped hidden risk factor, estimated from the trial itself. A test "flags" a trial when it rejects "the rate of aging did not change" at the 5 percent level.

At least 500 simulated trials per case, fixed seeds.

## Predictions

1. **P1. The repair still works in a 4-year trial.** In Part 1 with a 50 percent risk cut, the plain test flags 6.5 to 15 percent of trials and the new test 3 to 8 percent. With no effect the new test flags 3 to 8 percent.
2. **P2. The repair still costs no power in a 4-year trial, but there is little power to keep.** With aging truly slowed 20 percent, the new test flags 15 to 32 percent of trials when the spread is 0.3. With no hidden differences, the new test's figure is no more than 5 points below the plain test's.
3. **P3. The repair works with dropout and a mix of sizes.** In Part 2 with a 50 percent risk cut, the plain test flags 12 to 25 percent of trials and the new test 3 to 8 percent. With no effect the new test flags 3 to 8 percent. With aging truly slowed 20 percent the new test flags 60 to 80 percent.
4. **P4. A treatment that helps frail dogs more defeats both tests.** In Part 3's 10-year frail-helped-more case, the plain test flags 25 to 60 percent of trials and the new test more than 10 percent, although no dog's rate of aging changed. The error does not point the hopeful way: fewer than 3 percent of trials are flagged as slowed aging by either test.

Reported, not scored: the mirror case, the 4-year frail-helped-more case, the estimated spread, and how often the estimated spread sits at its lower limit.

## What would change our mind

If P1 or P3 fails because the new test flags far more than 8 percent, the repair of experiment 012 is limited to long, clean trials, and the open threat from experiments 004 and 008 stands for any realistic design. If P2 fails because power drops, the repair has a price in short trials. If P4 fails because the new test stays near 5 percent, the repair is sturdier than we think and limit 2 of experiment 012 can be dropped. If P4 fails the other way (many trials flagged as slowed aging), that is a new way for a trial to report slowed aging that is not there, and it would go at the top of `KNOWN_GAPS.md`.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
