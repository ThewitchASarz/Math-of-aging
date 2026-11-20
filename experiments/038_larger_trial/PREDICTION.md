# Experiment 038: prediction

Written before any code for this experiment was run.

## Question

When a treatment helps sturdy and frail dogs unequally but changes no dog's rate of aging, does enrolling more dogs make a slope test report a changed rate of aging more often, and how does that compare with the gain in power for a true slowing?

## Why it matters

A slope test asks whether treated dogs' risk of death rises more slowly, or faster, with age. Experiments 017, 033 and 037 found a trap. If a treatment cuts risk more in sturdy dogs than in frail ones, the test reports slowed aging, though no dog ages more slowly. If it helps frail dogs more, the test reports faster aging.

All of those experiments used 1,300 dogs. The usual answer to a weak trial is a bigger one. That answer may be wrong here. The trap is not noise. It is a steady error in the estimated rate. More dogs shrink the noise around the estimate, not the error in it, so a bigger trial should report the false finding more often and with more confidence.

If so, a large trial that reports slowed aging from death records alone is not safer than a small one, and the size of a trial says little about whether such a finding can be trusted.

This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works, or whether any real treatment helps frail or sturdy dogs more.

## Method, fixed in advance

**The model** is that of experiments 033 and 037. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r * exp(s * b * t)

where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread 0.3 among dogs alive at age 2; the spread is an ASSUMPTION), r is the treatment's cut in risk for that dog (1 in untreated dogs) and s is the change in the rate of aging (1 everywhere except in the "true slowing" cases, where it is 0.8 in treated dogs and r is 1). a and b are the UK fit in `dog_lifespan/results.json`.

**The trials:** 600, 1,300, 2,600 or 5,000 dogs, half treated from enrolment, one size, no dropout.

- 10 years, entry ages 7 to 12 (ASSUMPTION).
- 4 years, entry ages 10 to 14 (ASSUMPTION).

**Benefit patterns (all ASSUMPTIONS).** The sturdier half of dogs (lower hidden risk among dogs of the same entry age) gets one cut and the frailer half another. Written as "sturdier half and frailer half".

- 10 years, at 600, 2,600 and 5,000 dogs: 55 and 45; 60 and 40; 40 and 60 (frail helped more); 42 and 18 (a weaker treatment); and a true 20 percent slowing. The 1,300-dog figures for these are already known. Also an equal cut of 50 percent at 5,000 dogs, as a reference.
- 4 years, at 600, 1,300, 2,600 and 5,000 dogs: 55 and 45; 60 and 40. A true 20 percent slowing at 600, 2,600 and 5,000 dogs.

**Two analyses of every trial**, both reused unchanged from experiment 017: the plain slope test, and the hidden-differences slope test of experiment 012 (one extra setting for the spread of hidden risk). A trial "reports slowed aging" when the test rejects "the rate of aging did not change" at the 5 percent level and the estimated rate is below 1; it "reports faster aging" when the estimate is above 1. With a sound test each happens in about 2.5 percent of trials.

600 simulated trials per case with 600 or 1,300 dogs, 300 per case with 2,600 or 5,000 dogs. Fixed seeds. The script runs in parts of under 5 minutes each.

**Known before this experiment** (1,300 dogs, spread 0.3; experiments 012, 017, 033 and 037). Ten years, hidden-differences test: slowed aging reported in 9.8 percent of trials at 55 and 45, 28.7 at 60 and 40, and about 28 to 32 at 42 and 18; faster aging reported in 20.3 percent at 40 and 60; a true 20 percent slowing recognised in 81 percent (plain test 69). Equal 50 percent cut: the plain test reports faster aging in about 21 percent and the hidden-differences test flags about 5.5 percent in all. Four years: a true slowing is recognised in 26 percent (plain 23); 70 and 30 is reported as slowed aging in 23 percent (plain 12).

**How the ranges below were made.** If the error in the estimated rate stays the same and only the noise shrinks, the distance between the error and the noise grows with the square root of the number of dogs. The ranges are that rule applied to the 1,300-dog figures, with a margin.

## Predictions

1. **P1. In 10-year trials, more dogs mean more false reports of slowed aging.** The hidden-differences test reports slowed aging at 60 and 40 in 8 to 22 percent of trials with 600 dogs, 40 to 62 percent with 2,600 and 65 to 88 percent with 5,000. At 55 and 45: 3 to 11 percent with 600 dogs, 10 to 24 with 2,600 and 17 to 36 with 5,000.
2. **P2. The same holds for the other two patterns.** At 40 and 60 the hidden-differences test reports faster aging in 5 to 18 percent of trials with 600 dogs, 28 to 50 with 2,600 and 48 to 72 with 5,000. At 42 and 18 it reports slowed aging in 9 to 24 percent with 600 dogs, 42 to 65 with 2,600 and 68 to 90 with 5,000.
3. **P3. In 4-year trials the false reports grow about as fast as the power.** At 60 and 40 the hidden-differences test reports slowed aging in fewer than 10 percent of trials with 600 dogs, 4 to 14 percent with 1,300, 8 to 20 with 2,600 and 14 to 32 with 5,000. A true 20 percent slowing is recognised in 8 to 20 percent with 600 dogs, 36 to 58 with 2,600 and 62 to 84 with 5,000.
4. **P4. With an equal cut the hidden-differences test stays sound at 5,000 dogs and the plain test does not.** With a 50 percent cut for every dog in a 10-year trial of 5,000 dogs, the hidden-differences test flags 3 to 8 percent of trials in all, and the plain test reports faster aging in 48 to 75 percent.

Reported, not scored: the plain test in every case, the other direction, the median estimated rate at each size (expected not to move with the number of dogs), the share of dogs that die, 55 and 45 over 4 years, the 10-year power figures, and how well the square-root rule predicts each figure from the 1,300-dog one.

## What would change our mind

If P1 and P2 fail low (the false reports do not grow, or grow much less than the square-root rule says), the error is not a fixed error in the estimate, and a larger trial is not more exposed. If they fail high, the rule understates the problem. If P3 shows the power growing much faster than the false reports, a large 4-year trial would be a better bargain than this prediction says. If P4's first part fails, the hidden-differences test has a small error of its own that only a large trial reveals.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
