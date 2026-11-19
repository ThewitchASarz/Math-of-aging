# Experiment 037: prediction

Written before any code for this experiment was run.

## Question

When a treatment helps frail dogs a little more than sturdy ones, or helps sturdy dogs more but cuts risk by only 20 or 30 percent on average, how often does a slope test still report a change in the rate of aging that is not there?

## Why it matters

A slope test asks whether treated dogs' risk of death rises more slowly, or faster, with age. Experiments 017 and 033 found a trap. If a treatment cuts risk more in sturdy dogs than in frail ones, the test reports slowed aging, though no dog ages more slowly. Experiment 033 found that a modest difference is enough: cuts of 60 and 40 percent mislead the hidden-differences test in about 29 percent of 10-year trials.

Two things were left open.

- **The opposite direction.** If frail dogs are helped more, the tests report faster aging. A real trial that saw this might conclude the treatment speeds up aging and should be dropped. Experiment 017 measured it only for a large difference (30 and 70 percent). Nobody here has measured small differences.
- **A weaker treatment.** Every case so far cut risk by 50 percent on average. That is a very strong treatment. A cut of 20 or 30 percent is closer to what the running dog trials are designed around. If the trap shrinks with the size of the cut, it is a smaller worry for a realistic treatment.

This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works, or whether any real treatment helps frail or sturdy dogs more.

## Method, fixed in advance

**The model** is experiment 033's, with its code reused unchanged. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r * exp(b * t)

where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread 0.3 among dogs alive at age 2; the spread is an ASSUMPTION) and r is the treatment's cut in risk for that dog (1 in untreated dogs). No dog's rate of aging changes in any case. a and b are the UK fit in `dog_lifespan/results.json`.

**The trial:** 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years, one size, no dropout.

**Benefit patterns (all ASSUMPTIONS).** Two halves: the sturdier half of dogs (lower hidden risk among dogs of the same entry age) gets one cut and the frailer half another. Written below as "sturdier half and frailer half".

- Frail dogs helped more, average cut 50 percent: 45 and 55, 40 and 60, 35 and 65, and 30 and 70 (a repeat of experiment 017).
- Average cut 30 percent: 30 and 30 (no difference); 36 and 24, and 42 and 18 (sturdy helped more, the same ratios between the cuts as 60 and 40 and as 70 and 30); 18 and 42 (frail helped more).
- Average cut 20 percent: 20 and 20; 24 and 16, and 28 and 12; 12 and 28.

**Two analyses of every trial**, both reused unchanged from experiment 017: the plain slope test, and the hidden-differences slope test of experiment 012 (one extra setting for the spread of hidden risk). A trial "reports faster aging" when the test rejects "the rate of aging did not change" at the 5 percent level and the estimated rate is above 1; it "reports slowed aging" when the estimate is below 1. With a sound test each happens in about 2.5 percent of trials.

600 simulated trials per case, fixed seeds.

**Known before this experiment** (experiments 017 and 033, same trial). Sturdy helped more, average cut 50 percent: the hidden-differences test reports slowed aging in 9.8 percent of trials at 55 and 45, 28.7 at 60 and 40, 63.7 at 65 and 35 and 86.8 at 70 and 30. Frail helped more at 30 and 70: the plain test reports faster aging in 78 percent and the hidden-differences test in 30 percent. With an equal 50 percent cut the plain test reports faster aging in about 21 percent of trials.

## Predictions

1. **P1. With frail dogs helped more, the hidden-differences test is misled less often than in the mirror case.** It reports faster aging in fewer than 10 percent of trials at 45 and 55, in 5 to 18 percent at 40 and 60, and in 10 to 28 percent at 35 and 65. At 40 and 60 and at 35 and 65 this is less than half the mirror figure from experiment 033 (28.7 and 63.7 percent reporting slowed aging).
2. **P2. The plain test is misled often, even by a small difference.** With frail dogs helped more it reports faster aging in 22 to 42 percent of trials at 45 and 55, in 35 to 58 percent at 40 and 60, and in 48 to 72 percent at 35 and 65.
3. **P3. A weaker treatment makes the trap much smaller.** With sturdy dogs helped more, the hidden-differences test reports slowed aging in 4 to 14 percent of trials at 36 and 24 and in 14 to 32 percent at 42 and 18 (average cut 30 percent), and in fewer than 9 percent at 24 and 16 and 5 to 16 percent at 28 and 12 (average cut 20 percent).
4. **P4. The plain test's lean towards "faster aging" also shrinks with the cut.** With an equal cut for every dog it reports faster aging in 5 to 14 percent of trials at 30 percent and in 3 to 10 percent at 20 percent. With frail dogs helped more at 18 and 42 it reports faster aging in 12 to 32 percent.

Reported, not scored: the other test and the other direction in every case, the median estimated rate, the share of dogs that die, and the repeated reference case (30 and 70).

## What would change our mind

If P1 fails high, a treatment that helps frail dogs a little more would often be reported as speeding up aging even by the better test. If P2 fails low, the plain test is safer in this direction than experiment 033 suggested. If P3 fails high, the trap does not depend on a very strong treatment and matters for a treatment of realistic strength. If P3 fails low (near 2.5 percent everywhere), the trap is mostly a feature of the 50 percent cut used so far.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
