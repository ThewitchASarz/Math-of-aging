# Experiment 033: prediction

Written before any code for this experiment was run.

## Question

How large must the difference in benefit between sturdy and frail dogs be before a slope test reports slowed aging that is not there in more than 10 percent of trials?

## Why it matters

A slope test asks whether treated dogs' risk of death rises more slowly with age. Experiment 017 found a trap. If a treatment cuts risk more in sturdy dogs than in frail ones, the test reports slowed aging, though no dog ages more slowly. The reason: among treated dogs the frail ones die early and the well-protected sturdy ones remain, so the treated group's risk climbs slowly.

Experiment 017 used one large difference: risk cut by 70 percent in the sturdier half and 30 percent in the frailer half. Experiments 032 and 034 found that a frailty score measured at entry repairs this only if the score is nearly perfect. So the next thing a trial needs to know is how small a difference is already enough to cause trouble. If only large differences matter, the trap is a rare worry. If small ones matter, a finding of slowed aging needs support from outside the death records.

This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works, or whether any real treatment helps sturdy dogs more.

## Method, fixed in advance

**The model** is experiment 017's. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r * exp(b * t)

where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread V among dogs alive at age 2) and r is the treatment's cut in risk for that dog (1 in untreated dogs). No dog's rate of aging changes in any case. a and b are the UK fit in `dog_lifespan/results.json`.

**The trial**, unless stated: 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years, one size, no dropout.

**Benefit patterns (all ASSUMPTIONS).** The average cut is 50 percent in every case. "Sturdier half" means the half with lower hidden risk among dogs of the same entry age.

- Two halves: cuts of 50 and 50 (no difference), 55 and 45, 60 and 40, 65 and 35, 70 and 30 percent, the larger cut going to the sturdier half.
- Steady: the cut falls in a straight line with the dog's rank in hidden risk, from 90 percent in the sturdiest dog to 10 percent in the frailest (same half averages as 70 and 30), and from 70 to 30 percent (same half averages as 60 and 40).

**Cases.**

- Spread 0.3 (ASSUMPTION, as in experiments 012 and 017): all five two-halves patterns and both steady patterns.
- Spread 0.1 (about the size fitted to the UK life table in experiment 004): 60 and 40, and 70 and 30.
- 4-year trial, entry ages 10 to 14 (ASSUMPTION), spread 0.3: 70 and 30.
- 10-year trial with five size classes and 5 percent dropout a year (experiment 017's Part 2; size mix and dropout are ASSUMPTIONS), spread 0.3: 70 and 30.

**Two analyses of every trial**, both reused unchanged from experiment 017: the plain slope test, and the hidden-differences slope test of experiment 012 (one extra setting for the spread of hidden risk). A trial "reports slowed aging" when the test rejects "the rate of aging did not change" at the 5 percent level and the estimated rate is below 1. With a sound test this happens in about 2.5 percent of trials.

600 simulated trials per case, fixed seeds.

**Known before this experiment** (experiments 012 and 017, same trial, spread 0.3): with cuts of 50 and 50 the plain test reports slowed aging in about 0 percent of trials (it leans towards faster aging, flagging 19 percent) and the hidden-differences test in about 3 percent. With 70 and 30 the plain test reports slowed aging in 58 percent and the hidden-differences test in 88 percent.

## Predictions

1. **P1. A modest difference is enough for the hidden-differences test; the plain test needs a larger one.** At spread 0.3 the hidden-differences test reports slowed aging in 6 to 18 percent of trials with cuts of 55 and 45, and in 20 to 45 percent with 60 and 40. The plain test reports it in fewer than 12 percent with 60 and 40, and in 15 to 35 percent with 65 and 35.
2. **P2. A smaller spread of hidden risk weakens the trap but does not remove it.** At spread 0.1 with cuts of 70 and 30, the plain test reports slowed aging in 15 to 45 percent of trials and the hidden-differences test in 30 to 65 percent. With 60 and 40, both are below 25 percent.
3. **P3. A benefit that varies steadily is no safer than two halves.** At spread 0.3, the hidden-differences test reports slowed aging in at least 85 percent of trials with the steady 90-to-10 pattern. With the steady 70-to-30 pattern it reports it in 25 to 55 percent, and no more than 5 points less often than with two halves of 60 and 40.
4. **P4. The trap is weak in a 4-year trial and survives sizes and dropout in a 10-year one.** With cuts of 70 and 30 at spread 0.3: in the 4-year trial the plain test reports slowed aging in fewer than 8 percent of trials and the hidden-differences test in 5 to 20 percent. In the 10-year trial with sizes and dropout the plain test reports it in 35 to 60 percent and the hidden-differences test in 65 to 90 percent.

Reported, not scored: how often each test flags faster aging, the median estimated rate, the share of dogs that die, and the two repeated reference cases (50 and 50, 70 and 30).

## What would change our mind

If P1 fails low (few false reports at 60 and 40), the trap needs a large difference in benefit and is a smaller worry than experiments 017 to 034 suggest. If it fails high, even a 55-against-45 difference misleads often. If P2 fails low, the trap depends on a wide spread of hidden risk, which nobody has measured in dogs. If P3 fails, the two-halves pattern used so far overstated the problem. If P4's first part fails high, the trap also matters for trials of the length now running.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
