# Experiment 032: prediction

Written before any code for this experiment was run.

## Question

If each dog has a frailty score measured at entry that tracks its hidden risk imperfectly, how good must the score be before adding "benefit differs by score" to the slope test stops a treatment that helps sturdy dogs more from being read as slowed aging?

## Why it matters

Experiment 017 found a way for a trial to report slowed aging that is not there. If a treatment cuts risk more in sturdy dogs than in frail ones, the frail treated dogs die first and the treated group that remains looks better and better against the untreated group. Both slope tests read that as a slower rise of risk with age: the plain test in 58 percent of 10-year trials and the hidden-differences test in 88 percent, with no dog aging more slowly. Experiment 017 said a finding of slowed aging therefore needs a check of whether the benefit differs between frail and sturdy dogs. A real trial cannot see hidden risk. It can only measure something at entry (a frailty index, a vet's score, blood values) that tracks it partly. A trial planning such a check needs to know how good that measurement has to be.

This experiment only measures how one analysis behaves on simulated dogs. It says nothing about whether any treatment works or whether any real frailty score is good enough.

## Method, fixed in advance

**The model** is experiment 017's Part 3. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r * exp(b * s * t)

where z is the dog's hidden personal risk factor, r is the treatment's cut in risk (1 = none) and s is the treatment's change in the rate of aging (1 = none). Among dogs alive at age 2, z has average 1 and a gamma shape with spread (variance) 0.3 (ASSUMPTION, as in experiments 012 and 017; nobody has measured it). a and b are the UK fit in `dog_lifespan/results.json`.

**The trial** is experiment 012's: 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years of follow-up, one size, no dropout.

**The frailty score (ASSUMPTION throughout).** Rank each dog's hidden risk among dogs of its own entry age and turn the rank into a standard bell-curve number q (frail dogs high). The score is

    w = rho * q + sqrt(1 - rho^2) * noise

with bell-curve noise. rho is the quality of the score: 0 is useless, 1 is a perfect reading of the dog's rank. Tested: rho = 0.3, 0.5, 0.7, 0.9 and 1.0. The score is measured once, at entry, before treatment.

**True effects.**

- Sturdy dogs helped more: treated dogs in the sturdier half (hidden risk below the middle value for their entry age) have risk cut by 70 percent, the frailer half by 30 percent. No dog's rate of aging changes. This is experiment 017's case.
- Frail dogs helped more: the two cuts swapped.
- No effect; a 50 percent risk cut in every dog; aging truly slowed 20 percent in every dog.

**Analyses of every simulated trial.** The plain slope test of experiment 001 and the hidden-differences test of experiment 012, each as they were and each with two added settings: the score's effect on risk in every dog, and a "benefit differs by score" term (the treatment's risk cut may rise or fall in a straight line with the score). A test "flags" a trial when it rejects "the rate of aging did not change" at the 5 percent level; "flagged as slowed aging" means flagged with an estimated rate below the untreated one.

The sturdy-helped-more case is analysed at all five values of rho. The other cases at rho = 0.7. At least 400 simulated trials per case, fixed seeds.

## Predictions

1. **P1. Experiment 017's false finding is reproduced.** With sturdy dogs helped more and no score in the analysis, the plain test flags slowed aging in 50 to 65 percent of trials and the hidden-differences test in 82 to 93 percent.
2. **P2. A moderately good score is not enough.** With the score terms added to the plain test, slowed aging is still flagged in more than 30 percent of trials at rho = 0.5 and in 20 to 45 percent at rho = 0.7.
3. **P3. Only a nearly perfect score repairs it.** With the score terms added to the plain test, slowed aging is flagged in 8 to 25 percent of trials at rho = 0.9 and in 3 to 10 percent at rho = 1.0.
4. **P4. The score terms are cheap when they are not needed.** At rho = 0.7: with no effect, both tests with score terms flag 3 to 8 percent of trials; with aging truly slowed 20 percent, each test with score terms is no more than 5 points below the same test without them.

Reported, not scored: the hidden-differences test with score terms at every rho; the frail-helped-more and 50-percent-cut cases; how often the "benefit differs by score" term itself is found to differ from zero (which is how a trial would notice the problem); the estimated sizes.

## What would change our mind

If P2 fails because a score of quality 0.5 or 0.7 already brings false reports under 20 percent, the check that experiment 017 asked for is within reach of an ordinary measurement, and that is good news for trial planning. If P3 fails because even a perfect score leaves many false reports, the straight-line "benefit differs by score" term is the wrong repair and another analysis is needed. If P4 fails, the repair has a price even when the treatment helps every dog equally.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
