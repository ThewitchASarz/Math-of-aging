# Experiment 034: prediction

Written before any code for this experiment was run.

## Question

If a slope test that adjusts for a frailty score also gives entry age its own setting, does it stop leaning towards "faster aging", and how well does it then protect against a benefit that differs between dogs?

## Why it matters

Experiment 032 tested a repair for a trap found in experiment 017. The trap: a treatment that cuts risk more in sturdy dogs than in frail ones is read as slowed aging, though no dog ages more slowly. The repair: measure a frailty score at entry and let the analysis include it. Experiment 032 found the repair needs a very good score. It also found a fault in the repair itself. With a very good score and a treatment that does nothing, the adjusted test flagged 10 percent of trials, where 5 is expected, nearly all as faster aging.

The cause suggested then: the score ranks a dog among dogs of its own entry age. Once the score is known, risk rises with entry age more slowly than it rises with time in the trial. The analysis used one rate for both. The check of experiment 032 gave entry age its own setting in one very large simulated trial and the lean went away. That was one trial, not a test.

A second weakness: the "benefit differs by score" term was a straight line. The true difference in that experiment was a jump between two halves.

A trial that plans to adjust for a frailty score needs to know whether these two changes make the adjusted test safe to use. This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works or whether any real frailty score is good enough.

## Method, fixed in advance

**The model and the trial** are experiment 032's, unchanged. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r * exp(b * s * t)

where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread 0.3 among dogs alive at age 2: ASSUMPTION), r is the treatment's cut in risk and s its change in the rate of aging. a and b are the UK fit in `dog_lifespan/results.json`. The trial has 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

**Two kinds of frailty score (ASSUMPTION throughout).**

- Rank score, as in experiment 032: the dog's rank in hidden risk among dogs of its own entry age, turned into a bell-curve number q; score = rho * q + sqrt(1 - rho^2) * noise.
- Risk score: the logarithm of the dog's hidden risk itself, put on a standard scale across all dogs in the trial, mixed with noise in the same way.

rho is the quality of the score (0 useless, 1 perfect). Tested: 0.5, 0.7, 0.9 and 1.0 for the rank score; 0.7 and 1.0 for the risk score.

**True effects.** No effect. Sturdy dogs helped more: treated dogs in the sturdier half for their entry age have risk cut by 70 percent, the frailer half by 30 percent, and no dog's rate of aging changes (ASSUMPTION, as in experiments 017 and 032). Aging truly slowed 20 percent in every dog.

**Three analyses of every simulated trial.** All are the plain slope test with the score's effect on risk included. A test "flags" a trial when it rejects "the rate of aging did not change" at the 5 percent level.

- **Old**: experiment 032's. One rate for entry age and time in the trial. Benefit may differ by score in a straight line.
- **Entry-age**: the same, but the effect of entry age on risk is a separate setting from the rate at which risk rises during the trial.
- **Entry-age, five groups**: as entry-age, but the straight line is replaced by a separate risk cut for each fifth of the score (five groups of equal size).

The hidden-differences test of experiment 012 is not included: with entry age given its own setting its formula needs reworking, and that is left for later.

At least 600 simulated trials per case, fixed seeds.

## Predictions

1. **P1. The entry-age setting removes the lean.** With no effect and the rank score, the entry-age test flags 3 to 8 percent of trials at quality 0.9 and at quality 1.0 (the old test: about 10 percent at quality 1.0), and its median estimated rate is between 0.98 and 1.03 times the untreated rate at both.
2. **P2. Without the lean, the straight-line repair looks worse.** With sturdy dogs helped more and the rank score, the entry-age test reports slowed aging in 8 to 22 percent of trials at quality 1.0 (the old test: about 2 percent) and in 28 to 45 percent at quality 0.7.
3. **P3. Five groups help only when the score is very good.** With sturdy dogs helped more and the rank score, the five-group test reports slowed aging in 3 to 10 percent of trials at quality 1.0, and at quality 0.7 it is within 6 points of the entry-age test with a straight line.
4. **P4. The price in power is small.** With aging truly slowed 20 percent and the rank score: at quality 0.7 the entry-age test recognises it in 60 to 76 percent of trials; at quality 1.0 in 70 to 90 percent (the old test: about 64 percent); and at both qualities the five-group test is no more than 5 points below the entry-age test with a straight line.

Reported, not scored: every case with the risk score; qualities 0.5 and 0.9 where not named above; how often each test flags faster aging; the estimated sizes.

## What would change our mind

If P1 fails, the explanation given in experiment 032 for the lean is wrong or incomplete, and a score-adjusted slope test has a fault we do not understand. If P2 fails low (few false reports with a perfect score), the straight line was a better repair than experiment 032's extra run suggested. If P3 fails high, a benefit that differs between dogs cannot be repaired by grouping on a score, even a perfect one. If P4 fails, the entry-age setting costs too much to be the default.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
