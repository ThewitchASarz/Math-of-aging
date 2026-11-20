# Experiment 035: prediction

Written before any code for this experiment was written or run.

## Question

When a slope test is adjusted for a frailty score and gives entry age its own setting, does adding a setting for the hidden differences the score leaves over make false reports of slowed aging rarer or more common?

## Why it matters

Three earlier experiments lead here.

- Experiment 012 built a "hidden-differences" slope test. It has one extra setting for how much dogs differ in hidden risk. It stopped a plain cut in risk from being read as a change in the rate of aging.
- Experiment 017 found a trap. A treatment that cuts risk more in sturdy dogs than in frail ones is read as slowed aging, though no dog ages more slowly. The hidden-differences test fell into it more often than the plain test.
- Experiments 032 and 034 tried a repair: measure a frailty score at entry and include it in the analysis. Experiment 034 found the score-adjusted test must give entry age its own setting. With that done, a score of quality 0.7 still left about 40 percent of simulated trials reporting slowed aging that was not there.

Experiment 034 adjusted only the plain test. A real score will never be perfect, so dogs will still differ in ways the score does not see. A trial statistician would naturally add a setting for that left-over spread. This experiment asks whether doing so helps or harms. It also asks what each choice costs when the treatment is a plain cut in risk, the case the hidden-differences test was built for.

This experiment only measures how analyses behave on simulated dogs. It says nothing about whether any treatment works or whether any real frailty score is good enough.

## Method, fixed in advance

**The model and the trial** are those of experiments 032 and 034, unchanged. A dog that enters at age x0 has, t years later, a death risk of

    z * a * exp(b * x0) * r * exp(b * s * t)

where z is the dog's hidden personal risk factor (gamma shape, average 1 and spread 0.3 among dogs alive at age 2: ASSUMPTION), r is the treatment's cut in risk and s its change in the rate of aging. a and b are the UK fit in `dog_lifespan/results.json`. The trial has 1,300 dogs, half treated from enrolment, entry ages 7 to 12 (ASSUMPTION), 10 years, no dropout.

**Frailty scores (ASSUMPTION throughout)**, as in experiment 034: the rank score (the dog's rank in hidden risk among dogs of its own entry age, read with error) at quality 0.5, 0.7, 0.9 and 1.0, and the risk score (the dog's hidden risk itself, read with error) at quality 0.7 and 1.0. Also "no score": the analysis has no score terms at all.

**True effects (sizes are ASSUMPTIONS).** No effect. Sturdy dogs helped more: risk cut by 70 percent in the sturdier half and 30 percent in the frailer half, no change in any dog's rate of aging. Risk cut by 50 percent in every dog. Aging truly slowed 20 percent in every dog.

**Two analyses of every simulated trial.** A test "flags" a trial when it rejects "the rate of aging did not change" at the 5 percent level.

- **Entry-age** (experiment 034's): the plain slope test with the score's effect on risk, a straight-line "benefit differs by score" term, and entry age given its own setting. No hidden differences.
- **Entry-age with spread** (new here): the same, with one more setting v for hidden differences left over after the score. The left-over factor is taken to be gamma with average 1 and spread v among the dogs in the trial at entry, the same at every entry age and score. With L(t) the risk a dog has built up in the trial by time t under the plain formula, the chance of surviving to t is (1 + v * L(t)) ** (-1/v). The term that experiment 012 used for frail dogs dying before entry is dropped, because entry age's own setting takes over that job.

At least 600 simulated trials per case, fixed seeds.

## Predictions

1. **P1. The new test is sound when nothing happens.** With no effect, the entry-age-with-spread test flags 3 to 8 percent of trials with no score, with the rank score at quality 0.7 and with the rank score at quality 1.0.
2. **P2. The spread setting makes the trap worse, not better.** With sturdy dogs helped more, the entry-age-with-spread test reports slowed aging in 80 to 95 percent of trials with no score, 48 to 70 percent with the rank score at quality 0.7 (the entry-age test without spread: about 40) and 24 to 42 percent at quality 0.9 (without spread: about 21).
3. **P3. The spread setting is still needed for a plain cut in risk.** With risk cut by 50 percent in every dog, the entry-age test without spread reports faster aging in 12 to 28 percent of trials with no score and 8 to 22 percent with the rank score at quality 0.7. The entry-age-with-spread test flags 3 to 9 percent in both.
4. **P4. The spread setting costs little power.** With aging truly slowed 20 percent in every dog, the entry-age-with-spread test recognises it in 55 to 74 percent of trials with the rank score at quality 0.7, no more than 8 points below the entry-age test without spread on the same trials, and in 58 to 80 percent with no score.

Reported, not scored: every other score and quality; how often each test flags faster aging; the estimated rates; the fitted left-over spread.

## What would change our mind

If P1 fails, the reworked formula is faulty and nothing else here can be trusted until it is repaired. If P2 fails low (the spread setting brings false reports below those of the test without it), the setting is a partial repair and should be recommended with a score. If P3 fails because the test without spread is already sound, the spread setting is not needed once entry age has its own setting, and the simpler test is the better default. If P4 fails, the spread setting is too costly in a trial that learns the untreated rate only from time in the trial.

These are computational experiments. They show what the mathematics implies under stated assumptions. They do not replace a real trial, and they do not question any trial's own analysis plan, which answers its own questions.
