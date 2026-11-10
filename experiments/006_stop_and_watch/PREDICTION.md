# Experiment 006: prediction

Written before any code for this experiment was run.

## Question

Experiment 003 showed that when treatment stops after a year, the slope test of experiment 001 misreads a fading benefit as faster aging. Can a simple two-period analysis read such trials correctly, and what can it detect?

## Method, fixed in advance

- Same simulated trials as experiment 003: one year of treatment, then observation.
- New analysis: the treated group gets one risk ratio while on treatment and a separate one afterwards. Baseline mortality is estimated from the data as before.
- **Lasting-benefit test:** a likelihood-ratio test of "the after-treatment ratio is 1" at the 5 percent level. A benefit that survives the end of treatment is what slowed aging with a kept benefit would leave behind; a benefit that fades leaves none.
- Designs and effects as in experiment 003. At least 400 simulated trials per case.

## Predictions

1. **P1.** The new test is not fooled by a fading benefit: for risk cuts of 30 and 50 percent that fade, it reports a lasting benefit in 3 to 8 percent of trials in every design.
2. **P2.** For a 30 percent risk cut that persists, it reports a lasting benefit in at least 50 percent of 580-dog, 3-year trials.
3. **P3.** For one year of fully stopped aging with the benefit kept, it reports a lasting benefit in 30 to 60 percent of 580-dog, 3-year trials and in more than 85 percent of 580-dog, 10-year trials.
4. **P4.** For one year of 20 percent slowed aging, it reports a lasting benefit in fewer than 10 percent of 580-dog, 3-year trials.

## What would change our mind

If P1 is refuted, the two-period analysis is also fooled, and stop-and-watch trials need something more elaborate.
