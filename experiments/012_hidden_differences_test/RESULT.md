# Experiment 012: result

**Adding one setting for hidden differences between dogs repairs the slope test, and in these simulations it costs no power. When dogs differ widely in hidden risk (spread 0.3), the plain test reports a slope change in 19 percent of trials of a plain 50 percent risk cut; the new test does so in 5.5 percent. In a world with no hidden differences at all, the new test recognises a true 20 percent slowing of aging as often as the plain one (88.5 against 88.0 percent). The repair also works when the hidden differences have a different shape from the one the analysis assumes (4.5 percent).** Two predictions were confirmed, one was refuted in its main part, and one was half confirmed and half refuted.

This is a calculation under stated assumptions. It says how an analysis behaves on simulated dogs, not that any treatment works. Nobody has measured how large hidden differences between dogs really are; the spread of 0.3 is an assumption.

## What "hidden differences" means

Two dogs of the same age can carry different risks of dying for reasons no record shows. The frailest dogs die first, so the survivors look sturdier than they are, and the death rate of the group climbs more slowly at old ages than any single dog's risk does. A treatment that cuts risk keeps frail dogs alive longer, which changes this thinning-out. The plain test can mistake that for a change in the rate of aging. The new test adds one number, the spread of the hidden risk factor, and estimates it from the trial itself.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Spread 0.3, 50% risk cut: plain test reports a slope change in more than 15% of trials, new test in 3 to 8% | Plain 19% (check: 21%); new 5.5% (check: 6.1%) | Confirmed |
| P2 | New test reports a slope change in 3 to 8% of trials in all six cases with no effect or a risk cut (worlds A, B, C) | 3.5 to 6.9% | Confirmed |
| P3 | The repair costs power: with no hidden differences the plain test recognises a 20% slowing in more than 85%, the new test in 60 to 85% and at least 5 points lower; in worlds B and C the new test stays at 60% or more | Plain 88.0%, new 88.5% (check: 89.9 and 89.6). Worlds B and C: 84.5% and 81.4% | **Refuted** in its main part: no cost was found. The second part (60% or more in B and C) held |
| P4 | Different shape (log-normal, spread 0.3), 50% risk cut: plain test above 15%, new test below 10% | Plain 12.5% (check: 12.2%); new 4.5% (check: 5.0%) | **Half and half**: the plain test is disturbed less than predicted; the new test stays honest as predicted |

## Main table

Share of 800 simulated trials in which the test reports a slope change. Every trial has 1,300 dogs, entry ages 7 to 12, half treated from enrolment, and 10 years of follow-up. Each figure is uncertain by about 1 to 1.5 percentage points.

| World | True effect | Plain test | New test |
| --- | --- | --- | --- |
| A: no hidden differences | none | 6.9% | 6.9% |
| | 50% risk cut | 4.4% | 4.6% |
| | aging slowed 20% | 88.0% | 88.5% |
| B: spread 0.1 (about the size fitted to the UK table) | none | 5.8% | 5.0% |
| | 50% risk cut | 7.2% | 5.0% |
| | aging slowed 20% | 82.4% | 84.5% |
| C: spread 0.3 (ASSUMPTION) | none | 5.5% | 3.5% |
| | 50% risk cut | **18.9%** | 5.5% |
| | aging slowed 20% | 69.3% | 81.4% |
| D: spread 0.3, different shape (ASSUMPTION) | none | 4.8% | 4.6% |
| | 50% risk cut | **12.5%** | 4.5% |
| | aging slowed 20% | 75.5% | 81.0% |

In the first two rows of each world there is no slope change, so the figure is a false-alarm rate and should be near 5 percent. In the third row the slope really changes, so higher is better.

## What this means

1. **The open threat from experiments 004 and 008 has a repair.** With the extra setting, a plain risk cut is reported as a slope change about 5 percent of the time in every world tested, which is what the test is built for.
2. **The repair was free in these simulations.** We predicted it would cost at least 5 points of power. It cost none. When dogs have no hidden differences, the trial's own data say so: the estimated spread sits at zero in more than half the trials, and the new test then behaves like the plain one.
3. **The plain test's false alarms point the wrong way for a hopeful reading.** In world C with a risk cut, almost every false alarm reported *faster* aging in treated dogs (estimated slope about 10 percent steeper), not slower. So this weakness would not by itself create a false claim of slowed aging. It would make a real slowing look smaller, which is why the plain test recognises a true slowing less often in worlds C and D (69 and 76 percent) than the new test does (81 percent).
4. **The exact shape did not matter here.** In world D the analysis assumes one shape and the dogs follow another. The new test estimated a spread of about 0.16 where the true one was 0.3, and still gave 4.5 percent false alarms. This is one alternative shape, not a proof for all shapes.

For a real trial this suggests: a long trial that wants to compare slopes can include a hidden-differences setting in its analysis plan at little cost in a trial of this size and length. This is a statement about the mathematics, not a comment on any trial's own analysis plan, which answers its own questions.

## Independent check

A separate agent derived the survival formula itself, wrote its own simulation and fitting code from the description before reading ours, and ran 800 trials per case with its own seeds in six cases. All agreed with ours within simulation noise. The largest gaps were about 2 points in false-alarm rates and 3.5 points in one power figure, each under two standard errors. It then read our code. It checked the gradient numerically (error about 4 parts in a billion), confirmed the selection of frail dogs at entry in both simulators, found no failed fits in 360 trials run through our code, and confirmed that `results.json` matches every number quoted here. It confirmed the inputs: the level and slope match `dog_lifespan/results.json`, the spread of 0.1 matches experiment 004's 0.102, and experiment 008's control was 21.5 percent. It found no departure from the method fixed in `PREDICTION.md` except that 800 trials were run where at least 300 were promised.

What the check changed:

- **The power gain in worlds C and D is not described as a clean gain.** Our first reading was "the new test is more powerful when hidden differences are large". The check pointed out that the plain test is the wrong model there, so its figure is not a fair baseline. We now say only that the new test's 81 percent is trustworthy and that the plain test's lower figure reflects its bias.
- **False alarms are stated as "about 5 percent" and no finer.** With no hidden differences and no effect we drew 6.9 percent and the check drew 4.8 percent. With spread 0.3 and no effect the plain test gave 5.5 percent here and 7.4 percent in the check. These gaps are about 1.7 standard errors each; neither run can say whether either test is slightly off.
- It added limits 2, 3, 5 and 6 below.

**Not checked:** worlds B, world A with a risk cut, and world D with no effect or slowed aging were not re-run by the checker; those figures rest on our code, which the checker read and found no error in.

## How this compares with published work

That hidden differences between individuals bend a mortality curve is long published: `papers/literature-check.md` credits it to Beard (1959) and to Vaupel, Manton and Stallard (1979), read there through Missov et al. (2015); the originals were not opened. The analysis added here is that published idea applied to this project's slope test. The literature check found no hidden-differences model fitted to dogs and no measured spread for dogs (Kraus et al. 2013 discuss the idea without fitting it). So the spreads used here have no published figure to be compared with, and no published error rates for a slope test with and without such a setting were found in the searches run.

## Limits

1. **ASSUMPTION: spreads of 0.3 in worlds C and D.** No measurement exists. Only world B's 0.1 is tied to a fit (0.102 in experiment 004), and that fit is itself tentative: it comes from a deaths-only table that mixes breeds.
2. **Only two shapes were tried, both at one spread.** A wider spread, two distinct groups of dogs, or a treatment that helps frail dogs more than sturdy ones were not simulated. The new test assumes the hidden factor multiplies a dog's whole risk and is the same in both arms.
3. **This is a long trial in which 93 to 99 percent of dogs die.** That is what lets the trial estimate the spread. In a 3 or 4 year trial far fewer dogs die and the spread would be much harder to estimate. The result should not be carried over to short trials without testing.
4. **The analysis has the right shape for the simulated dogs** in worlds A to C: a doubling law, an effect constant from enrolment, no dropout, one size class. A real analysis will fit less well. Experiment 011 added dropout and sizes to the plain test; the new test has not been run with them.
5. **When there are no hidden differences the estimated spread sits on its lower limit** in over half the trials. The usual 5 percent threshold is not exact in that situation. The simulated false-alarm rates (4.6 to 6.9 percent here, 4.8 in the check) give no sign of a problem, but they cannot rule out a small one.
6. **One effect size each**: a 50 percent risk cut and a 20 percent slowing. World B barely tests the repair, because the plain test was only mildly off there (7.2 percent).
7. ASSUMPTION: entry ages 7 to 12, as in experiments 001, 010 and 011.

## Inputs

- Level and slope of dog mortality: `dog_lifespan/results.json` (fit to the UK life table of Teng et al. 2022).
- Spread 0.1: rounded from 0.102 in `experiments/004_old_dogs/results.json`.
- Simulation of gamma-shaped hidden differences: `experiments/004_old_dogs/old_dogs.py`. Plain test: `experiments/001_trial_fingerprint/fingerprint.py`.

## Files

`hidden_differences.py` (about 5 minutes), `results.json`, `PREDICTION.md`.
