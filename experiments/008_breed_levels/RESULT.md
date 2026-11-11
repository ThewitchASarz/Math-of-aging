# Experiment 008: result

**Differences in the level of risk between breeds, of the size published, do not break the slope test.** Even with the widest published spread, a plain 50 percent cut in risk was misread as a change in the rate of aging in about 5 percent of simulated trials, which is the rate the test is designed to have. Giving each breed its own baseline also works and costs no power. One prediction was refuted and three confirmed.

This is a computational experiment. It shows what the mathematics and the published breed figures imply under stated assumptions. It does not replace a real trial.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | At spread 0.15, ignoring breed, a 50% risk cut is misread in fewer than 10% of 1,300-dog, 10-year trials | 4% | Confirmed |
| P2 | Ignoring breed, the same misreading happens in more than 10% of trials at spread 0.29 and more than 15% at spread 0.71 | 5% and 6% | **Refuted** |
| P3 | With a baseline per breed, a slope change is reported in 3 to 8% of trials in every case with no true change of slope | 4 to 6% in all 12 cases | Confirmed |
| P4 | With a baseline per breed, a true 20% slowing is recognised in more than 80% of 1,300-dog, 10-year trials at every spread | 88 to 91% | Confirmed |

P2 was our main expectation and it was wrong. The reason is below.

## Main table

Share of 1,000 simulated trials in which the slope test reports a change. With 1,000 trials, each figure is uncertain by about 1 percentage point (2 points for the slowing rows).

| Breed list (spread at birth) | Design | Truth | Ignoring breed | Breed-aware |
| --- | --- | --- | --- | --- |
| All 19 groups (0.71) | 4 years | No effect | 6% | 6% |
| | 4 years | 50% risk cut | 6% | 5% |
| | 4 years | 20% slowing | 26% | 25% |
| | 10 years | No effect | 6% | 5% |
| | 10 years | 50% risk cut | 6% | 4% |
| | 10 years | 20% slowing | 88% | 91% |
| Without the French Bulldog (0.29) | 4 years | No effect | 5% | 6% |
| | 4 years | 50% risk cut | 5% | 5% |
| | 4 years | 20% slowing | 29% | 28% |
| | 10 years | No effect | 6% | 5% |
| | 10 years | 50% risk cut | 5% | 5% |
| | 10 years | 20% slowing | 87% | 89% |
| Without the four shortest-lived breeds (0.15) | 4 years | No effect | 5% | 6% |
| | 4 years | 50% risk cut | 6% | 6% |
| | 4 years | 20% slowing | 28% | 27% |
| | 10 years | No effect | 5% | 5% |
| | 10 years | 50% risk cut | 4% | 5% |
| | 10 years | 20% slowing | 87% | 88% |

All designs have 1,300 dogs. Read every "no effect" and "risk cut" figure as "about 4 to 6 percent": the differences between them are within simulation noise.

## Why the prediction failed

Experiment 004's check found that hidden differences with a spread of 0.3 push false alarms to about 20 percent. We re-ran that case here as a control and got 22 percent, so that finding stands. Yet a breed spread of 0.3 among dogs entering the trial does nothing. Two things explain it.

1. **The spread sits in a handful of rare breeds.** In the full list, five groups (French Bulldog, English Bulldog, Chihuahua, Pug, American Bulldog) carry about 82 percent of the spread but make up about 3 percent of the dogs entering a trial at ages 7 to 12. The French Bulldog alone would supply about 2 dogs out of 1,300. The other 14 groups have a spread of only about 0.06.
2. **Those dogs die early in the trial.** The slope test is disturbed when a large part of the trial population is thinned gradually over the whole follow-up. A few high-risk dogs that die in the first year or two do not bend the overall curve.

The checker measured this directly with one very large simulated trial: ignoring breed biases the estimated slope multiplier by about +0.02 (true value 1.00), against +0.10 for the control that breaks the test.

**So the single "spread" number used in experiments 004 and 007 is a poor yardstick.** What matters is how risk is spread among the bulk of the dogs, not the variance, which a few extreme breeds can dominate.

## What this means

1. **Published breed differences in level of risk are not a threat to the slope test**, with or without the doubtful French Bulldog figure. The open question from experiment 007 closes with "no" for this kind of difference.
2. **A baseline per breed is safe to use.** It keeps false alarms at about 5 percent and power at 88 to 91 percent in the 10-year design. It is not needed for level differences of this size, but it does no harm.
3. **The threat that remains is hidden differences inside breeds, and breeds that age at different rates.** Neither was simulated here. The control shows that differences spread widely across all dogs, of spread 0.3, still break the test.

## Independent check

A separate agent wrote its own simulation and both analyses before reading ours. It reproduced the spreads to four digits and every tested cell within simulation noise (for example 6.0 against our 5.5 percent, 88 against 91 percent, and 22 against 22 percent for the control). It confirmed, with a separate fit of all 22 settings, that our shortcut for the breed-aware fit (solving each breed's level exactly) gives the true best fit. It opened the source paper and matched all 19 rows of the breed table. It found no errors in the code.

What the check changed:

- It supplied the explanation above (the breakdown of the spread by breed and the large-trial bias figures), which we then adopted.
- It asked that false-alarm rates be stated as "about 4 to 6 percent", not as exact figures.
- It added limits, listed below.

## Limits

- **ASSUMPTION: breeds differ only in level of risk and share one rate of aging.** Experiment 007 suggests larger dogs age faster. A baseline per breed would not repair differences in the rate of aging.
- **No differences inside a breed.** Every dog of a breed has the same risk here. Real dogs of one breed differ, and "crossbred" is a mixed group.
- **ASSUMPTION: breed levels come from life expectancy at birth.** That includes deaths of young dogs, which do not follow the doubling law. The French Bulldog's risk at ages 7 to 12 is probably overstated. This does not matter for the conclusion, because the test is insensitive to those breeds.
- **ASSUMPTION: the breed mix at birth follows the counts in the published table**, which are the dogs in that study, not a census.
- A real trial records breed from owners, with errors. Breed is taken as known here.
- As in experiment 010: no dropout, treatment effect from the day of enrolment, and mortality that follows the doubling law at all ages.
- The control with widely spread hidden differences, and the table of breed shares at entry, were added after the main result was seen. They are not predictions.

## Files

`breed_levels.py` (about 3 minutes), `results.json`, `PREDICTION.md`.
