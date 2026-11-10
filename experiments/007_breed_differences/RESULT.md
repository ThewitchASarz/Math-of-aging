# Experiment 007: result

**A mix of dog sizes does not break the slope test, but ignoring size raises its false-alarm rate in short trials from 5 to about 9 percent; a separate baseline for each size class fixes that.** Whether differences between breeds are large enough to break the test is still open: the published figures give a spread anywhere from 0.15 to 0.71 depending on one unreliable breed estimate. Three predictions were confirmed and one refuted.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Spread of risk between breeds corresponds to a variance above 0.2 | 0.71, but 0.29 without one breed | Confirmed, with a strong caveat |
| P2 | Ignoring size, a 50% risk cut is misread as a slope change in over 10% of 1,300-dog, 10-year trials | 6% | **Refuted** |
| P3 | With a baseline per size class, false alarms are 3 to 8% in every case | 4 to 8% | Confirmed |
| P4 | The size-aware analysis recognises a true 20% slowing in over 80% of 10-year trials | 94% | Confirmed |

## Part 1: how different are breeds?

Converting published life expectancies into risk multipliers gives a spread (variance divided by squared mean) of:

| Group | Spread |
| --- | --- |
| 19 UK breeds, as published | 0.71 |
| Without the French Bulldog | 0.29 |
| Without the four shortest-lived breeds | 0.15 |
| Five US size classes | 0.14 |
| For comparison: frailty fitted in experiment 004 | 0.10 |
| For comparison: frailty that broke the slope test | 0.3 |

**The top figure rests on one number that should not be trusted.** The French Bulldog's published life expectancy of 4.5 years comes from a table built from deaths in 2016 to 2020, while the breed's registrations rose fourteenfold from 2011 to 2020. A fast-growing breed is mostly young dogs, so its recorded deaths are mostly young deaths. The source paper itself says the data are biased toward younger dogs and that such tables are sensitive to population change. The conversion from life expectancy to a risk multiplier is also only exact if breeds differ in level and share one rate of aging, which is not established.

## Part 2: trials with a mix of sizes

Share of 300 simulated trials in which the slope test reports a change.

| Design | Truth | Ignoring size | Size-aware |
| --- | --- | --- | --- |
| 1,300 dogs, 4 years | No effect | 9% | 7% |
| 1,300 dogs, 4 years | 30% risk cut | 6% | 5% |
| 1,300 dogs, 4 years | 20% slowing | 42% | 32% |
| 1,300 dogs, 10 years | No effect | 7% | 4% |
| 1,300 dogs, 10 years | 50% risk cut | 6% | 7% |
| 1,300 dogs, 10 years | 20% slowing | 90% | 94% |

The checker re-ran the first row with 2,000 trials: 9.2 percent ignoring size, 5.2 percent size-aware. The inflation is real, and so is the fix. Ignoring size also biases the estimated slope in short trials (a true multiplier of 0.8 is estimated at 0.70; size-aware, 0.78), which is why the "ignoring size" analysis appears to recognise more slowing in the 4-year row.

## What this means

1. **Adjust for size class.** It costs nothing in long trials and removes a real inflation in short ones.
2. **Size is not the same as breed.** Part 2 models size classes only, with a spread of 0.14. It says nothing about differences between breeds, which may be several times larger.
3. **Whether breed differences break the test is unanswered.** Experiment 004's check found that a spread of 0.3 breaks it. The breed spread is somewhere between 0.15 and 0.71. This goes to the queue.

## Independent check

A separate agent recomputed Part 1 and wrote its own simulation and both analyses for Part 2 before reading ours. Everything reproduced. It found no errors, three minor code points (one fixed here: a column was read by position, now by name), and the two findings that reshaped this result: the dependence of the breed spread on the French Bulldog figure, which it traced to the source paper, and the confirmation with 2,000 trials that the 9 percent false-alarm rate is real.

One limitation it noted remains: Part 2 uses the size-class calibration without its background death-rate term, so dogs in the simulation die somewhat less often than the calibration implies.

## Files

`breeds.py` (about 70 seconds), `results.json`, `PREDICTION.md`.
