# Experiment 003: result

**A one-year treatment is very hard to read in a three-year trial, and a benefit that fades can be mistaken for faster aging.** In a trial shaped like the rapamycin trial, slowing aging by 20 percent for the treatment year is detected in 8 percent of simulated trials, barely above the false-alarm rate. All four predictions were confirmed.

These are simulated trials under stated assumptions. They say nothing about whether rapamycin or any drug works.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | 580 dogs, 3 years: one year of 20% slowed aging, benefit kept, detected in under 10% of trials | 8% | Confirmed |
| P2 | Same design: aging fully stopped for the year, detected in 40 to 75% | 50% | Confirmed |
| P3 | A 30% risk cut that fades is flagged as a slope change in over 8% of 580-dog trials and over 25% of 5,000-dog trials, mostly as faster aging | 23% and 93%; always as faster aging | Confirmed |
| P4 | One year of 20% slowed aging is flagged as a slope change in under 10% of trials in every design | 4 to 8% | Confirmed |

## Main table

Share of 400 simulated trials. Treatment lasts one year; dogs enter at ages 7 to 12.

| What the treatment really does | 580 dogs, 3 years: detected | 580 dogs, 10 years: detected | 5,000 dogs, 3 years: detected | 580 dogs, 3 years: misread as faster aging |
| --- | --- | --- | --- | --- |
| Nothing | 6% | 4% | 4% | 3% |
| Slows aging 20% for the year, benefit kept | 8% | 9% | 23% | 3% |
| Slows aging 50% for the year, benefit kept | 17% | 43% | 86% | 1% |
| Stops aging for the year, benefit kept | 50% | 94% | 100% | 1% |
| Cuts risk 30% for the year, then fades | 10% | 10% | 51% | 23% |
| Cuts risk 50% for the year, then fades | 25% | 12% | 92% | 52% |
| Cuts risk 30%, and the cut persists | 70% | 98% | 100% | 2% |

## What this means

1. **One year of modestly slowed aging leaves almost no trace.** Slowing aging by 20 percent for a year makes a dog 0.2 years younger for good. After the treatment year that equals a lasting risk cut of about 6 percent, far below what 580 dogs can show.
2. **Even a year of completely stopped aging is a coin flip** in a 580-dog, 3-year trial. Following the same dogs for 10 years would detect it 94 percent of the time.
3. **A benefit that fades looks like harm to the slope test.** If a drug cuts risk only while it is taken, treated dogs catch up with the placebo group afterwards. The fingerprint test of experiment 001 reads that catch-up as faster aging: in 23 percent of 580-dog trials for a 30 percent cut, and in 93 percent of 5,000-dog trials.
4. **The fingerprint test is valid only while treatment continues.** Prediction 2 of the theory paper now says so. A trial with a treatment period and a follow-up period needs an analysis that models the two periods separately.

## Independent check

A separate agent wrote its own simulation and likelihood code before reading ours and reproduced every number within Monte Carlo noise. It ran our code and its own on the same simulated trials and got identical decisions on every trial. It found no errors. It noted that the 6 percent equivalence holds only after the treatment year; during the year the advantage builds up gradually.

## Assumptions

Entry ages, all-dog UK mortality, no dropout and exact doubling-law mortality, as in experiment 001. "Benefit kept" assumes a treated dog stays biologically younger by the amount gained and then ages at the normal rate. Whether any real treatment behaves that way is unknown.

## Files

`fading.py` (about 1 minute), `results.json`, `PREDICTION.md`.
