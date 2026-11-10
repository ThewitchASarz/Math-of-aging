# Experiment 006: result

**A two-period analysis reads stop-and-watch trials correctly: it is not fooled by a benefit that fades, and it detects a benefit that lasts.** A fading risk cut is reported as a lasting benefit in 4 to 6 percent of simulated trials, which is the false-alarm rate the test is designed for. All four predictions were confirmed. The test has clear limits, listed below.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Fading risk cuts are reported as lasting in 3 to 8% of trials in every design | 4 to 6% | Confirmed |
| P2 | A 30% cut that persists is reported as lasting in at least 50% of 580-dog, 3-year trials | 65% | Confirmed |
| P3 | A year of stopped aging, benefit kept: 30 to 60% (3 years), over 85% (10 years) | 56% and 94% | Confirmed |
| P4 | A year of 20% slowed aging: under 10% of 580-dog, 3-year trials | 8% | Confirmed |

## Main table

Share of 400 simulated trials in which the test reports a lasting benefit. Treatment lasts one year.

| What the treatment really does | 580 dogs, 3 years | 580 dogs, 10 years | 5,000 dogs, 3 years |
| --- | --- | --- | --- |
| Nothing | 5% | 7% | 4% |
| Cuts risk 30% for the year, then fades | 4% | 5% | 6% |
| Cuts risk 50% for the year, then fades | 4% | 4% | 6% |
| Cuts risk 30%, and the cut persists | 65% | 98% | 100% |
| Stops aging for the year, benefit kept | 56% | 94% | 100% |
| Slows aging 50% for the year, benefit kept | 22% | 47% | 87% |
| Slows aging 20% for the year, benefit kept | 8% | 9% | 27% |

The two-period analysis also sees a fading benefit that the usual single-ratio test mostly misses: a 30 percent cut that fades is detected as "some effect" in 32 percent of 580-dog, 3-year trials, against 10 percent for the usual test in experiment 003.

## What this means

1. **The misreading found in experiment 003 is fixed.** A trial with a treatment period and an observation period should estimate one risk ratio for each. Under that analysis a fading benefit is no longer mistaken for faster aging.
2. **A benefit that outlasts treatment is the thing to look for.** It is what slowed aging with a kept benefit would leave behind.
3. **"No lasting benefit" would not be evidence against slowed aging.** For a year of 20 percent slowed aging the test finds a lasting benefit in only 8 percent of 580-dog trials. Absence of a signal says little.

## Independent check

A separate agent wrote its own simulation and likelihood code before reading ours and reproduced every number within Monte Carlo noise. It found no errors. It added four limits, all adopted here:

- The test cannot tell slowed aging from a risk cut that simply persists. Both leave a lasting lower risk.
- Only a complete fade was tested. A partial fade would be reported as lasting.
- The designed 5 percent false-alarm rate depends on the doubling law. With strong hidden differences between dogs (frailty variance 1) the checker measured 12 percent in a 5,000-dog trial; with variance 0.5, 7 percent.
- Power for modest slowing is low at every trial size tested.

## Files

`two_period.py` (about 40 seconds), `results.json`, `PREDICTION.md`.
